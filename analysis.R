## ======================================================================
## 0. Setup and Arguments
## ======================================================================
options(warn=-1)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript analysis.R <p_0_value> <output_pdf_path>")
}

# Parse arguments
p0_input <- as.numeric(args[1])
output_pdf <- args[2]

cat("Running simulation with p_0 =", p0_input, "\n")
cat("Saving output to:", output_pdf, "\n")

## ======================================================================
## 1. Libraries
## ======================================================================
library(pomp)
library(deSolve)
library(tidyverse)
library(latex2exp)

## ======================================================================
## 2. Wright–Fisher POMP definition
## ======================================================================

rproc_code <- Csnippet("
  // Current allele frequency p_t = X / N
  double p_t = (N > 0) ? (X / N) : 0;

  // Selection: Δp = s * (p - 0.5) * p * (1 - p)
  double delta_p = s * (p_t - 0.5) * p_t * (1.0 - p_t);
  double p_star  = p_t + delta_p;
  if (p_star < 0) p_star = 0;
  if (p_star > 1) p_star = 1;

  // Symmetric mutation: p_next = p_star * (1 - 2μ) + μ
  double p_next = p_star * (1.0 - 2.0 * mu) + mu;
  if (p_next < 0) p_next = 0;
  if (p_next > 1) p_next = 1;

  // Wright–Fisher sampling
  X = rbinom(N_next, p_next);
")

rinit_code <- Csnippet("
  X = rbinom(N, p_0);
")

rmeas_code <- Csnippet("
  double p = (N > 0) ? (X / N) : 0;
  y = rbinom(n_samp, p);
")

dmeas_code <- Csnippet("
  double p = (N > 0) ? (X / N) : 0;
  if (ISNA(y)) {
    lik = (give_log) ? 0 : 1;
  } else {
    lik = dbinom(y, n_samp, p, give_log);
  }
")

## ======================================================================
## 3. Covariates and parameters
## ======================================================================

t_max <- 100L

true_params <- c(
  s   = 0.0,    # selection coefficient
  mu  = 1e-7,   # mutation rate
  p_0 = p0_input   # initial frequency (DYNAMIC from args)
)

time_seq <- 0:t_max

# Constant population size with a bottleneck in the middle
N0       <- 5000L
pop_size <- rep(N0, length(time_seq))
#pop_size[time_seq >= 40 & time_seq <= 60] <- 100L  # bottleneck from generation 40 to 60

covar_data <- tibble(time = time_seq, N = pop_size) %>%
  mutate(N_next = lead(N, default = last(N)))

## ======================================================================
## 4. Simulate one “true” Wright–Fisher trajectory and an observation
## ======================================================================

pomp_sim <- pomp(
  data     = tibble(time = time_seq, y = NA_integer_, n_samp = 0L),
  times    = "time",
  t0       = 0,
  rprocess = discrete_time(step.fun = rproc_code, delta.t = 1),
  rinit    = rinit_code,
  covar    = covariate_table(covar_data, times = "time"),
  statenames = "X",
  paramnames = c("s", "mu", "p_0")
)

set.seed(42)
true_sim <- simulate(pomp_sim, params = true_params, nsim = 1, format = "data.frame")

final_X <- true_sim$X[nrow(true_sim)]
final_N <- covar_data$N[covar_data$time == t_max]

sample_size <- 100L
obs_y <- rbinom(1, size = sample_size, prob = final_X / final_N)

# cat("True final frequency:", round(final_X / final_N, 3), "\n")
# cat("Observed count:", obs_y, "/", sample_size, "\n")

## ======================================================================
## 5. Full POMP object (for PF sanity check)
## ======================================================================

obs_data <- tibble(time = time_seq, y = NA_integer_, n_samp = 0L)
obs_data$y[obs_data$time == t_max]      <- obs_y
obs_data$n_samp[obs_data$time == t_max] <- sample_size

wf_pomp <- pomp(
  data     = obs_data,
  times    = "time",
  t0       = 0,
  rprocess = discrete_time(step.fun = rproc_code, delta.t = 1),
  rmeasure = rmeas_code,
  dmeasure = dmeas_code,
  rinit    = rinit_code,
  covar    = covariate_table(covar_data, times = "time"),
  statenames = "X",
  paramnames = c("s", "mu", "p_0")
)

# set.seed(10)
# pf <- pfilter(wf_pomp, params = true_params, Np = 10000)
# cat("\nPF log-likelihood at true parameters:", logLik(pf), "\n")

## ======================================================================
## 6. BDI ODE approximation (equations (5)–(6) in the note)
## ======================================================================

K_max   <- 40L               # truncate: keep k = 0..K_max
n_BDI   <- sample_size       # BDI sample size n
rho_fun <- approxfun(time_seq, pop_size / N0, rule = 2)

theta_eff <- true_params["mu"] * 4 * N0   # θ = 4 N0 μ
gamma_eff <- true_params["s"]  * 4 * N0   # γ = 4 N0 s (here h = 1)

bdi_parms <- list(
  n      = n_BDI,
  theta  = theta_eff,
  gamma  = gamma_eff,
  rho_fun = rho_fun
)

bdi_rhs <- function(t, p, parms) {
  # p[j] corresponds to p_{j-1}, j = 1..(K+1)
  K <- length(p) - 1L

  n      <- parms$n
  theta  <- parms$theta
  gamma  <- parms$gamma
  rho_t  <- parms$rho_fun(t)

  # As in the note: f = nθ/2, b(t) = n/(2ρ(t)), d(t) = n/(2ρ(t)) - γ/2
  double_f <- n * theta / 2
  double_b <- n / (2 * rho_t)
  double_d <- double_b - gamma / 2

  dp <- numeric(length(p))

  # k = 0:
  dp[1] <- -double_f * p[1] + double_d * p[2]

  # 1 <= k <= K-1
  if (K >= 2) {
    for (k in 1:(K - 1)) {
      pk_minus <- p[k]     # p_{k-1}
      pk       <- p[k + 1] # p_k
      pk_plus  <- p[k + 2] # p_{k+1}

      dp[k + 1] <- (double_f + (k - 1) * double_b) * pk_minus -
        (double_f + k * (double_b + double_d)) * pk +
        (k + 1) * double_d * pk_plus
    }
  }

  # k = K
  dp[K + 1] <- (double_f + (K - 1) * double_b) * p[K] -
    (double_f + K * (double_b + double_d)) * p[K + 1]

  list(dp)
}

## Initial distribution: Bin(n_BDI, p_0) truncated to 0..K_max
p0_full <- dbinom(0:K_max,
                  size = n_BDI,
                  prob = true_params["p_0"])
# Note: sum(p0_full) < 1 because of truncation at K_max
p0_init <- p0_full

# Integrate the ODE from t = 0 to t = t_max
t_max_ode <- t_max / (N0)
step_size <- 1 / (N0)
bdi_sol <- ode(
  y     = p0_init,
  times = seq(0, t_max_ode, by = step_size),
  func  = bdi_rhs,
  parms = bdi_parms,
  method = "lsoda"
)

## Final-time distribution at t = t_max, truncated and renormalized on 0..K_max
bdi_final_raw <- as.numeric(bdi_sol[nrow(bdi_sol), -1]) # drop time column
bdi_final_raw[bdi_final_raw < 0] <- 0                   # numerical safeguard

bdi_final_trunc <- bdi_final_raw / sum(bdi_final_raw)   # truncate 0..K and renormalize

## ======================================================================
## 7. Monte Carlo: empirical truncated distribution of Y under WF
## ======================================================================

set.seed(2025)
n_rep <- 5000L
ys <- integer(n_rep)

for (i in 1:n_rep) {
  sim_i <- simulate(pomp_sim, params = true_params, nsim = 1, format = "data.frame")
  X_T_i <- sim_i$X[nrow(sim_i)]
  ys[i] <- rbinom(1, size = sample_size, prob = X_T_i / final_N)
}

# Keep only 0..K_max and treat as a truncated sample
tab_y <- table(factor(ys, levels = 0:K_max))
emp_counts <- as.numeric(tab_y)
emp_prob_trunc <- emp_counts / sum(emp_counts)  # renormalize on 0..K_max

## ======================================================================
## 8. Truncated BDI likelihood and comparison plots
## ======================================================================

# Compare probabilities
compare_prob_df <- tibble(
  k       = 0:K_max,
  WF_emp  = emp_prob_trunc,
  BDI_ode = bdi_final_raw
)

compare_prob_long <- compare_prob_df %>%
  filter(k >= 1) %>% 
  pivot_longer(cols = c("WF_emp", "BDI_ode"),
               names_to = "source", values_to = "prob")

p_prob <- ggplot(compare_prob_long,
                 aes(x = k, y = prob, fill = source)) +
  geom_col(position = "dodge") +
  labs(
    title = TeX(paste0("WF vs BDI $ p_0 = ", p0_input, "$")),
    subtitle = TeX(paste0("Population Size $N = ", N0, "$, Sample Size $n = ", sample_size, "$, Particles $n_{rep} = ", n_rep, "$")),
    x = "Sample count k",
    y = "Probability"
  ) +
  theme_bw(base_size = 13)

# Save the plot to the path specified in arguments
ggsave(output_pdf, plot = p_prob, width = 8, height = 6)