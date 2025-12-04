# List of p_0 values to iterate over
P0_VALUES = ["1e-3", "5e-3", "1e-2", "5e-2"]

# The default rule (first rule) defines what we want to produce
rule all:
    input:
        expand("results/p_prob_p0_{p0}.pdf", p0=P0_VALUES)

# Rule to run the analysis script
rule run_analysis:
    input:
        script = "analysis.R"
    output:
        plot = "results/p_prob_p0_{p0}.pdf"  # Fixed: use {p0}, not {wildcards.p0}
    params:
        p0 = "{p0}"
    shell:
        """
        Rscript {input.script} {params.p0} {output.plot}
        """