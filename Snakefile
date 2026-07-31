# =============================================================================
# phymap-workflow  |  Snakemake phylogenetic transmission pipeline
#
# Stages:
#   1. prep      – parse FASTA headers + metadata → dates file + phymapr metadata
#   2. raxml     – RAxML-ng maximum-likelihood tree
#   3. treetime  – TreeTime molecular-clock tree (annotated NEXUS)
#   4. phymapr   – ancestral state reconstruction → animated transmission map
#
# Usage:
#   snakemake --use-conda -j <threads>
# =============================================================================

configfile: "config/config.yaml"

RESULTS_BASE = config.get("results_dir", "results")
if config.get("prune", False):
    RESULTS = f"{RESULTS_BASE}_pruned"
else:
    RESULTS = f"{RESULTS_BASE}_non_pruned"

import sys
import shutil

ON_WINDOWS = sys.platform.startswith("win")

def get_raxml_command():
    if shutil.which("raxml-ng"):
        return "raxml-ng"
    elif ON_WINDOWS:
        return "wsl ~/.local/bin/raxml-ng"
    else:
        return "raxml-ng"



# -----------------------------------------------------------------------------
# Target rule – declare final outputs
# -----------------------------------------------------------------------------
rule all:
    input:
        f"{RESULTS}/transmission_map.gif",
        f"{RESULTS}/tree_plot.png",


# -----------------------------------------------------------------------------
# Stage 1: Metadata preparation
#   - Validates FASTA header ↔ metadata specimen column alignment
#   - Builds TreeTime dates file  (name, date)
#   - Builds phymapr metadata TSV (specimen, date, location, extras)
# -----------------------------------------------------------------------------
rule prep:
    input:
        fasta    = config["fasta"],
        metadata = config["metadata"],
    output:
        dates        = f"{RESULTS}/dates.tsv",
        phymapr_meta = f"{RESULTS}/phymapr_metadata.tsv",
        fasta_out    = f"{RESULTS}/prep_aligned.fasta",
    params:
        provider                  = config.get("provider", "pathoplexus"),
        custom_col_specimen       = config.get("custom_col_specimen", ""),
        custom_col_date           = config.get("custom_col_date", ""),
        custom_col_location_parts = config.get("custom_col_location_parts", []),
        exclude_accessions        = config.get("exclude_accessions", []),
        prune                     = config.get("prune", False),
        prune_snp_cutoff          = config.get("prune_snp_cutoff", 5),
        prune_date_resolution     = config.get("prune_date_resolution", "month"),
        pruning_method            = config.get("pruning_method", "fps"),
        prune_max_representatives = config.get("prune_max_representatives", 3),
        prune_min_snp_difference  = config.get("prune_min_snp_difference", 2),
        prune_clade_snp_cutoff    = config.get("prune_clade_snp_cutoff", 5),
    conda:
        "envs/prep.yaml"
    script:
        "scripts/prep_metadata.py"


# -----------------------------------------------------------------------------
# Stage 2: Maximum-likelihood tree (RAxML-ng)
# -----------------------------------------------------------------------------
rule raxml:
    input:
        fasta = f"{RESULTS}/prep_aligned.fasta",
    output:
        tree = f"{RESULTS}/raxml/tree.raxml.bestTree",
    params:
        model     = config.get("raxml_model", "GTR+G"),
        prefix    = f"{RESULTS}/raxml/tree",
        seed      = config.get("raxml_seed", 42),
        threads   = config.get("raxml_threads", "auto"),
        raxml_bin = get_raxml_command()
    shell:
        "{params.raxml_bin} --msa {input.fasta} --model {params.model} --prefix {params.prefix} --seed {params.seed} --threads {params.threads} --force model_lh_impr --redo"


# -----------------------------------------------------------------------------
# Stage 3: Molecular-clock time tree (TreeTime)
#   Produces an annotated NEXUS readable by treeio::read.beast()
# -----------------------------------------------------------------------------
rule treetime:
    input:
        tree  = f"{RESULTS}/raxml/tree.raxml.bestTree",
        fasta = f"{RESULTS}/prep_aligned.fasta",
        dates = f"{RESULTS}/dates.tsv",
    output:
        nexus = f"{RESULTS}/treetime/timetree.nexus",
    params:
        outdir         = f"{RESULTS}/treetime",
        clock_filter   = config.get("clock_filter_iqd", 4),
        coalescent     = config.get("treetime_coalescent", "skyline"),
        clock_rate_opt = lambda wildcards: f"--clock-rate {config['treetime_clock_rate']}" if config.get("treetime_clock_rate") else ""
    conda:
        "envs/treetime.yaml"
    shell:
        "treetime --tree {input.tree} --aln {input.fasta} --dates {input.dates} --outdir {params.outdir} --clock-filter {params.clock_filter} --coalescent {params.coalescent} --confidence --verbose 2 {params.clock_rate_opt}"


# -----------------------------------------------------------------------------
# Stage 4: phymapr  –  ancestral state reconstruction + animated map
# -----------------------------------------------------------------------------
rule phymapr:
    input:
        timetree = f"{RESULTS}/treetime/timetree.nexus",
        metadata = f"{RESULTS}/phymapr_metadata.tsv",
    output:
        map        = f"{RESULTS}/transmission_map.gif",
        tree_plot  = f"{RESULTS}/tree_plot.png",
    params:
        map_style     = config.get("map_style", "polygon"),
        n_prune_early = config.get("n_prune_early", 4),
        # Column names in the phymapr_metadata.tsv produced by the prep step.
        # These are set automatically by the prep script; only change if you
        # modified prep_metadata.py to use different output column names.
        col_specimen  = "accessionVersion",
        col_date      = "sampleCollectionDate",
        col_location  = "location",
        threshold_direct_snp    = config.get("transmission_threshold_direct_snp", 2),
        threshold_indirect_snp  = config.get("transmission_threshold_indirect_snp", 5),
        distant_rendering_style = config.get("distant_rendering_style", "import_arrow"),
        animation_pace_balance  = config.get("animation_pace_balance", 0.5),
    shell:
        "conda run -p .snakemake/conda/d23c926e285bb4f5325d40134a9268ea_ Rscript --vanilla scripts/run_phymapr.R "
        "--timetree {input.timetree} "
        "--metadata {input.metadata} "
        "--map {output.map} "
        "--tree_plot {output.tree_plot} "
        "--map_style {params.map_style} "
        "--n_prune_early {params.n_prune_early} "
        "--col_specimen {params.col_specimen} "
        "--col_date {params.col_date} "
        "--col_location {params.col_location} "
        "--threshold_direct_snp {params.threshold_direct_snp} "
        "--threshold_indirect_snp {params.threshold_indirect_snp} "
        "--distant_rendering_style {params.distant_rendering_style} "
        "--animation_pace_balance {params.animation_pace_balance}"
