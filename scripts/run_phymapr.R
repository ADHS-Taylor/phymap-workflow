# =============================================================================
# run_phymapr.R
# Stage 4 of phymap-workflow: ancestral state reconstruction + animated map
#
# Called by Snakemake via the `script:` directive.
# Accesses inputs/outputs/params through the `snakemake` S4 object.
# =============================================================================

# ---------------------------------------------------------------------------
# 1. Ensure phymapr is installed
#    All heavy dependencies (ggtree, treeio, phytools, etc.) are provided by
#    the conda environment (envs/r_phymapr.yaml).  phymapr itself is a
#    lightweight R-only package installed once from GitHub.
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# 1. Load dependencies and source phymapr directly to bypass compiler issues
# ---------------------------------------------------------------------------
library(ape)
library(treeio)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tidygeocoder)
library(gganimate)
library(lubridate)
library(maps)
library(gifski)
library(tidyr)

if (dir.exists("C:/Code/phymapr")) {
  message("[phymapr] Loading local phymapr package by sourcing R files directly (compiler-free)...")
  r_files <- list.files("C:/Code/phymapr/R", pattern = "\\.[Rr]$", full.names = TRUE)
  for (f in r_files) {
    source(f)
  }
} else {
  if (!requireNamespace("phymapr", quietly = TRUE)) {
    message("[phymapr] phymapr not found – installing from GitHub...")
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = "https://cloud.r-project.org")
    }
    remotes::install_github("ADHS-Taylor/phymapr", upgrade = "never")
  }
  library(phymapr)
}

# ---------------------------------------------------------------------------
# 2. Read parameters from Command Line Arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# Default values
timetree_file <- NULL
metadata_file <- NULL
map_out <- NULL
tree_out <- NULL
map_style <- "polygon"
n_prune_early <- 4
col_specimen <- "accessionVersion"
col_date <- "sampleCollectionDate"
col_location <- "location"
threshold_direct_snp <- 2
threshold_indirect_snp <- 5
distant_rendering_style <- "import_arrow"
animation_pace_balance <- 0.5

# Parse arguments
i <- 1
while (i <= length(args)) {
  arg <- args[i]
  if (arg == "--timetree") {
    timetree_file <- args[i+1]; i <- i + 2
  } else if (arg == "--metadata") {
    metadata_file <- args[i+1]; i <- i + 2
  } else if (arg == "--map") {
    map_out <- args[i+1]; i <- i + 2
  } else if (arg == "--tree_plot") {
    tree_out <- args[i+1]; i <- i + 2
  } else if (arg == "--map_style") {
    map_style <- args[i+1]; i <- i + 2
  } else if (arg == "--n_prune_early") {
    n_prune_early <- as.integer(args[i+1]); i <- i + 2
  } else if (arg == "--col_specimen") {
    col_specimen <- args[i+1]; i <- i + 2
  } else if (arg == "--col_date") {
    col_date <- args[i+1]; i <- i + 2
  } else if (arg == "--col_location") {
    col_location <- args[i+1]; i <- i + 2
  } else if (arg == "--threshold_direct_snp") {
    threshold_direct_snp <- as.integer(args[i+1]); i <- i + 2
  } else if (arg == "--threshold_indirect_snp") {
    threshold_indirect_snp <- as.integer(args[i+1]); i <- i + 2
  } else if (arg == "--distant_rendering_style") {
    distant_rendering_style <- args[i+1]; i <- i + 2
  } else if (arg == "--animation_pace_balance") {
    animation_pace_balance <- as.numeric(args[i+1]); i <- i + 2
  } else {
    i <- i + 1
  }
}

message("[phymapr] timetree : ", timetree_file)
message("[phymapr] metadata : ", metadata_file)
message("[phymapr] map_style: ", map_style)
message("[phymapr] cols     : specimen=", col_specimen,
        "  date=", col_date, "  location=", col_location)
message("[phymapr] advanced : threshold_direct_snp=", threshold_direct_snp,
        "  threshold_indirect_snp=", threshold_indirect_snp,
        "  distant_rendering_style=", distant_rendering_style,
        "  animation_pace_balance=", animation_pace_balance)

# ---------------------------------------------------------------------------
# 3. Ensure output directories exist
# ---------------------------------------------------------------------------
dir.create(dirname(map_out),  showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(tree_out), showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 4. Run the full phymapr pipeline
# ---------------------------------------------------------------------------
results <- generate_phylo_transmission(
  timetree_file = timetree_file,
  metadata_file = metadata_file,
  col_specimen  = col_specimen,
  col_date      = col_date,
  col_location  = col_location,
  map_style     = map_style,
  n_prune_early = n_prune_early,
  threshold_direct_snp   = threshold_direct_snp,
  threshold_indirect_snp = threshold_indirect_snp,
  distant_rendering_style = distant_rendering_style,
  animation_pace_balance  = animation_pace_balance
)

# ---------------------------------------------------------------------------
# 5. Save outputs
# ---------------------------------------------------------------------------

# Tree plot  →  PNG
message("[phymapr] Saving tree plot: ", tree_out)
ggsave(
  filename = tree_out,
  plot     = results$tree_plot,
  width    = 14,
  height   = 9,
  dpi      = 150
)

# Animated map  →  GIF
message("[phymapr] Rendering and saving animated map: ", map_out)
anim <- animate(
  results$map_animation,
  width  = 1200,
  height = 900,
  res    = 120,
  fps    = 10,
  renderer = gifski_renderer()
)
anim_save(map_out, animation = anim)

message("[phymapr] Done.")
message("  Tree : ", tree_out)
message("  Map  : ", map_out)
