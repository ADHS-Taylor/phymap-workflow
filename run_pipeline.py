#!/usr/bin/env python3
"""
run_pipeline.py - Standalone phymap-workflow runner (no Snakemake/conda required).

Runs the full phylogenetic transmission mapping pipeline:
  1. Prep metadata (validate FASTA <-> metadata, produce dates + phymapr metadata)
  2. Build ML tree (IQ-TREE -> RAxML-ng -> R NJ fallback)
  3. Run TreeTime (molecular clock / timetree)
  4. Run phymapr (ancestral reconstruction + animated map)
  5. Summary of outputs

Requirements:
  - Python 3.10+ with pandas, biopython installed
  - R with ape, treeio, ggtree, phymapr (or local source)
  - TreeTime (pip install phylo-treetime)
  - Optional: IQ-TREE on PATH for ML tree (best option, native Windows binary)
  - Optional: RAxML-ng on PATH (alternative ML tree)
  - If neither ML tool found, uses R ape::njs() neighbor-joining fallback

Usage:
  python run_pipeline.py --fasta aligned.fasta --metadata meta.tsv
  python run_pipeline.py --fasta data.fasta.gz --metadata meta.tsv --outdir my_results
  python run_pipeline.py --config my_config.yaml
"""

from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# =============================================================================
# Default configuration values
# =============================================================================
DEFAULTS = {
    "provider": "pathoplexus",
    "outdir": "./results",
    "raxml_model": "GTR+G",
    "clock_filter_iqd": 4,
    "treetime_coalescent": "skyline",
    "map_style": "polygon",
    "n_prune_early": 4,
    "threshold_direct_snp": 2,
    "threshold_indirect_snp": 5,
    "distant_rendering_style": "import_arrow",
    "animation_pace_balance": 0.5,
    "prune": False,
}


# =============================================================================
# Utility functions
# =============================================================================

def log_step(step: int, total: int, title: str) -> float:
    """Print a step header and return the start time."""
    print(f"\n{'='*70}")
    print(f"[step {step}/{total}] {title}")
    print(f"{'='*70}")
    return time.time()


def log_done(start: float) -> None:
    """Print elapsed time for a step."""
    elapsed = time.time() - start
    print(f"  [OK] Done in {elapsed:.1f}s")


def fatal(msg: str) -> None:
    """Print error and exit."""
    print(f"\n[ERROR] {msg}", file=sys.stderr)
    sys.exit(1)


def find_rscript() -> str | None:
    """Find Rscript executable: PATH first, then common Windows locations."""
    r = shutil.which("Rscript")
    if r:
        return r
    # Check common Windows R install paths
    program_files = Path("C:/Program Files/R")
    if program_files.exists():
        candidates = sorted(program_files.glob("*/bin/Rscript.exe"), reverse=True)
        if candidates:
            return str(candidates[0])
    program_files_x86 = Path("C:/Program Files (x86)/R")
    if program_files_x86.exists():
        candidates = sorted(program_files_x86.glob("*/bin/Rscript.exe"), reverse=True)
        if candidates:
            return str(candidates[0])
    return None


def decompress_gz(gz_path: Path, dest_dir: Path) -> Path:
    """Decompress a .gz file into dest_dir, return path to decompressed file."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    out_name = gz_path.stem  # removes .gz
    out_path = dest_dir / out_name
    print(f"  Decompressing {gz_path.name} -> {out_path}")
    with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    return out_path


def validate_alignment(fasta_path: Path) -> int:
    """Check all sequences are the same length. Returns sequence count."""
    lengths = set()
    count = 0
    current_len = 0
    in_seq = False

    with open(fasta_path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if in_seq:
                    lengths.add(current_len)
                    count += 1
                current_len = 0
                in_seq = True
            else:
                current_len += len(line)
        if in_seq:
            lengths.add(current_len)
            count += 1

    if len(lengths) == 0:
        fatal("FASTA file contains no sequences.")
    if len(lengths) > 1:
        sorted_lens = sorted(lengths)
        fatal(
            f"Sequences are NOT aligned (found {len(lengths)} different lengths: "
            f"min={sorted_lens[0]}, max={sorted_lens[-1]}). "
            f"Please provide a pre-aligned FASTA file."
        )
    print(f"  Validated: {count} sequences, all {lengths.pop()} bp (aligned OK)")
    return count


def load_yaml_config(config_path: str) -> dict:
    """Load a YAML config file. Tries PyYAML, falls back to basic parsing."""
    path = Path(config_path)
    if not path.exists():
        fatal(f"Config file not found: {config_path}")

    try:
        import yaml
        with open(path, "r") as fh:
            data = yaml.safe_load(fh)
        return data if isinstance(data, dict) else {}
    except ImportError:
        print("  [WARN] PyYAML not installed, using basic YAML parser (key: value only)")
        result = {}
        with open(path, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if ":" in line:
                    key, _, val = line.partition(":")
                    key = key.strip()
                    val = val.strip().strip('"').strip("'")
                    # Try to parse booleans and numbers
                    if val.lower() == "true":
                        result[key] = True
                    elif val.lower() == "false":
                        result[key] = False
                    else:
                        try:
                            result[key] = int(val)
                        except ValueError:
                            try:
                                result[key] = float(val)
                            except ValueError:
                                result[key] = val
        return result


def format_size(size_bytes: int) -> str:
    """Human-readable file size."""
    if size_bytes < 1024:
        return f"{size_bytes} B"
    elif size_bytes < 1024 * 1024:
        return f"{size_bytes / 1024:.1f} KB"
    else:
        return f"{size_bytes / (1024 * 1024):.1f} MB"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Standalone phymap-workflow pipeline runner. "
                    "Runs metadata prep -> tree building -> TreeTime -> phymapr "
                    "without requiring Snakemake or conda.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  python run_pipeline.py --fasta aligned.fasta --metadata meta.tsv
  python run_pipeline.py --fasta data.fasta.gz --metadata meta.tsv --prune
  python run_pipeline.py --config config.yaml --outdir custom_results
  python run_pipeline.py --fasta a.fasta --metadata m.tsv --timetree existing.nexus
""",
    )
    parser.add_argument("--fasta", type=str, help="Path to aligned FASTA file (supports .gz)")
    parser.add_argument("--metadata", type=str, help="Path to metadata TSV/CSV file (supports .gz)")
    parser.add_argument("--outdir", type=str, default=None, help="Output directory (default: ./results)")
    parser.add_argument("--provider", type=str, default=None,
                        help="Metadata provider: 'pathoplexus' or 'custom' (default: pathoplexus)")
    parser.add_argument("--prune", action="store_true", default=None,
                        help="Enable sequence pruning in prep step")
    parser.add_argument("--config", type=str, default=None,
                        help="Path to YAML config file (overrides CLI defaults)")
    parser.add_argument("--timetree", type=str, default=None,
                        help="Path to pre-existing timetree file (skips tree building + TreeTime)")
    return parser.parse_args()



# =============================================================================
# Step 1: Prep metadata
# =============================================================================

def run_prep(fasta: Path, metadata: Path, outdir: Path, cfg: dict) -> tuple[Path, Path, Path]:
    """Run prep_metadata.py by importing it directly."""
    scripts_dir = Path(__file__).resolve().parent / "scripts"
    if str(scripts_dir) not in sys.path:
        sys.path.insert(0, str(scripts_dir))

    dates_out = outdir / "dates.tsv"
    meta_out = outdir / "phymapr_metadata.tsv"
    fasta_out = outdir / "prep_aligned.fasta"

    dates_out.parent.mkdir(parents=True, exist_ok=True)

    prep_script = scripts_dir / "prep_metadata.py"
    if not prep_script.exists():
        fatal(f"prep_metadata.py not found at {prep_script}")

    # Read the script source and strip the final main() call so we control execution
    script_source = prep_script.read_text(encoding="utf-8")

    # Remove the trailing `main()` invocation so we can call it ourselves
    # after injecting the correct variable values.
    # The script ends with `main()` at module level (after function definitions).
    script_source = script_source.replace("\nmain()", "\n# main() call removed by run_pipeline.py")

    # Execute the script to define all functions and the try/except variable block.
    # The try block will fail on `snakemake.input.fasta` (NameError), and the
    # except block will set placeholder values. That's fine - we override below.
    exec_ns = {"__builtins__": __builtins__, "__name__": "__not_main__"}
    exec(compile(script_source, str(prep_script), "exec"), exec_ns)

    # Override variables with our actual values
    exec_ns["fasta_file"] = str(fasta)
    exec_ns["metadata_file"] = str(metadata)
    exec_ns["dates_out"] = str(dates_out)
    exec_ns["meta_out"] = str(meta_out)
    exec_ns["fasta_out"] = str(fasta_out)
    exec_ns["provider"] = cfg.get("provider", DEFAULTS["provider"])
    exec_ns["custom_col_specimen"] = cfg.get("custom_col_specimen", "")
    exec_ns["custom_col_date"] = cfg.get("custom_col_date", "")
    exec_ns["custom_col_location_parts"] = cfg.get("custom_col_location_parts", [])
    exec_ns["exclude_accessions"] = cfg.get("exclude_accessions", [])
    exec_ns["prune"] = cfg.get("prune", DEFAULTS["prune"])
    exec_ns["prune_snp_cutoff"] = cfg.get("prune_snp_cutoff", 5)
    exec_ns["prune_date_resolution"] = cfg.get("prune_date_resolution", "month")
    exec_ns["pruning_method"] = cfg.get("pruning_method", "fps")
    exec_ns["prune_max_representatives"] = cfg.get("prune_max_representatives", 3)
    exec_ns["prune_min_snp_difference"] = cfg.get("prune_min_snp_difference", 2)
    exec_ns["prune_clade_snp_cutoff"] = cfg.get("prune_clade_snp_cutoff", 5)

    # Now call main() once with the correct values
    exec_ns["main"]()

    # Verify outputs
    if not dates_out.exists():
        fatal("Prep step failed: dates.tsv was not created.")
    if not meta_out.exists():
        fatal("Prep step failed: phymapr_metadata.tsv was not created.")
    if not fasta_out.exists():
        fatal("Prep step failed: prep_aligned.fasta was not created.")

    return fasta_out, dates_out, meta_out


# =============================================================================
# Step 2: Build tree (RAxML-ng or R fallback)
# =============================================================================

def run_tree_building(fasta: Path, outdir: Path, cfg: dict) -> Path:
    """Build ML tree. Try iqtree2 -> raxml-ng -> R ape::njs() fallback."""
    tree_dir = outdir / "tree"
    tree_dir.mkdir(parents=True, exist_ok=True)

    model = cfg.get("raxml_model", DEFAULTS["raxml_model"])

    # --- Try IQ-TREE first (native Windows/Mac/Linux binaries available) ---
    iqtree = shutil.which("iqtree2") or shutil.which("iqtree")
    if iqtree:
        print(f"  Using IQ-TREE: {iqtree}")
        prefix = str(tree_dir / "iqtree")
        cmd = [
            iqtree,
            "-s", str(fasta),
            "-m", model,
            "--prefix", prefix,
            "-seed", str(cfg.get("raxml_seed", 42)),
            "-T", str(cfg.get("raxml_threads", "AUTO")),
            "--redo",
            "-fast",  # Fast ML search (good enough for time trees)
        ]
        print(f"  Command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        iqtree_out = tree_dir / "iqtree.treefile"
        if result.returncode == 0 and iqtree_out.exists():
            print(f"  Tree written: {iqtree_out}")
            return iqtree_out
        else:
            print(f"  [WARN] IQ-TREE failed (exit {result.returncode}):")
            if result.stderr:
                print(f"  {result.stderr[:500]}")
            print("  Trying next method...")

    # --- Try RAxML-ng ---
    raxml = shutil.which("raxml-ng")
    if raxml:
        print(f"  Using RAxML-ng: {raxml}")
        raxml_dir = tree_dir / "raxml"
        raxml_dir.mkdir(parents=True, exist_ok=True)
        prefix = str(raxml_dir / "tree")
        tree_out = raxml_dir / "tree.raxml.bestTree"
        cmd = [
            raxml, "--msa", str(fasta),
            "--model", model,
            "--prefix", prefix,
            "--seed", str(cfg.get("raxml_seed", 42)),
            "--threads", str(cfg.get("raxml_threads", "auto")),
            "--force", "model_lh_impr",
            "--redo",
        ]
        print(f"  Command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0 and tree_out.exists():
            print(f"  Tree written: {tree_out}")
            return tree_out
        else:
            print(f"  [WARN] RAxML-ng failed (exit {result.returncode}):")
            if result.stderr:
                print(f"  {result.stderr[:500]}")
            print("  Falling back to R-based tree building...")

    # Fallback: R-based tree using ape::njs()
    rscript = find_rscript()
    if not rscript:
        fatal(
            "Cannot build tree: neither raxml-ng nor Rscript found. "
            "Install RAxML-ng or R with the 'ape' package."
        )

    print(f"  Using R fallback (ape::njs): {rscript}")
    tree_out_r = tree_dir / "nj_tree.nwk"

    r_script_content = f'''\
library(ape)
cat("[R] Reading FASTA...\\n")
aln <- read.FASTA("{fasta.as_posix()}")
cat("[R] Computing distance matrix (TN93)...\\n")
dm <- dist.dna(aln, model = "TN93", pairwise.deletion = TRUE)
# Replace NA/Inf with max finite distance
dm_mat <- as.matrix(dm)
max_d <- max(dm_mat[is.finite(dm_mat)], na.rm = TRUE)
dm_mat[!is.finite(dm_mat)] <- max_d
dm <- as.dist(dm_mat)
cat("[R] Building NJ tree (njs - tolerates missing)...\\n")
tree <- njs(dm)
# Fix negative branch lengths
tree$edge.length[tree$edge.length < 0] <- 0.0
cat("[R] Writing newick tree...\\n")
write.tree(tree, file = "{tree_out_r.as_posix()}")
cat("[R] Done. Tree tips:", Ntip(tree), "\\n")
'''

    r_temp = tree_dir / "_build_tree.R"
    r_temp.write_text(r_script_content, encoding="utf-8")

    result = subprocess.run(
        [rscript, "--vanilla", str(r_temp)],
        capture_output=True, text=True,
    )
    if result.stdout:
        print(result.stdout)
    if result.returncode != 0:
        print(result.stderr)
        fatal(f"R tree building failed (exit {result.returncode})")

    # Clean up temp R script
    r_temp.unlink(missing_ok=True)

    if not tree_out_r.exists():
        fatal("R tree building completed but tree file not found.")

    print(f"  Tree written: {tree_out_r}")
    return tree_out_r


# =============================================================================
# Step 3: TreeTime
# =============================================================================

def run_treetime(tree: Path, fasta: Path, dates: Path, outdir: Path, cfg: dict) -> Path:
    """Run TreeTime to produce a time-scaled tree."""
    tt_outdir = outdir / "treetime"
    tt_outdir.mkdir(parents=True, exist_ok=True)

    clock_filter = str(cfg.get("clock_filter_iqd", DEFAULTS["clock_filter_iqd"]))
    coalescent = cfg.get("treetime_coalescent", DEFAULTS["treetime_coalescent"])

    cmd = [
        sys.executable, "-m", "treetime",
        "--tree", str(tree),
        "--aln", str(fasta),
        "--dates", str(dates),
        "--outdir", str(tt_outdir),
        "--clock-filter", clock_filter,
        "--coalescent", coalescent,
        "--confidence",
        "--verbose", "2",
    ]

    # Optional clock rate
    clock_rate = cfg.get("treetime_clock_rate", "")
    if clock_rate:
        cmd.extend(["--clock-rate", str(clock_rate)])

    print(f"  Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.stdout:
        # Print last 20 lines of stdout for summary
        lines = result.stdout.strip().split("\n")
        for line in lines[-20:]:
            print(f"  {line}")

    if result.returncode != 0:
        print(f"\n  TreeTime stderr (last 30 lines):")
        stderr_lines = result.stderr.strip().split("\n") if result.stderr else []
        for line in stderr_lines[-30:]:
            print(f"  {line}")
        fatal(f"TreeTime failed (exit {result.returncode})")

    # Find timetree output
    timetree_nexus = tt_outdir / "timetree.nexus"
    if not timetree_nexus.exists():
        # Try alternative names
        for alt_name in ["timetree.nex", "timetree.nwk", "annotated_tree.nexus"]:
            alt = tt_outdir / alt_name
            if alt.exists():
                timetree_nexus = alt
                break
        else:
            # List what's in the directory
            contents = list(tt_outdir.iterdir())
            fatal(
                f"TreeTime completed but timetree.nexus not found in {tt_outdir}. "
                f"Contents: {[f.name for f in contents]}"
            )

    print(f"  Timetree written: {timetree_nexus}")
    return timetree_nexus



# =============================================================================
# Step 4: Run phymapr (R-based ancestral reconstruction + animated map)
# =============================================================================

def run_phymapr(timetree: Path, metadata: Path, outdir: Path, cfg: dict) -> tuple[Path, Path]:
    """Run phymapr via Rscript subprocess."""
    rscript = find_rscript()
    if not rscript:
        fatal("Rscript not found. Install R to run phymapr.")

    map_out = outdir / "transmission_map.gif"
    tree_out = outdir / "tree_plot.png"

    map_style = cfg.get("map_style", DEFAULTS["map_style"])
    n_prune_early = cfg.get("n_prune_early", DEFAULTS["n_prune_early"])
    threshold_direct_snp = cfg.get("threshold_direct_snp", DEFAULTS["threshold_direct_snp"])
    threshold_indirect_snp = cfg.get("threshold_indirect_snp", DEFAULTS["threshold_indirect_snp"])
    distant_rendering_style = cfg.get("distant_rendering_style", DEFAULTS["distant_rendering_style"])
    animation_pace_balance = cfg.get("animation_pace_balance", DEFAULTS["animation_pace_balance"])

    # Check if the existing run_phymapr.R script exists and use it via CLI args
    run_phymapr_script = Path(__file__).resolve().parent / "scripts" / "run_phymapr.R"

    if run_phymapr_script.exists():
        print(f"  Using existing run_phymapr.R: {run_phymapr_script}")
        cmd = [
            rscript, "--vanilla", str(run_phymapr_script),
            "--timetree", str(timetree),
            "--metadata", str(metadata),
            "--map", str(map_out),
            "--tree_plot", str(tree_out),
            "--map_style", str(map_style),
            "--n_prune_early", str(n_prune_early),
            "--col_specimen", "accessionVersion",
            "--col_date", "sampleCollectionDate",
            "--col_location", "location",
            "--threshold_direct_snp", str(threshold_direct_snp),
            "--threshold_indirect_snp", str(threshold_indirect_snp),
            "--distant_rendering_style", str(distant_rendering_style),
            "--animation_pace_balance", str(animation_pace_balance),
        ]
    else:
        # Write a temporary R script
        print("  Writing temporary phymapr R script...")
        r_script_content = f'''\
# Temporary phymapr runner generated by run_pipeline.py

# Try loading phymapr as an installed package first
phymapr_loaded <- FALSE
tryCatch({{
  library(phymapr)
  phymapr_loaded <- TRUE
  message("[phymapr] Loaded installed phymapr package")
}}, error = function(e) {{
  message("[phymapr] Package not installed, trying local source...")
}})

if (!phymapr_loaded) {{
  # Try sourcing from companion repo
  local_path <- "C:/Users/MARTINTA/Code/phymapr/R"
  if (dir.exists(local_path)) {{
    message("[phymapr] Sourcing from: ", local_path)
    # Load required dependencies first
    required_pkgs <- c("ape", "treeio", "ggtree", "ggplot2", "dplyr",
                       "tidygeocoder", "gganimate", "lubridate", "maps",
                       "gifski", "tidyr")
    for (pkg in required_pkgs) {{
      if (requireNamespace(pkg, quietly = TRUE)) {{
        library(pkg, character.only = TRUE)
      }}
    }}
    r_files <- list.files(local_path, pattern = "\\\\.[Rr]$", full.names = TRUE)
    for (f in r_files) {{
      source(f)
    }}
    phymapr_loaded <- TRUE
  }} else {{
    stop("[phymapr] Cannot find phymapr package or local source at: ", local_path)
  }}
}}

# Run the pipeline
message("[phymapr] Running generate_phylo_transmission...")
results <- generate_phylo_transmission(
  timetree_file = "{timetree.as_posix()}",
  metadata_file = "{metadata.as_posix()}",
  col_specimen  = "accessionVersion",
  col_date      = "sampleCollectionDate",
  col_location  = "location",
  map_style     = "{map_style}",
  n_prune_early = {n_prune_early},
  threshold_direct_snp   = {threshold_direct_snp},
  threshold_indirect_snp = {threshold_indirect_snp},
  distant_rendering_style = "{distant_rendering_style}",
  animation_pace_balance  = {animation_pace_balance}
)

# Save tree plot
message("[phymapr] Saving tree plot: {tree_out.as_posix()}")
dir.create(dirname("{tree_out.as_posix()}"), showWarnings = FALSE, recursive = TRUE)
ggsave(
  filename = "{tree_out.as_posix()}",
  plot     = results$tree_plot,
  width    = 14,
  height   = 9,
  dpi      = 150
)

# Animated map
message("[phymapr] Rendering animated map: {map_out.as_posix()}")
dir.create(dirname("{map_out.as_posix()}"), showWarnings = FALSE, recursive = TRUE)
anim <- animate(
  results$map_animation,
  width  = 1200,
  height = 900,
  res    = 120,
  fps    = 10,
  renderer = gifski_renderer()
)
anim_save("{map_out.as_posix()}", animation = anim)

# Optional: HTML summary with plotly if available
tryCatch({{
  if (requireNamespace("plotly", quietly = TRUE)) {{
    message("[phymapr] Generating interactive HTML summary...")
    library(plotly)
    html_out <- "{(outdir / 'transmission_summary.html').as_posix()}"
    if (!is.null(results$tree_plot)) {{
      p <- ggplotly(results$tree_plot)
      htmlwidgets::saveWidget(p, file = html_out, selfcontained = TRUE)
      message("[phymapr] HTML summary saved: ", html_out)
    }}
  }}
}}, error = function(e) {{
  message("[phymapr] Note: plotly HTML summary skipped (", conditionMessage(e), ")")
}})

message("[phymapr] Done.")
'''
        r_temp = outdir / "_run_phymapr.R"
        r_temp.write_text(r_script_content, encoding="utf-8")
        cmd = [rscript, "--vanilla", str(r_temp)]

    print(f"  Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.stdout:
        for line in result.stdout.strip().split("\n"):
            print(f"  {line}")
    if result.stderr:
        # R prints messages to stderr
        for line in result.stderr.strip().split("\n"):
            print(f"  {line}")

    if result.returncode != 0:
        fatal(f"phymapr R script failed (exit {result.returncode})")

    # Clean up temp script if we wrote one
    temp_r = outdir / "_run_phymapr.R"
    if temp_r.exists():
        temp_r.unlink(missing_ok=True)

    # Verify outputs
    if not map_out.exists():
        fatal(f"phymapr completed but {map_out.name} was not created.")
    if not tree_out.exists():
        fatal(f"phymapr completed but {tree_out.name} was not created.")

    print(f"  Map: {map_out}")
    print(f"  Tree plot: {tree_out}")
    return map_out, tree_out


# =============================================================================
# Step 5: Summary
# =============================================================================

def print_summary(outdir: Path) -> None:
    """Print a summary table of all output files with sizes."""
    print(f"\n{'='*70}")
    print(f"  OUTPUT SUMMARY")
    print(f"{'='*70}")
    print(f"  {'File':<50} {'Size':>10}")
    print(f"  {'-'*50} {'-'*10}")

    total_size = 0
    for f in sorted(outdir.rglob("*")):
        if f.is_file() and not f.name.startswith("_"):
            size = f.stat().st_size
            total_size += size
            rel = f.relative_to(outdir)
            print(f"  {str(rel):<50} {format_size(size):>10}")

    print(f"  {'-'*50} {'-'*10}")
    print(f"  {'TOTAL':<50} {format_size(total_size):>10}")
    print(f"\n  Output directory: {outdir.resolve()}")
    print(f"{'='*70}\n")


# =============================================================================
# Main entry point
# =============================================================================

def main() -> None:
    """Main pipeline orchestrator."""
    args = parse_args()

    # ---- Load config file if provided ----
    file_cfg = {}
    if args.config:
        print(f"Loading config from: {args.config}")
        file_cfg = load_yaml_config(args.config)

    # ---- Merge config: DEFAULTS < file_cfg < CLI args ----
    cfg = dict(DEFAULTS)
    cfg.update({k: v for k, v in file_cfg.items() if v is not None})

    # CLI args override (only if explicitly provided)
    if args.provider is not None:
        cfg["provider"] = args.provider
    if args.prune is not None and args.prune:
        cfg["prune"] = True
    if args.outdir is not None:
        cfg["outdir"] = args.outdir
    elif "outdir" not in file_cfg and "results_dir" not in file_cfg:
        cfg["outdir"] = DEFAULTS["outdir"]
    elif "results_dir" in file_cfg:
        cfg["outdir"] = file_cfg["results_dir"]

    # Resolve fasta and metadata from CLI or config
    fasta_path = args.fasta or file_cfg.get("fasta")
    metadata_path = args.metadata or file_cfg.get("metadata")

    if not fasta_path:
        fatal("--fasta is required (or set 'fasta' in config YAML)")
    if not metadata_path:
        fatal("--metadata is required (or set 'metadata' in config YAML)")

    # Resolve paths relative to config file location if they're relative
    if args.config and not Path(fasta_path).is_absolute():
        config_dir = Path(args.config).resolve().parent
        fasta_candidate = config_dir / fasta_path
        if fasta_candidate.exists():
            fasta_path = str(fasta_candidate)
    if args.config and not Path(metadata_path).is_absolute():
        config_dir = Path(args.config).resolve().parent
        meta_candidate = config_dir / metadata_path
        if meta_candidate.exists():
            metadata_path = str(meta_candidate)

    fasta_path = Path(fasta_path)
    metadata_path = Path(metadata_path)
    outdir = Path(cfg["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)

    # Validate inputs exist
    if not fasta_path.exists():
        fatal(f"FASTA file not found: {fasta_path}")
    if not metadata_path.exists():
        fatal(f"Metadata file not found: {metadata_path}")

    timetree_path = None
    if args.timetree:
        timetree_path = Path(args.timetree)
        if not timetree_path.exists():
            fatal(f"Timetree file not found: {timetree_path}")

    # ---- Pipeline banner ----
    total_steps = 5 if timetree_path is None else 4
    print(f"\n{'='*70}")
    print(f"  PHYMAP-WORKFLOW STANDALONE PIPELINE")
    print(f"{'='*70}")
    print(f"  FASTA:    {fasta_path}")
    print(f"  Metadata: {metadata_path}")
    print(f"  Output:   {outdir.resolve()}")
    print(f"  Provider: {cfg['provider']}")
    print(f"  Prune:    {cfg['prune']}")
    if timetree_path:
        print(f"  Timetree: {timetree_path} (pre-existing, skipping build+TreeTime)")
    print(f"{'='*70}")

    pipeline_start = time.time()

    # ---- Step 0: Decompress .gz inputs if needed ----
    temp_dir = outdir / "_temp"
    if fasta_path.suffix == ".gz":
        temp_dir.mkdir(parents=True, exist_ok=True)
        print(f"\n  Auto-decompressing FASTA (.gz detected)...")
        fasta_path = decompress_gz(fasta_path, temp_dir)
    if metadata_path.suffix == ".gz":
        temp_dir.mkdir(parents=True, exist_ok=True)
        print(f"  Auto-decompressing metadata (.gz detected)...")
        metadata_path = decompress_gz(metadata_path, temp_dir)

    # ---- Step 1: Validate alignment ----
    step_num = 1
    t = log_step(step_num, total_steps, "Validate alignment & prep metadata")
    validate_alignment(fasta_path)

    try:
        prep_fasta, dates_file, phymapr_meta = run_prep(
            fasta_path, metadata_path, outdir, cfg
        )
    except SystemExit as e:
        fatal(f"Prep step failed: {e}")
    except Exception as e:
        fatal(f"Prep step failed with exception: {e}")
    log_done(t)

    # ---- Step 2: Build tree (or use provided timetree) ----
    step_num += 1
    if timetree_path is None:
        t = log_step(step_num, total_steps, "Build phylogenetic tree")
        try:
            tree_file = run_tree_building(prep_fasta, outdir, cfg)
        except Exception as e:
            fatal(f"Tree building failed: {e}")
        log_done(t)

        # ---- Step 3: TreeTime ----
        step_num += 1
        t = log_step(step_num, total_steps, "Run TreeTime (molecular clock)")
        try:
            timetree_path = run_treetime(tree_file, prep_fasta, dates_file, outdir, cfg)
        except Exception as e:
            fatal(f"TreeTime failed: {e}")
        log_done(t)
    else:
        t = log_step(step_num, total_steps, "Using pre-existing timetree (skipping build + TreeTime)")
        print(f"  Timetree: {timetree_path}")
        log_done(t)

    # ---- Step 4: phymapr ----
    step_num += 1
    t = log_step(step_num, total_steps, "Run phymapr (ancestral reconstruction + map)")
    try:
        map_file, tree_plot = run_phymapr(timetree_path, phymapr_meta, outdir, cfg)
    except Exception as e:
        fatal(f"phymapr failed: {e}")
    log_done(t)

    # ---- Step 5: Summary ----
    step_num += 1
    t = log_step(step_num, total_steps, "Output summary")
    print_summary(outdir)

    # Clean up temp directory
    if temp_dir.exists():
        try:
            shutil.rmtree(temp_dir)
        except Exception:
            pass

    elapsed_total = time.time() - pipeline_start
    print(f"  Pipeline completed in {elapsed_total:.1f}s")
    print(f"  All outputs in: {outdir.resolve()}\n")


if __name__ == "__main__":
    main()
