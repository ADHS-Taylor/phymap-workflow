# phymap-workflow

A reproducible Snakemake pipeline that takes an aligned FASTA and sample metadata from end to end: maximum-likelihood tree → molecular-clock time tree → animated phylogeographic transmission map.

```
aligned FASTA + metadata
        │
        ▼
   [1] prep        parse headers, build dates file & phymapr metadata
        │
        ▼
   [2] RAxML-ng    maximum-likelihood tree (GTR+G)
        │
        ▼
   [3] TreeTime    molecular-clock time-scaled tree (annotated NEXUS)
        │
        ▼
   [4] phymapr     ancestral state reconstruction → tree plot + animated map GIF
```

---

## Installation (2 commands)

You need [conda / miniforge](https://github.com/conda-forge/miniforge) once on your system.  Everything else is handled automatically.

```bash
# 1. Install Snakemake into your base environment
conda install -n base -c conda-forge snakemake

# 2. Run the pipeline (Snakemake creates all tool environments on first run)
snakemake --use-conda -j 8
```

That's it.  Each stage runs in its own isolated conda environment — no manual tool installation needed.

---

## Quick start: Pathoplexus data

Pathoplexus exports pre-aligned FASTA + TSV metadata with consistent column names, making it the fastest path to a result.

1. Download an **aligned nucleotide FASTA** and its **metadata TSV** from [Pathoplexus](https://pathoplexus.org).
2. Place both files in the `data/` folder of this repo.
3. Edit `config/config.yaml` — change only the two `fasta:` and `metadata:` paths:

```yaml
fasta:    "data/your_virus_aligned.fasta"
metadata: "data/your_virus_metadata.tsv"
provider: "pathoplexus"   # leave as-is
```

4. Run:

```bash
snakemake --use-conda -j 8
```

Results land in `results/`:

| File | Description |
|------|-------------|
| `results/transmission_map.gif` | Animated geographic transmission map |
| `results/tree_plot.png` | Time-scaled phylogenetic tree |
| `results/treetime/timetree.nexus` | Annotated NEXUS tree (BEAST-compatible) |
| `results/raxml/tree.raxml.bestTree` | ML tree (newick) |
| `results/dates.tsv` | Dates file used by TreeTime |
| `results/phymapr_metadata.tsv` | Cleaned metadata passed to phymapr |

---

## Using your own data (custom provider)

If your metadata doesn't come from Pathoplexus, set `provider: "custom"` and specify your column names:

```yaml
fasta:    "data/my_sequences_aligned.fasta"
metadata: "data/my_metadata.csv"

provider: "custom"
custom_col_specimen: "sample_id"          # must match FASTA header IDs exactly
custom_col_date:     "collection_date"    # ISO 8601: YYYY-MM-DD (or YYYY-MM, YYYY)
custom_col_location_parts:
  - "city"
  - "state"
  - "country"
```

Location columns are joined into a single geocodable string (e.g. `"Phoenix, Arizona, USA"`) and passed to [tidygeocoder](https://jessecambon.github.io/tidygeocoder/) (OpenStreetMap) to resolve coordinates.  An internet connection is required for this step.

---

## Configuration reference

All options live in `config/config.yaml`.

| Key | Default | Description |
|-----|---------|-------------|
| `fasta` | — | Path to aligned FASTA |
| `metadata` | — | Path to metadata TSV/CSV |
| `provider` | `pathoplexus` | `pathoplexus` or `custom` |
| `raxml_model` | `GTR+G` | RAxML-ng substitution model |
| `raxml_seed` | `42` | Random seed for reproducibility |
| `raxml_threads` | `auto` | Thread count for RAxML-ng |
| `clock_filter_iqd` | `4` | TreeTime outlier filter (IQDs from clock regression) |
| `treetime_coalescent` | `skyline` | TreeTime coalescent model: `skyline`, `opt`, or a numeric Ne |
| `treetime_clock_rate` | `""` | Fixed clock rate (leave blank to auto-infer) |
| `map_style` | `polygon` | `polygon` (fast/light) or `tile` (dark/satellite) |
| `n_prune_early` | `4` | Number of early ancestor nodes to prune from animation |
| `results_dir` | `results` | Output directory |

---

## Supported Pathoplexus metadata columns

The prep step automatically maps these Pathoplexus columns:

| Role | Column |
|------|--------|
| Specimen ID (= FASTA header) | `accessionVersion` |
| Collection date | `sampleCollectionDate` |
| Location (geocoded) | `geoLocCity` + `geoLocAdmin1` + `geoLocCountry` |

Additional columns (`genotype`, `geoLocCountry`, etc.) are carried forward into the phymapr metadata TSV.

---

## Re-running a specific stage

Snakemake skips stages whose outputs already exist.  To force a specific rule to re-run:

```bash
snakemake --use-conda -j 8 --forcerun raxml
```

To re-run everything from scratch:

```bash
snakemake --use-conda -j 8 --forceall
```

---

## Dependencies

Each stage runs in its own auto-managed conda environment:

| Stage | Environment | Key tools |
|-------|-------------|-----------|
| prep | `envs/prep.yaml` | Python, pandas |
| raxml | `envs/raxml.yaml` | RAxML-ng ≥ 1.2 |
| treetime | `envs/treetime.yaml` | TreeTime ≥ 0.11 |
| phymapr | `envs/r_phymapr.yaml` | R ≥ 4.3, ggtree, treeio, phytools, phymapr |

[phymapr](https://github.com/ADHS-Taylor/phymapr) is installed automatically from GitHub the first time the R environment is used.

---

## Project structure

```
phymap-workflow/
├── Snakefile               # pipeline definition
├── config/
│   └── config.yaml         # all user-configurable settings
├── envs/
│   ├── prep.yaml           # Python environment for metadata prep
│   ├── raxml.yaml          # RAxML-ng environment
│   ├── treetime.yaml       # TreeTime environment
│   └── r_phymapr.yaml      # R + phymapr environment
├── scripts/
│   ├── prep_metadata.py    # Stage 1: FASTA/metadata parser
│   └── run_phymapr.R       # Stage 4: phymapr wrapper
├── data/                   # place your input files here (gitignored)
└── results/                # pipeline outputs (gitignored)
```
