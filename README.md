# phymap-workflow

A standalone pipeline that takes aligned FASTA sequences and sample metadata through to an animated phylogeographic transmission map. Handles tree building, molecular clock fitting, and visualization in a single command.

```
aligned FASTA + metadata
        |
        v
   [1] prep       validate IDs, geocode locations, optional pruning
        |
        v
   [2] tree       IQ-TREE or RAxML-ng (ML), or R neighbor-joining fallback
        |
        v
   [3] treetime   molecular clock, time-scaled annotated NEXUS
        |
        v
   [4] phymapr    ancestral reconstruction, transmission map + tree plot
```

## Requirements

**Python 3.10+**
```bash
pip install pandas biopython phylo-treetime numpy
```

**R 4.x+**
```R
install.packages(c("BiocManager", "remotes"))
BiocManager::install(c("ggtree", "treeio"))
remotes::install_github("ADHS-Taylor/phymapr")
```

**Optional (recommended): IQ-TREE for ML tree building**

Download the binary for your platform from [iqtree.github.io](https://iqtree.github.io/) and add it to your PATH. If IQ-TREE isn't found, the pipeline falls back to neighbor-joining in R (works fine for smaller datasets).

No conda, Snakemake, or WSL required.

## Quick Start

### Pathoplexus data

Download an **aligned nucleotide FASTA** and **metadata TSV** from [pathoplexus.org](https://pathoplexus.org) for any supported pathogen, then:

```bash
python run_pipeline.py --fasta rsv-b_aligned-nuc.fasta.gz --metadata rsv-b_metadata.tsv.gz --outdir results/rsv-b --prune
```

The pipeline handles `.gz` decompression, validates alignment, and runs all four stages. Results land in `results/rsv-b/`:

| File | Description |
|------|-------------|
| `tree_plot.png` | Time-scaled phylogenetic tree |
| `transmission_map_static.png` | All pathways at once (easiest to read) |
| `transmission_map.gif` | Animated version |
| `treetime/timetree.nexus` | Annotated NEXUS for downstream use |

### Custom data

If your data isn't from Pathoplexus, use `--provider custom` and a config file:

```bash
python run_pipeline.py --fasta my_alignment.fasta --metadata my_metadata.csv --provider custom --config config/my_config.yaml --outdir results/my_run
```

Your config specifies column mappings:
```yaml
provider: "custom"
custom_col_specimen: "sample_id"
custom_col_date: "collection_date"
custom_col_location_parts:
  - "city"
  - "state"
  - "country"
```

Specimen IDs must match FASTA headers exactly. Dates should be ISO format (YYYY-MM-DD or YYYY-MM). Location columns get joined into a string and geocoded via OpenStreetMap.

### Already have a timetree?

Skip tree building entirely:

```bash
python run_pipeline.py --fasta aligned.fasta --metadata meta.tsv --timetree my_tree.nexus --outdir results
```

## Simulated Examples

The `simulated_data/` folder contains three test datasets (Virus X, Y, Z) that run the full pipeline without needing external downloads. These are described in detail in the [phymapr](https://github.com/ADHS-Taylor/phymapr) README.

```bash
python run_pipeline.py --fasta simulated_data/virus_x/virus_x_sequences.fasta --metadata simulated_data/virus_x/virus_x_metadata.csv --provider custom --outdir results/virus_x
```

## Options

```
python run_pipeline.py --help

  --fasta FASTA        Aligned FASTA file (supports .gz)
  --metadata METADATA  Metadata TSV/CSV (supports .gz)
  --outdir OUTDIR      Output directory (default: ./results)
  --provider PROVIDER  "pathoplexus" or "custom" (default: pathoplexus)
  --prune              Enable sequence pruning (reduces redundancy)
  --config CONFIG      YAML config file (overrides CLI defaults)
  --timetree TIMETREE  Pre-existing timetree (skips tree building + TreeTime)
```

## Pruning

When `--prune` is enabled, the prep step computes pairwise SNP distances and clusters sequences using complete-linkage. One representative per location/month/clade is kept. This reduces dataset size substantially (often 50-80% reduction) while preserving geographic and genetic diversity.

SNP distances are computed in Python (no external tools required). If `pairsnp` is installed (`pip install pairsnp`) it's used for speed; otherwise falls back to numpy.

## Snakemake (advanced)

The `Snakefile` is still included for users who prefer Snakemake orchestration with conda environments. See the config examples in `config/`. This requires conda and is not necessary for most use cases — `run_pipeline.py` does the same thing without the dependency overhead.

## Notes

- Sequences must be **aligned** (all same length). If you download from Pathoplexus, choose "Aligned nucleotide" not just "Nucleotide."
- For slow-evolving pathogens like mpox, TreeTime may fail to fit a clock automatically. Set `treetime_clock_rate` in your config (e.g., `6e-5` for mpox).
- Datasets with multiple divergent lineages (e.g., measles D8 + B3) should be split by genotype before running. Mixing lineages that diverged centuries ago produces meaningless deep-root trees.
- Location geocoding requires internet on first run for each unique location string.

## Project Structure

```
phymap-workflow/
├── run_pipeline.py          # main entry point (no Snakemake needed)
├── Snakefile                # alternative Snakemake orchestration
├── config/
│   ├── example_config.yaml  # annotated config template
│   └── virus_*_config.yaml  # configs for simulated data
├── envs/                    # conda env definitions (for Snakemake path)
├── scripts/
│   ├── prep_metadata.py     # metadata validation + pruning
│   └── run_phymapr.R        # R wrapper for phymapr
├── simulated_data/          # Virus X, Y, Z test datasets
│   ├── virus_x/
│   ├── virus_y/
│   └── virus_z/
└── data/                    # place your input files here
```
