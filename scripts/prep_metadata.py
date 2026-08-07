"""
prep_metadata.py
================
Stage 1 of phymap-workflow.

Reads the aligned FASTA and metadata file, validates that sequence IDs match,
then produces two output files:

  dates.tsv           – (name, date) for TreeTime
  phymapr_metadata.tsv – (accessionVersion, sampleCollectionDate, location, …)
                          for phymapr's generate_phylo_transmission()

Provider-aware:
  provider="pathoplexus"  →  uses known Pathoplexus column names
  provider="custom"       →  uses columns from config (custom_col_* params)
"""

import sys
import pandas as pd
from pathlib import Path
import datetime
import calendar


# ---------------------------------------------------------------------------
# Snakemake injects a `snakemake` object when running via `script:` directive.
# For interactive testing, fall back to dummy values.
# ---------------------------------------------------------------------------
try:
    fasta_file    = snakemake.input.fasta          # noqa: F821
    metadata_file = snakemake.input.metadata       # noqa: F821
    dates_out     = snakemake.output.dates         # noqa: F821
    meta_out      = snakemake.output.phymapr_meta  # noqa: F821
    fasta_out     = snakemake.output.fasta_out     # noqa: F821

    provider                  = snakemake.params.provider                   # noqa: F821
    custom_col_specimen       = snakemake.params.custom_col_specimen        # noqa: F821
    custom_col_date           = snakemake.params.custom_col_date            # noqa: F821
    custom_col_location_parts = snakemake.params.custom_col_location_parts  # noqa: F821
    exclude_accessions        = snakemake.params.exclude_accessions         # noqa: F821
    prune                     = snakemake.params.prune                      # noqa: F821
    prune_snp_cutoff          = snakemake.params.prune_snp_cutoff           # noqa: F821
    prune_date_resolution     = snakemake.params.prune_date_resolution      # noqa: F821
    pruning_method            = snakemake.params.pruning_method             # noqa: F821
    prune_max_representatives = snakemake.params.prune_max_representatives  # noqa: F821
    prune_min_snp_difference  = snakemake.params.prune_min_snp_difference   # noqa: F821
    prune_clade_snp_cutoff    = snakemake.params.prune_clade_snp_cutoff     # noqa: F821
except NameError:
    print("Running outside Snakemake – using placeholder values for testing.")
    fasta_file    = "data/test.fasta"
    metadata_file = "data/test_metadata.tsv"
    dates_out     = "results/dates.tsv"
    meta_out      = "results/phymapr_metadata.tsv"
    fasta_out     = "results/prep_aligned.fasta"
    provider      = "pathoplexus"
    custom_col_specimen       = ""
    custom_col_date           = ""
    custom_col_location_parts = []
    exclude_accessions        = []
    prune                     = False
    prune_snp_cutoff          = 5
    prune_date_resolution     = "month"
    pruning_method            = "fps"
    prune_max_representatives = 3
    prune_min_snp_difference  = 2
    prune_clade_snp_cutoff    = 5


# ---------------------------------------------------------------------------
# Provider presets
# ---------------------------------------------------------------------------
PATHOPLEXUS = {
    "col_specimen": "accessionVersion",
    "col_date":     "sampleCollectionDate",
    "location_parts": ["geoLocCity", "geoLocAdmin1", "geoLocCountry"],
    # Extra columns to carry forward into phymapr metadata (if present)
    "extra_cols": ["genotype", "geoLocCountry", "geoLocAdmin1", "geoLocCity"],
}


def resolve_columns(provider, meta_columns):
    """Return (col_specimen, col_date, location_parts) based on provider."""
    if provider == "pathoplexus":
        cfg = PATHOPLEXUS
        col_specimen = cfg["col_specimen"]
        col_date     = cfg["col_date"]
        # Only keep location parts that actually exist in the metadata
        location_parts = [c for c in cfg["location_parts"] if c in meta_columns]
    else:
        # custom provider
        col_specimen = custom_col_specimen
        col_date     = custom_col_date
        location_parts = [
            c for c in custom_col_location_parts if c in meta_columns
        ]

    return col_specimen, col_date, location_parts


def read_metadata(path):
    """Auto-detect delimiter from extension."""
    p = Path(path)
    sep = "\t" if p.suffix.lower() in (".tsv", ".txt") else ","
    df = pd.read_csv(path, sep=sep, low_memory=False)
    print(f"[prep] Loaded metadata: {len(df)} rows, {len(df.columns)} columns")
    return df


def extract_fasta_ids(fasta_path):
    """Return list of sequence IDs (everything after '>' up to first space)."""
    ids = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    print(f"[prep] Found {len(ids)} sequences in FASTA")
    return ids


def build_location_string(row, location_parts):
    """Concatenate non-empty location fields into a geocodable string."""
    parts = []
    for col in location_parts:
        val = row.get(col, "")
        if pd.notna(val) and str(val).strip():
            parts.append(str(val).strip())
    return ", ".join(parts)


def validate_alignment(fasta_ids, meta, col_specimen):
    fasta_set = set(fasta_ids)
    meta_set  = set(meta[col_specimen].dropna().astype(str))

    in_both   = fasta_set & meta_set
    only_fasta = fasta_set - meta_set
    only_meta  = meta_set - fasta_set

    print(f"[prep] Sequence <-> metadata match: {len(in_both)} matched")
    if only_fasta:
        print(
            f"[prep] WARNING: {len(only_fasta)} FASTA IDs have no metadata row "
            f"(first 5: {list(only_fasta)[:5]})"
        )
    if only_meta:
        print(
            f"[prep] INFO: {len(only_meta)} metadata rows have no FASTA sequence "
            f"(will be excluded from dates file)"
        )
    if not in_both:
        sys.exit(
            "[prep] ERROR: Zero sequences matched between FASTA headers and "
            f"metadata column '{col_specimen}'. "
            "Check provider setting and column names."
        )
    return in_both


def parse_midpoint_date(date_str):
    if pd.isna(date_str):
        return None
    val = str(date_str).strip()
    if not val or val.lower() in ("nan", "null", "none"):
        return None
    
    # Normalize delimiter
    val = val.replace("/", "-")
    parts = val.split("-")
    
    try:
        if len(parts) == 3:
            # YYYY-MM-DD
            year = int(parts[0])
            month = int(parts[1])
            day = int(parts[2])
            dt = datetime.date(year, month, day)
            return dt.strftime("%Y-%m-%d")
        elif len(parts) == 2:
            # YYYY-MM
            year = int(parts[0])
            month = int(parts[1])
            # Calculate midpoint of the month
            last_day = calendar.monthrange(year, month)[1]
            mid_day = last_day // 2 + 1
            dt = datetime.date(year, month, mid_day)
            return dt.strftime("%Y-%m-%d")
        elif len(parts) == 1 and len(val) == 4:
            # YYYY
            year = int(val)
            # July 2nd is the midpoint of a 365-day year
            dt = datetime.date(year, 7, 2)
            return dt.strftime("%Y-%m-%d")
        else:
            return val
    except Exception:
        # Fallback to original value if parsing fails
        return val



def filter_fasta(input_path, output_path, keep_ids):
    """Filter the input FASTA to only keep records matching keep_ids."""
    count = 0
    with open(input_path, "r") as inf, open(output_path, "w") as outf:
        current_seq = []
        current_id = None
        for line in inf:
            if line.startswith(">"):
                if current_id is not None:
                    if current_id in keep_ids:
                        outf.write(f">{current_id}\n" + "".join(current_seq))
                        count += 1
                # Extract new sequence ID (up to first space)
                current_id = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        # Don't forget last sequence!
        if current_id is not None and current_id in keep_ids:
            outf.write(f">{current_id}\n" + "".join(current_seq))
            count += 1
    print(f"[prep] Written filtered FASTA to {output_path} ({count} sequences)")


def get_date_group(date_str, resolution="month"):
    if pd.isna(date_str):
        return "unknown"
    val = str(date_str).strip()
    if not val or val.lower() in ("nan", "null", "none"):
        return "unknown"
    
    # Normalize delimiter
    val = val.replace("/", "-")
    parts = val.split("-")
    
    if len(parts) >= 1 and len(parts[0]) == 4:
        year = parts[0]
        month = parts[1] if len(parts) >= 2 else "07"
        day = parts[2] if len(parts) >= 3 else "15"
        
        if resolution == "year":
            return year
        elif resolution == "month":
            return f"{year}-{month.zfill(2)}"
        elif resolution == "week":
            try:
                dt = datetime.date(int(year), int(month), int(day))
                return f"{year}-W{dt.isocalendar()[1]:02d}"
            except Exception:
                return f"{year}-{month.zfill(2)}"
        elif resolution == "day":
            return f"{year}-{month.zfill(2)}-{day.zfill(2)}"
            
    return "unknown"


def compute_snp_distance_matrix(fasta_path, seq_ids=None):
    """
    Compute pairwise SNP distance matrix from an aligned FASTA.
    
    Strategy:
      1. Try pairsnp (fast C/scipy sparse backend)
      2. Fall back to numpy vectorized (pure Python, handles any numpy version)
    
    Returns:
      (list_of_ids, dict_of_dicts) where dict[id1][id2] = SNP distance
    """
    from Bio import SeqIO
    import numpy as np

    records = list(SeqIO.parse(fasta_path, "fasta"))
    if seq_ids is not None:
        seq_ids_set = set(seq_ids)
        records = [r for r in records if r.id in seq_ids_set]

    ids = [r.id for r in records]
    n = len(ids)

    if n == 0:
        return ids, {}

    # Try pairsnp first (fast sparse approach)
    try:
        import pairsnp
        # pairsnp.calculate_snp_matrix may fail on newer numpy; we'll patch in-memory
        # Write a temp fasta with just our subset if needed
        import tempfile, os
        tmp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
        for r in records:
            tmp_fasta.write(f">{r.id}\n{str(r.seq)}\n")
        tmp_fasta.close()

        # Attempt the call — if pairsnp's numpy usage is broken, catch it
        sparse_mat = pairsnp.calculate_snp_matrix(tmp_fasta.name)
        os.unlink(tmp_fasta.name)

        # Convert sparse to dense dict
        dense = sparse_mat.toarray() if hasattr(sparse_mat, 'toarray') else np.array(sparse_mat)
        # Make symmetric (pairsnp stores upper triangle)
        dense = dense + dense.T
        dist_dict = {s1: {s2: int(dense[i][j]) for j, s2 in enumerate(ids)} for i, s1 in enumerate(ids)}
        print(f"[prune] SNP distances computed via pairsnp ({n} sequences)")
        return ids, dist_dict

    except Exception as e:
        print(f"[prune] pairsnp unavailable or failed ({e}), using numpy fallback...")

    # Numpy fallback: vectorized pairwise distance
    print(f"[prune] Computing SNP distances with numpy ({n} sequences)...")

    # Encode sequences as uint8 arrays
    seq_arrays = []
    for r in records:
        s = np.frombuffer(str(r.seq).upper().encode('ascii'), dtype=np.uint8)
        seq_arrays.append(s)
    
    seqs = np.array(seq_arrays)  # shape: (n, alignment_length)
    
    # Valid positions: A=65, C=67, G=71, T=84
    valid_mask = np.isin(seqs, [65, 67, 71, 84])  # True where base is ACGT

    dist_dict = {s: {t: 0 for t in ids} for s in ids}
    
    # Compute pairwise distances (vectorized per-pair, not per-position)
    for i in range(n):
        for j in range(i + 1, n):
            both_valid = valid_mask[i] & valid_mask[j]
            d = int(np.sum((seqs[i] != seqs[j]) & both_valid))
            dist_dict[ids[i]][ids[j]] = d
            dist_dict[ids[j]][ids[i]] = d
    
    print(f"[prune] SNP distances computed via numpy fallback ({n} sequences)")
    return ids, dist_dict


def prune_sequences(meta_matched, fasta_file, col_specimen, col_date, location_parts, resolution, dates_out,
                    method, max_reps, min_snp_diff, clade_cutoff):
    """
    Cross-platform sequence pruning. Supports two methods:
      1. "fps" (Furthest-Point Sampling): groups by location/date, computes SNP distances,
         and keeps up to max_reps genetically diverse representatives.
      2. "clade" (Clade-Based Local Pruning): computes global SNP distances, performs 
         complete-linkage clustering, and keeps 1 representative per location/date/clade.
    
    SNP distances are computed using pairsnp (fast) or numpy (fallback).
    No external tools (snp-dists, WSL) required.
    """
    from Bio import SeqIO
    
    print(f"[prune] Starting sequence pruning using method='{method}' (date resolution='{resolution}')")
    
    # 1. Build location strings and date groups for all rows
    meta_matched = meta_matched.copy()
    meta_matched['_loc_str'] = meta_matched.apply(
        lambda r: build_location_string(r, location_parts), axis=1
    )
    meta_matched['_date_group'] = meta_matched[col_date].apply(
        lambda d: get_date_group(d, resolution)
    )
    
    # 2. Read FASTA sequences into memory
    print(f"[prune] Reading input FASTA {fasta_file} ...")
    seq_records = {r.id: r for r in SeqIO.parse(fasta_file, "fasta")}
    
    # Matched sequence IDs in metadata and FASTA
    all_matched_ids = [str(sid) for sid in meta_matched[col_specimen].dropna().unique() if str(sid) in seq_records]
    
    if not all_matched_ids:
        print("[prune] WARNING: No matched sequence IDs found. Nothing to prune.")
        return []

    # Helper function to compute number of Ns
    def get_n_count(sid):
        return str(seq_records[sid].seq).upper().count('N')
    
    # Write subset FASTA for distance computation
    import tempfile, os
    temp_dir = Path(dates_out).parent / "temp_prune"
    temp_dir.mkdir(parents=True, exist_ok=True)
    temp_fasta = temp_dir / "prune_temp.fasta"

    retained_ids = []
    
    if method == "fps":
        print(f"[prune] Running Furthest-Point Sampling (max_reps={max_reps}, min_snp_diff={min_snp_diff})")
        # Group by location and date
        groups = meta_matched.groupby(['_loc_str', '_date_group'])
        
        for (loc, dt_grp), df_grp in groups:
            group_ids = [str(sid) for sid in df_grp[col_specimen].dropna().unique() if str(sid) in seq_records]
            if not group_ids:
                continue
            if len(group_ids) <= max_reps:
                retained_ids.extend(group_ids)
                continue
                
            # Write group to temp FASTA
            from Bio import SeqIO as SeqIOWriter
            group_records = [seq_records[sid] for sid in group_ids]
            SeqIOWriter.write(group_records, str(temp_fasta), "fasta")
            
            # Compute SNP distance matrix
            try:
                _, dist_matrix = compute_snp_distance_matrix(str(temp_fasta))
                
                # Run FPS loop
                selected = []
                best_start = min(group_ids, key=get_n_count)
                selected.append(best_start)
                
                while len(selected) < max_reps:
                    max_min_d = -1
                    best_candidate = None
                    for cand in group_ids:
                        if cand in selected:
                            continue
                        min_d = min(dist_matrix[cand][sel] for sel in selected)
                        if min_d > max_min_d:
                            max_min_d = min_d
                            best_candidate = cand
                            
                    if best_candidate is not None and max_min_d >= min_snp_diff:
                        selected.append(best_candidate)
                    else:
                        break
                        
                retained_ids.extend(selected)
                
            except Exception as e:
                print(f"[prune] WARNING: Distance computation failed for group {loc} | {dt_grp} ({e}). Retaining all.")
                retained_ids.extend(group_ids)
                
    elif method == "clade":
        print(f"[prune] Running Clade-Based Complete-Linkage clustering (clade_cutoff={clade_cutoff})")
        # Write all matched sequences to temp FASTA
        from Bio import SeqIO as SeqIOWriter
        all_records = [seq_records[sid] for sid in all_matched_ids]
        SeqIOWriter.write(all_records, str(temp_fasta), "fasta")
        
        try:
            _, dist_matrix = compute_snp_distance_matrix(str(temp_fasta))
            
            # Complete linkage clustering
            clusters = [[name] for name in all_matched_ids]
            while True:
                min_dist = float('inf')
                to_merge = None
                for i in range(len(clusters)):
                    for j in range(i + 1, len(clusters)):
                        # Complete linkage: max distance between any pair across clusters
                        max_d = 0
                        for n1 in clusters[i]:
                            for n2 in clusters[j]:
                                d = dist_matrix.get(n1, {}).get(n2, 0)
                                if d > max_d:
                                    max_d = d
                        if max_d < min_dist:
                            min_dist = max_d
                            to_merge = (i, j)
                
                if min_dist <= clade_cutoff and to_merge is not None:
                    i, j = to_merge
                    clusters[i].extend(clusters[j])
                    clusters.pop(j)
                else:
                    break
                    
            print(f"[prune] Clustered {len(all_matched_ids)} sequences into {len(clusters)} distinct clades")
            
            # Map each sequence ID to its clade ID
            seq_to_clade = {}
            for cid, members in enumerate(clusters):
                for member in members:
                    seq_to_clade[member] = cid
                    
            # Add clade ID to metadata
            meta_matched = meta_matched.copy()
            meta_matched['_clade_id'] = meta_matched[col_specimen].apply(lambda sid: seq_to_clade.get(str(sid), -1))
            
            # Group by location, date, and clade ID — keep 1 representative per combo
            groups = meta_matched.groupby(['_loc_str', '_date_group', '_clade_id'])
            for keys, df_grp in groups:
                group_ids = [str(sid) for sid in df_grp[col_specimen].dropna().unique() if str(sid) in seq_records]
                if not group_ids:
                    continue
                best_member = min(group_ids, key=get_n_count)
                retained_ids.append(best_member)
                
        except Exception as e:
            print(f"[prune] ERROR: SNP distance/clustering failed ({e}). Retaining all sequences.")
            retained_ids = all_matched_ids
            
    else:
        print(f"[prune] Unknown method '{method}' - retaining all sequences")
        retained_ids = all_matched_ids
        
    # Clean up temp files
    try:
        if temp_fasta.exists():
            temp_fasta.unlink()
        if temp_dir.exists():
            temp_dir.rmdir()
    except Exception:
        pass
        
    print(f"[prune] Retained {len(retained_ids)} of {len(all_matched_ids)} matched sequences (pruned {len(all_matched_ids) - len(retained_ids)})")
    return retained_ids


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    meta      = read_metadata(metadata_file)
    fasta_ids = extract_fasta_ids(fasta_file)

    col_specimen, col_date, location_parts = resolve_columns(
        provider, list(meta.columns)
    )

    # Validate required columns exist
    for label, col in [("specimen", col_specimen), ("date", col_date)]:
        if col not in meta.columns:
            sys.exit(
                f"[prep] ERROR: {label} column '{col}' not found in metadata. "
                f"Available columns: {list(meta.columns)}"
            )

    matched_ids = validate_alignment(fasta_ids, meta, col_specimen)

    # Filter out excluded accessions
    matched_ids = [x for x in matched_ids if x not in exclude_accessions]

    # Keep only rows that have a matching FASTA sequence and are not excluded
    meta_matched = meta[meta[col_specimen].astype(str).isin(matched_ids)].copy()

    # Calculate midpoint date for each entry in the date column
    meta_matched[col_date] = meta_matched[col_date].apply(parse_midpoint_date)

    if prune:
        retained_ids = prune_sequences(
            meta_matched, fasta_file, col_specimen, col_date, location_parts,
            prune_date_resolution, dates_out,
            pruning_method, prune_max_representatives, prune_min_snp_difference, prune_clade_snp_cutoff
        )
        meta_matched = meta_matched[meta_matched[col_specimen].astype(str).isin(retained_ids)].copy()
        matched_ids = [x for x in matched_ids if x in retained_ids]


    # -------------------------------------------------------------------------
    # Output 1: dates.tsv  (name, date)  for TreeTime
    # -------------------------------------------------------------------------
    Path(dates_out).parent.mkdir(parents=True, exist_ok=True)
    dates_df = meta_matched[[col_specimen, col_date]].copy()
    dates_df.columns = ["name", "date"]
    dates_df = dates_df.dropna(subset=["date"])
    dates_df.to_csv(dates_out, sep="\t", index=False)
    print(f"[prep] Written dates file: {dates_out}  ({len(dates_df)} rows)")

    # -------------------------------------------------------------------------
    # Output 2: phymapr_metadata.tsv
    #   Standard output columns:
    #     accessionVersion   – specimen ID (matches tree tip labels)
    #     sampleCollectionDate – date
    #     location           – geocodable string for tidygeocoder
    #   Plus any provider-specific extras present in the data
    # -------------------------------------------------------------------------
    meta_matched = meta_matched.copy()

    # Build location column
    if location_parts:
        meta_matched["location"] = meta_matched.apply(
            lambda row: build_location_string(row, location_parts), axis=1
        )
    else:
        print(
            "[prep] WARNING: No location columns found/configured. "
            "'location' column will be empty – phymapr geocoding will fail."
        )
        meta_matched["location"] = ""

    # Normalise output column names so phymapr always sees the same names
    # regardless of provider.  The Snakefile hardcodes these names in params.
    rename_map = {col_specimen: "accessionVersion", col_date: "sampleCollectionDate"}
    meta_out_df = meta_matched.rename(columns=rename_map)

    # Core columns always written
    out_cols = ["accessionVersion", "sampleCollectionDate", "location"]

    # Carry forward provider-specific extras if present
    if provider == "pathoplexus":
        for extra in PATHOPLEXUS["extra_cols"]:
            mapped = rename_map.get(extra, extra)
            if mapped in meta_out_df.columns and mapped not in out_cols:
                out_cols.append(mapped)
    else:
        # For custom providers, carry everything except the renamed cols
        for col in meta_out_df.columns:
            if col not in out_cols:
                out_cols.append(col)

    # Only include columns that exist
    out_cols = [c for c in out_cols if c in meta_out_df.columns]

    Path(meta_out).parent.mkdir(parents=True, exist_ok=True)
    meta_out_df[out_cols].to_csv(meta_out, sep="\t", index=False)
    print(f"[prep] Written phymapr metadata: {meta_out}  ({len(meta_out_df)} rows)")
    print(f"[prep] Columns: {out_cols}")

    # Output 3: Filtered FASTA file
    filter_fasta(fasta_file, fasta_out, set(matched_ids))


main()
