import pandas as pd
from Bio import SeqIO

# Load accession versions
meta_df = pd.read_csv("C:/Code/data/me_adhs/final_meta.csv")
accessions = set(meta_df["accessionVersion"].dropna().astype(str))

# Filter FASTA
input_fasta = "C:/Code/data/pathoplexus/measles_aligned-nuc_2026-06-24T0256.fasta"
output_fasta = "C:/Code/data/me_adhs/me_adhs_aligned.fasta"

filtered_records = []
for record in SeqIO.parse(input_fasta, "fasta"):
    seq_id = record.id.split()[0]
    if seq_id in accessions:
        filtered_records.append(record)

# Write output FASTA
SeqIO.write(filtered_records, output_fasta, "fasta")
print(f"Written {len(filtered_records)} sequences to {output_fasta}")
