import random
import os
import pandas as pd

def generate_virus_z_data():
    # Set seed for reproducibility
    random.seed(67890)
    
    # Define samples representing real-world surveillance imperfections:
    # 1. Crowded clusters (identical/near-identical genomes in one place/time)
    # 2. Unsampled intermediates (long time gaps and high mutation counts on branches)
    # 3. Unlinked global importations (highly drifted sequences with no intermediate link)
    # 4. Dead-end singletons (sampling failure / lost chains)
    
    # Columns: id, state, city, date, parent_id, lat, long, mutations_from_parent
    samples_info = [
        # --- NYC Root outbreak center (ground travel rail links) ---
        ("S00", "New York", "New York City", "2025-01-15", None, 40.7128, -74.0060, 0),
        ("S01", "New York", "New York City", "2025-01-17", "S00", 40.7128, -74.0060, 1), # dead-end singleton
        ("S02", "New Jersey", "Newark", "2025-01-20", "S00", 40.7357, -74.1724, 2),        # dead-end singleton
        
        # --- Atlanta Outbreak: Crowded Sequencing (superspreader event) ---
        # S_ATL_hub is the lineage root for Georgia.
        ("S_ATL_hub", "Georgia", "Atlanta", "2025-01-22", "S00", 33.7490, -84.3880, 2),
        # 9 sequences collected in a tight cluster (crowd), almost identical (0-1 SNPs)
        ("S_ATL_1", "Georgia", "Atlanta", "2025-02-15", "S_ATL_hub", 33.7490, -84.3880, 1),
        ("S_ATL_2", "Georgia", "Atlanta", "2025-02-15", "S_ATL_hub", 33.7490, -84.3880, 0),
        ("S_ATL_3", "Georgia", "Atlanta", "2025-02-16", "S_ATL_hub", 33.7490, -84.3880, 0),
        ("S_ATL_4", "Georgia", "Atlanta", "2025-02-16", "S_ATL_hub", 33.7490, -84.3880, 1),
        ("S_ATL_5", "Georgia", "Atlanta", "2025-02-17", "S_ATL_hub", 33.7490, -84.3880, 0),
        ("S_ATL_6", "Georgia", "Atlanta", "2025-02-17", "S_ATL_hub", 33.7490, -84.3880, 0),
        ("S_ATL_7", "Georgia", "Atlanta", "2025-02-18", "S_ATL_hub", 33.7490, -84.3880, 1),
        ("S_ATL_8", "Georgia", "Atlanta", "2025-02-18", "S_ATL_hub", 33.7490, -84.3880, 0),
        ("S_ATL_9", "Georgia", "Atlanta", "2025-02-19", "S_ATL_hub", 33.7490, -84.3880, 0),
        
        # --- Los Angeles Outbreak: Another Crowded Cluster ---
        # Seeded by flight from Atlanta
        ("S_LA_hub", "California", "Los Angeles", "2025-01-25", "S_ATL_hub", 34.0522, -118.2437, 2),
        # 7 sequences collected in a tight cluster (crowd), almost identical (0-1 SNPs)
        ("S_LA_1", "California", "Los Angeles", "2025-02-22", "S_LA_hub", 34.0522, -118.2437, 1),
        ("S_LA_2", "California", "Los Angeles", "2025-02-22", "S_LA_hub", 34.0522, -118.2437, 0),
        ("S_LA_3", "California", "Los Angeles", "2025-02-24", "S_LA_hub", 34.0522, -118.2437, 0),
        ("S_LA_4", "California", "Los Angeles", "2025-02-24", "S_LA_hub", 34.0522, -118.2437, 1),
        ("S_LA_5", "California", "Los Angeles", "2025-02-26", "S_LA_hub", 34.0522, -118.2437, 0),
        ("S_LA_6", "California", "Los Angeles", "2025-02-26", "S_LA_hub", 34.0522, -118.2437, 0),
        ("S_LA_7", "California", "Los Angeles", "2025-02-28", "S_LA_hub", 34.0522, -118.2437, 1),
        
        # --- Unsampled Intermediates (Gaps in surveillance) ---
        # 1. Chicago to Denver: Seeded from Chicago, but we miss IL/IA/NE ground travel.
        # Denver appears suddenly after 2 months with 14 mutations (drifted)
        ("S_CHI", "Illinois", "Chicago", "2025-01-22", "S00", 41.8781, -87.6298, 1),
        ("S_DEN", "Colorado", "Denver", "2025-03-25", "S_CHI", 39.7392, -104.9903, 14), 
        ("S_DEN_1", "Colorado", "Denver", "2025-03-28", "S_DEN", 39.7392, -104.9903, 1),
        ("S_DEN_2", "Colorado", "Denver", "2025-03-28", "S_DEN", 39.7392, -104.9903, 0),
        
        # 2. Seattle to Salt Lake City: Seeded from Seattle, but Idaho is completely unsampled.
        # Salt Lake City appears with 11 mutations
        ("S_SEA", "Washington", "Seattle", "2025-02-05", "S_LA_hub", 47.6062, -122.3321, 2),
        ("S_SLC", "Utah", "Salt Lake City", "2025-03-20", "S_SEA", 40.7608, -111.8910, 11),
        
        # --- Unlinked Global Importations (Independent introductions) ---
        # Florida samples appear with highly drifted genomes (26 mutations from NYC root),
        # showing an importation from outside the US surveillance pool (no direct US parent).
        ("S_MIA_import", "Florida", "Miami", "2025-03-01", "S00", 25.7617, -80.1918, 26),
        ("S_MIA_1", "Florida", "Miami", "2025-03-05", "S_MIA_import", 25.7617, -80.1918, 1),
        ("S_MIA_2", "Florida", "Miami", "2025-03-10", "S_MIA_import", 25.7617, -80.1918, 2),
        ("S_TPA", "Florida", "Tampa", "2025-03-12", "S_MIA_import", 27.9506, -82.4572, 4), # local spread of import
        
        # --- Dead-end / Sparse branches (Loss to follow-up) ---
        ("S_HOU", "Texas", "Houston", "2025-03-02", "S_ATL_hub", 29.7604, -95.3698, 8),
        ("S_DAL", "Texas", "Dallas", "2025-03-15", "S_HOU", 32.7767, -96.7970, 5),
    ]
    
    # Generate sequences
    seq_len = 5000
    bases = ['A', 'C', 'G', 'T']
    root_seq = [random.choice(bases) for _ in range(seq_len)]
    
    sequences = {}
    sequences["S00"] = "".join(root_seq)
    
    for sid, state, city, date, parent_id, lat, lon, branch_muts in samples_info:
        if parent_id is None:
            continue
        
        parent_seq = list(sequences[parent_id])
        mutated_seq = parent_seq.copy()
        
        # If branch_muts is greater than 0, introduce mutations
        if branch_muts > 0:
            mutation_positions = random.sample(range(seq_len), branch_muts)
            for pos in mutation_positions:
                orig_base = parent_seq[pos]
                new_base = random.choice([b for b in bases if b != orig_base])
                mutated_seq[pos] = new_base
                
        sequences[sid] = "".join(mutated_seq)
        
    # Write FASTA file
    output_dir = "C:/Code/data/virus_z"
    os.makedirs(output_dir, exist_ok=True)
    
    fasta_path = os.path.join(output_dir, "virus_z_sequences.fasta")
    with open(fasta_path, "w") as fh:
        for sid, seq in sequences.items():
            fh.write(f">{sid}\n{seq}\n")
            
    print(f"Generated FASTA: {fasta_path} ({len(sequences)} sequences)")
    
    # Write Metadata
    rows = []
    for sid, state, city, date, parent_id, lat, lon, branch_muts in samples_info:
        rows.append({
            "sample_id": sid,
            "collection_date": date,
            "city": city,
            "state": state,
            "country": "USA",
            "lat": lat,
            "long": lon
        })
        
    df = pd.DataFrame(rows)
    metadata_path = os.path.join(output_dir, "virus_z_metadata.csv")
    df.to_csv(metadata_path, index=False)
    print(f"Generated Metadata: {metadata_path} ({len(df)} rows)")

if __name__ == "__main__":
    generate_virus_z_data()
