import random
import os
import pandas as pd

def generate_virus_y_data():
    # Set seed for reproducibility
    random.seed(12345)
    
    # Define samples capturing hub-and-spoke (airline) and local (ground) spread
    # id, state, city, date, parent_id, lat, long, mutations_from_parent
    samples_info = [
        # Root outbreak center (New York City)
        ("S00", "New York", "New York City", "2025-01-15", None, 40.7128, -74.0060, 0),
        
        # ---------------------------------------------------------
        # Primary Airline Hubs (rapid, long-distance jumps from NYC)
        # ---------------------------------------------------------
        ("S01", "Georgia", "Atlanta", "2025-01-20", "S00", 33.7490, -84.3880, 2),
        ("S02", "Illinois", "Chicago", "2025-01-22", "S00", 41.8781, -87.6298, 1),
        ("S03", "California", "Los Angeles", "2025-01-25", "S00", 34.0522, -118.2437, 2),
        ("S04", "Texas", "Dallas", "2025-01-28", "S01", 32.7767, -96.7970, 2),
        ("S05", "Colorado", "Denver", "2025-02-02", "S02", 39.7392, -104.9903, 1),
        ("S06", "Washington", "Seattle", "2025-02-05", "S03", 47.6062, -122.3321, 2),
        ("S07", "Florida", "Miami", "2025-02-08", "S01", 25.7617, -80.1918, 2),
        
        # ---------------------------------------------------------
        # Regional Spreads (local diffusion via ground travel)
        # ---------------------------------------------------------
        # Local spread from Atlanta (GA) hub
        ("S08", "Alabama", "Birmingham", "2025-02-12", "S01", 33.5186, -86.8104, 5),
        ("S09", "Tennessee", "Nashville", "2025-02-18", "S01", 36.1627, -86.7816, 6),
        ("S10", "North Carolina", "Charlotte", "2025-02-25", "S01", 35.2271, -80.8431, 4),
        
        # Local spread from Chicago (IL) hub
        ("S11", "Wisconsin", "Milwaukee", "2025-02-15", "S02", 43.0389, -87.9065, 4),
        ("S12", "Indiana", "Indianapolis", "2025-02-22", "S02", 39.7684, -86.1581, 5),
        ("S13", "Michigan", "Detroit", "2025-03-01", "S02", 42.3314, -83.0458, 6),
        
        # Local spread from Los Angeles (CA) hub
        ("S14", "California", "San Diego", "2025-02-10", "S03", 32.7157, -117.1611, 4),
        ("S15", "California", "Riverside", "2025-02-14", "S03", 33.9533, -117.3962, 5),
        ("S16", "Nevada", "Las Vegas", "2025-02-20", "S03", 36.1716, -115.1398, 6),
        ("S17", "Arizona", "Phoenix", "2025-02-28", "S03", 33.4484, -112.0740, 7),
        
        # Local spread from Dallas (TX) hub
        ("S18", "Texas", "Houston", "2025-02-18", "S04", 29.7604, -95.3698, 5),
        ("S19", "Texas", "Austin", "2025-02-24", "S04", 30.2672, -97.7431, 4),
        ("S20", "Oklahoma", "Oklahoma City", "2025-03-05", "S04", 35.4676, -97.5164, 6),
        
        # Local spread from Denver (CO) hub
        ("S21", "Colorado", "Colorado Springs", "2025-02-20", "S05", 38.8339, -104.8214, 4),
        ("S22", "Utah", "Salt Lake City", "2025-03-02", "S05", 40.7608, -111.8910, 5),
        
        # Local spread from Seattle (WA) hub
        ("S23", "Oregon", "Portland", "2025-02-25", "S06", 45.5152, -122.6784, 4),
        ("S24", "Washington", "Spokane", "2025-03-10", "S06", 47.6588, -117.4260, 5),
        
        # Local spread from Miami (FL) hub
        ("S25", "Florida", "Orlando", "2025-02-28", "S07", 28.5383, -81.3792, 4),
        ("S26", "Florida", "Tampa", "2025-03-08", "S07", 27.9506, -82.4572, 5),
        
        # Local East Coast spreads from NYC directly (regional rail / road travel)
        ("S27", "Pennsylvania", "Philadelphia", "2025-01-25", "S00", 39.9526, -75.1652, 3),
        ("S28", "New Jersey", "Newark", "2025-01-18", "S00", 40.7357, -74.1724, 2),
        ("S29", "Massachusetts", "Boston", "2025-02-05", "S00", 42.3601, -71.0589, 4),
    ]
    
    # 1. Generate sequences based on the topology
    seq_len = 5000
    bases = ['A', 'C', 'G', 'T']
    root_seq = [random.choice(bases) for _ in range(seq_len)]
    
    sequences = {}
    sequences["S00"] = "".join(root_seq)
    
    # Generate descendant sequences
    for sid, state, city, date, parent_id, lat, lon, branch_muts in samples_info:
        if parent_id is None:
            continue
        
        parent_seq = list(sequences[parent_id])
        mutated_seq = parent_seq.copy()
        
        # Introduce the specified number of mutations
        mutation_positions = random.sample(range(seq_len), branch_muts)
        for pos in mutation_positions:
            orig_base = parent_seq[pos]
            new_base = random.choice([b for b in bases if b != orig_base])
            mutated_seq[pos] = new_base
            
        sequences[sid] = "".join(mutated_seq)
        
    # Write FASTA file
    output_dir = "C:/Code/data/virus_y"
    os.makedirs(output_dir, exist_ok=True)
    
    fasta_path = os.path.join(output_dir, "virus_y_sequences.fasta")
    with open(fasta_path, "w") as fh:
        for sid, seq in sequences.items():
            fh.write(f">{sid}\n{seq}\n")
            
    print(f"Generated FASTA: {fasta_path} ({len(sequences)} sequences)")
    
    # 2. Generate Metadata
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
    metadata_path = os.path.join(output_dir, "virus_y_metadata.csv")
    df.to_csv(metadata_path, index=False)
    print(f"Generated Metadata: {metadata_path} ({len(df)} rows)")

if __name__ == "__main__":
    generate_virus_y_data()
