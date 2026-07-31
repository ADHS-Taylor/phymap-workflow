import random
import os
import pandas as pd

def generate_hypothetical_data():
    # Set seed for reproducibility
    random.seed(42)
    
    # 24 unique US states and coordinates (capturing capitals/major cities)
    # The order below defines a realistic progression: East Coast -> South -> Midwest -> West -> West Coast
    samples_info = [
        # id, state, city, date, parent_id, lat, long
        ("S00", "New York", "New York City", "2025-01-15", None, 40.7128, -74.0060),
        ("S01", "Massachusetts", "Boston", "2025-02-10", "S00", 42.3601, -71.0589),
        ("S02", "Pennsylvania", "Philadelphia", "2025-02-28", "S00", 39.9526, -75.1652),
        ("S03", "New Jersey", "Newark", "2025-03-15", "S02", 40.7357, -74.1724),
        ("S04", "Maryland", "Baltimore", "2025-04-05", "S03", 39.2904, -76.6122),
        ("S05", "Virginia", "Richmond", "2025-04-20", "S04", 37.5407, -77.4360),
        ("S06", "North Carolina", "Raleigh", "2025-05-15", "S05", 35.7796, -78.6382),
        ("S07", "Georgia", "Atlanta", "2025-06-10", "S06", 33.7490, -84.3880),
        ("S08", "Florida", "Miami", "2025-07-01", "S07", 25.7617, -80.1918),
        ("S09", "Alabama", "Birmingham", "2025-07-20", "S07", 33.5186, -86.8104),
        ("S10", "Tennessee", "Nashville", "2025-08-05", "S06", 36.1627, -86.7816),
        ("S11", "Louisiana", "New Orleans", "2025-08-25", "S09", 29.9511, -90.0715),
        ("S12", "Texas", "Houston", "2025-09-15", "S11", 29.7604, -95.3698),
        ("S13", "Florida", "Orlando", "2025-08-15", "S08", 28.5383, -81.3792),
        ("S14", "New York", "Buffalo", "2025-03-20", "S01", 42.8864, -78.8784),
        ("S15", "Ohio", "Columbus", "2025-03-01", "S00", 39.9612, -82.9988),
        ("S16", "Michigan", "Detroit", "2025-03-25", "S15", 42.3314, -83.0458),
        ("S17", "Illinois", "Chicago", "2025-04-15", "S16", 41.8781, -87.6298),
        ("S18", "Missouri", "St. Louis", "2025-05-01", "S17", 38.6270, -90.1994),
        ("S19", "Colorado", "Denver", "2025-06-01", "S18", 39.7392, -104.9903),
        ("S20", "Utah", "Salt Lake City", "2025-06-20", "S19", 40.7608, -111.8910),
        ("S21", "Arizona", "Phoenix", "2025-07-10", "S20", 33.4484, -112.0740),
        ("S22", "Nevada", "Las Vegas", "2025-07-30", "S21", 36.1716, -115.1398),
        ("S23", "California", "Los Angeles", "2025-08-15", "S22", 34.0522, -118.2437),
        ("S24", "Oregon", "Portland", "2025-09-01", "S23", 45.5152, -122.6784),
        ("S25", "Washington", "Seattle", "2025-09-20", "S24", 47.6062, -122.3321),
        ("S26", "California", "San Francisco", "2025-09-10", "S23", 37.7749, -122.4194),
        ("S27", "Texas", "Dallas", "2025-10-05", "S12", 32.7767, -96.7970),
        ("S28", "Washington", "Spokane", "2025-10-15", "S25", 47.6588, -117.4260),
        ("S29", "Illinois", "Springfield", "2025-05-10", "S17", 39.7817, -89.6501),
        ("S30", "Colorado", "Colorado Springs", "2025-07-05", "S19", 38.8339, -104.8214),
        ("S31", "Arizona", "Tucson", "2025-08-20", "S21", 32.2226, -110.9747)
    ]
    
    # 1. Generate sequences
    seq_len = 5000
    bases = ['A', 'C', 'G', 'T']
    root_seq = [random.choice(bases) for _ in range(seq_len)]
    
    sequences = {}
    sequences["S00"] = "".join(root_seq)
    
    # Generate descendant sequences by introducing mutations
    # The list is already sorted in topological order (parents defined before children)
    for sid, state, city, date, parent_id, lat, lon in samples_info:
        if parent_id is None:
            continue
        
        # Mutate from parent
        parent_seq = list(sequences[parent_id])
        mutated_seq = parent_seq.copy()
        
        # Introduce 12 random single-nucleotide substitutions per branch
        num_mutations = 12
        mutation_positions = random.sample(range(seq_len), num_mutations)
        for pos in mutation_positions:
            orig_base = parent_seq[pos]
            new_base = random.choice([b for b in bases if b != orig_base])
            mutated_seq[pos] = new_base
            
        sequences[sid] = "".join(mutated_seq)
    
    # Write FASTA file
    output_dir = "C:/Code/data/virus_x"
    os.makedirs(output_dir, exist_ok=True)
    
    fasta_path = os.path.join(output_dir, "virus_x_sequences.fasta")
    with open(fasta_path, "w") as fh:
        for sid, seq in sequences.items():
            fh.write(f">{sid}\n{seq}\n")
            
    print(f"Generated FASTA: {fasta_path} ({len(sequences)} sequences)")
    
    # 2. Generate Metadata
    rows = []
    for sid, state, city, date, parent_id, lat, lon in samples_info:
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
    metadata_path = os.path.join(output_dir, "virus_x_metadata.csv")
    df.to_csv(metadata_path, index=False)
    print(f"Generated Metadata: {metadata_path} ({len(df)} rows)")

if __name__ == "__main__":
    generate_hypothetical_data()
