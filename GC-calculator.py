def calculate_gc_content(dna_sequence):
    # Convert the sequence to uppercase to ensure consistency
    dna_sequence = dna_sequence.upper()
    
    # Count the number of G, C, A, and T
    g_count = dna_sequence.count('G')
    c_count = dna_sequence.count('C')
    a_count = dna_sequence.count('A')
    t_count = dna_sequence.count('T')
    
    # Calculate the total count of nucleotides
    total_count = a_count + t_count + g_count + c_count
    
    # Calculate the GC content
    gc_content = (g_count + c_count) / total_count * 100
    
    return gc_content

# Example usage
dna_sequence = "AGCTATAG"
gc_content = calculate_gc_content(dna_sequence)
print(f"The GC content of the sequence is {gc_content:.2f}%")
