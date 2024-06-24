import os

def gc_content(kmer):
    return (kmer.count('G') + kmer.count('C')) / len(kmer) * 100

def calculate_gc_content(dna_sequence, k):
    dna_sequence = dna_sequence.upper()
    kmers = [dna_sequence[i:i+k] for i in range(len(dna_sequence) - k + 1)]
    gc_contents = [gc_content(kmer) for kmer in kmers]
    overall_gc_content = sum(gc_contents) / len(gc_contents) if gc_contents else 0
    return overall_gc_content

def read_dna_sequence_from_file(file_path):
    with open(file_path, 'r') as file:
        dna_sequence = file.read().replace('\n', '').replace('\r', '')
    return dna_sequence

def process_files_in_directory(directory, k):
    for filename in os.listdir(directory):
        if filename.endswith('.fna'):
            file_path = os.path.join(directory, filename)
            dna_sequence = read_dna_sequence_from_file(file_path)
            overall_gc_content = calculate_gc_content(dna_sequence, k)
            
            output_filename = f"{os.path.splitext(filename)[0]}-GC.txt"
            output_file_path = os.path.join(directory, output_filename)
            
            with open(output_file_path, 'w') as output_file:
                output_file.write(f"Overall GC content: {overall_gc_content}%\n")
            
            print(f"Processed {filename}: Overall GC content = {overall_gc_content}%")

# Script configuration
directory = 'Genomes'
k = 4

# Process files in the specified directory
process_files_in_directory(directory, k)
