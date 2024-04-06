# Read Document A and extract the IDs
with open('Document_A.txt', 'r') as file_a:
    ids_a = [line.strip() for line in file_a]

# Read Document B and extract the IDs
with open('Final-gftB.fasta', 'r') as file_b:
    ids_b = [line.strip()[1:].split('_')[0] for line in file_b if line.startswith('>')]

# Find missing IDs
missing_ids = set(ids_a) - set(ids_b)

# Find repeated IDs
repeated_ids = [id for id in ids_b if ids_b.count(id) > 1]

# Print missing IDs
print("Missing IDs in Document B:")
for missing_id in missing_ids:
    print(missing_id)

# Print repeated IDs
print("\nRepeated IDs in Document B:")
for repeated_id in set(repeated_ids):
    print(repeated_id)
