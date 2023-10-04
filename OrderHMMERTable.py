# Read the table from a text file
with open('HMMER-1-3.txt', 'r') as file:
    # Create a list of tuples from the lines in the file
    table = [line.strip().split() for line in file]

# Sort the table by the second column in descending order
sorted_table = sorted(table, key=lambda x: int(x[1]), reverse=True)

# Print the sorted table
for row in sorted_table:
    print(row[0], row[1])

