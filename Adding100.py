# Define the input and output file paths
input_file = 'your_table.txt'
output_file = 'updated_table.txt'

# Read the file and update the values
with open(input_file, 'r') as input, open(output_file, 'w') as output:
    for line in input:
        column1, column2 = line.strip().split('\t')
        if not column2:
            column2 = '100'
        output.write(f'{column1}\t{column2}\n')

print(f'The updated table has been saved to {output_file}.')
