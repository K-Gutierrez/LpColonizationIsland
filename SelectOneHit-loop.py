import pandas as pd
import os

# Directory containing input files
input_directory = '/mnt/sequence/kgutierrez/HuntingColIslands/'

# Directory where output files will be saved
output_directory = '/mnt/sequence/kgutierrez/HuntingColIslands/gftA/OneHitPerGenera/'

# List of files in the input directory
input_files = os.listdir(input_directory)

# Iterate over each input file
for input_file in input_files:
    # Build the full path of the input file
    input_path = os.path.join(input_directory, input_file)

    # Load data from the input file
    df = pd.read_csv(input_path, delimiter='\t')

    # Handle the case of a single record
    if len(df) == 1:
        # Simply copy that single record
        selected_row = df.iloc[0]
    else:
        # Find the minimum number in the third column if present
        if 'genome.contigs' in df.columns:
            min_contigs = df['genome.contigs'].min()
            min_contigs_rows = df[df['genome.contigs'] == min_contigs]

            # If there are more than one row with the minimum number, choose the one with the maximum value in the fourth column if present
            if 'genome.genome_length' in df.columns:
                selected_row = min_contigs_rows.loc[min_contigs_rows['genome.genome_length'].idxmax()]
            else:
                # If the column is not present, simply select the first row
                selected_row = min_contigs_rows.iloc[0]
        else:
            # If the column is not present, simply select the first row
            selected_row = df.iloc[0]

    # Create the name of the output file
    output_file_name = f"{os.path.splitext(input_file)[0]}-OneHit.txt"
    output_path = os.path.join(output_directory, output_file_name)

    # Write the selected row to the output file on a single line
    with open(output_path, 'w') as f:
        f.write('\t'.join(selected_row.astype(str).tolist()))
