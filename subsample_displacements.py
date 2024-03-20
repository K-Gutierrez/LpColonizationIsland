# uses Python3
# By Ren Dodge, Ludington Lab at the Carnegie Institute for Science, 2023
# This script will batch process all .csv files in a designated subfolder and generate rose plots for the x-y coordinates
# There must be at least the number 'tracks_per_sample' with track length > 'min_track_length' and those tracks will be randomly subsampled
# so each data set contributes an equal number of tracks to total sample size
# .CSVs must have appropriate column naming, as from TrackMate 7 output 
# note that extraneous rows 2-4 in Trackmate output must be manually deleted prior to anlysis

import os
import subprocess
import pandas as pd
import sys
import random
import csv
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

#************** Define the folder path containing the CSV files
folder_path = '/Users/path/'
title = 'LpWFdelta'
min_track_length = 12
tracks_per_sample = 100

#=======================
#Functions
#=======================

def csv_df(csv_file, i):
    #import data - first argument
    df = pd.read_csv(csv_file)

    particles = pd.DataFrame(df)

    #Adding 10000*i to each track id so they become unique between .csvs, i must be passed from the for loop which calls this function. 
    particles['TRACK_ID']+= 10000*i

    # Create/reset a dictionary
    particles_dict = {}

    ## Iterate through unique TRACK_IDs
    for track_id in particles['TRACK_ID'].unique():
        # Filter rows for the current TRACK_ID and sort by POSITION_T
        track_data = particles[particles['TRACK_ID'] == track_id].sort_values(by='POSITION_T')

        # Extract POSITION_X and POSITION_Y coordinates
        x_values, y_values, t_values = zip(*track_data[['POSITION_X', 'POSITION_Y', 'POSITION_T']].values)

        # Subtract the first x and y values
        x_values = [x - x_values[0] for x in x_values]
        y_values = [y - y_values[0] for y in y_values]

        # Assign to the dictionary with TRACK_ID as the key
        particles_dict[track_id] = (x_values, y_values, t_values)

    return particles_dict

def subsample(particles_dict, min_len, nsamples):
    # Filter track IDs based on the length condition
    track_ids = []

    for track_id, data in particles_dict.items():
        x_values, y_values, t_values = data  # Unpack the tuple

        if len(x_values) >= min_len:
            track_ids.append(track_id)

    print("Total track IDs:", len(track_ids))  # Print the length of track IDs

    # Ensure that there are at least nsamples track IDs in the dictionary
    if len(track_ids) >= nsamples:
        # Subsample nsamples random track IDs
        random_track_ids = random.sample(track_ids, nsamples)
        print("Sampled track IDs:", len(random_track_ids))  # Print the length of track IDs

        # Access the corresponding data using the sampled track IDs
        sampled_data = {track_id: particles_dict[track_id] for track_id in random_track_ids}
        # 'sampled_data' now contains a dictionary with nsamples random track IDs and their corresponding data
        return sampled_data
    else:
        print(f"Subsampling error: There are fewer than {nsamples} track IDs in the dictionary.")

#END subsample function
        
def rose_plot (coordinates_dict):  #This function makes a ROSE PLOT accepts a dictionary as input 
    
    from matplotlib import cm
    # Create a figure of size 
    plt.figure(figsize=(3, 2.5), dpi=300)

    # Plot each track      
    for track_id, (x_values, y_values, _) in coordinates_dict.items():

        #Plot the x and y values with the calculated color
        plt.plot(x_values, y_values, label=f'TRACK_ID {track_id}', color='black', alpha=0.3, linewidth=1.0)

    # Add labels and legend (update title to argument variable someday)
    plt.title(f'{title}', loc='center')

    # Set limits
    limits = (-30, 30)

    # Draw axes and set limits
    for spine in ['right', 'top', 'bottom', 'left']:
        plt.gca().spines[spine].set_color('black')

    plt.axhline(0, color='grey', linewidth=0.5)  # Draw x-axis
    plt.axvline(0, color='grey', linewidth=0.5)  # Draw y-axis

    plt.xlim(limits)

    # Set x-axis ticks
    ticks = np.linspace(*limits, 3)
    labels = [f'{int(tick)}' for tick in ticks]
    plt.xticks(ticks, labels)

    # Add 'µm' label for x=0 tick
    labels[1] += ' µm'  # Assuming x=0 is the second tick
    plt.xticks(ticks, labels)

    plt.ylim(limits)
    plt.yticks(np.linspace(*limits, 3), [f'{int(tick)}' for tick in np.linspace(*limits, 3)])


    # Set the aspect ratio to be equal
    plt.gca().set_aspect('equal', adjustable='box')

    # Set ticks and labels color
    #ax = plt.gca()
    #ax.tick_params(axis='both', colors='black')

    # Set the aspect ratio to be equal
    plt.gca().set_aspect('equal', adjustable='box')

    # Adjust layout to prevent tick label overflow
    plt.tight_layout()


    # Save the figure
    output_path = Path("./output-png")
    output_path.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path / f'{title}_RosePlot.pdf')

    # if you want to show the result on the screen - turn off for batch processing
    #plt.show()

    #END rose_plot function

#=======================

#MAIN CODE
#=======================
        
        
# Change the current working directory
os.chdir(folder_path)

# Get a list of all CSV files in the folder
csv_files = [file for file in os.listdir(folder_path) if file.endswith('.csv')]

sampled_dict = {}

# Loop through each CSV file and process it using script1.py
i=1
for csv_file in csv_files:
    particles_dict = csv_df(csv_file, i)
   # Apply subsample function to particles_dict
    sampled_data = subsample(particles_dict, min_track_length, tracks_per_sample)
    # Update the sampled_dict with the subsampled data
    sampled_dict.update(sampled_data)
    i+=1

#make a graph
rose_plot (sampled_dict)

#output the .csv
# Specify the output file path
output_csv_path = f'{title}_subsampled_tracks.csv'

# Open the CSV file in write mode
with open(output_csv_path, 'w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)

    # Write the header row
    csv_writer.writerow(['TRACK_ID', 'POSITION_X', 'POSITION_Y', 'POSITION_T', 'FRAME'])

    # Iterate through particles_dict and write data rows
    for track_id, (x_values, y_values, t_values) in sampled_dict.items():
        # Repeat track_id for each corresponding x, y, t value
        rows = zip([track_id] * len(x_values), x_values, y_values, t_values, range(1, len(x_values) + 1))
        csv_writer.writerows(rows)


    
