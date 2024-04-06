# uses Python3
# By Ren Dodge, Ludington Lab at the Carnegie Institute for Science, 2023
# This script generates rose plots from track displacement data
# .CSVs must have appropriate column naming, as from TrackMate 7 output 
# note that extraneous rows 2-4 in Trackmate output must be manually deleted prior to analyis

# ==========================
# SECTION: Data Processing
# =========================

#import libraries and functions
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import random
from pathlib import Path


# ==========================
# SECTION2 : FUNTIONS
# =========================

def subsample (particles_dict, min_len, nsamples):
    #subsamples a dictionary with track length limit 
    # first limit the length of track in the dictionary to min_len (like 10 is a good number)
    for track_id, data in list(particles_dict.items()):
        #'x_values' is a key in the data dictionary
        if 'x_values' in data and len(data['x_values']) < min_len:
            # Remove the entry if the length of 'x_values' is less than min_len
    
    #Next, particles_dict contains only entries where the length of 'x_values' is min_len or greater
            del particles_dict[track_id]

    #Then will be a random subsampling of the data, to select nsamples random tracks
    import random
    track_ids = list(particles_dict.keys())

    # Ensure that there are at least nsamples track IDs in the dictionary
    if len(track_ids) >= nsamples:
        # Subsample nsamples random track IDs
        random_track_ids = random.sample(track_ids, nsamples)
        # Access the corresponding data using the sampled track IDs
        sampled_data = {track_id: particles_dict[track_id] for track_id in random_track_ids}
        # 'sampled_data' now contains a dictionary with 100 random track IDs and their corresponding data
        return sampled_data
    else:
        print("Subsampling error, there are fewer than min length:", min_len, " track IDs in the dictionary.")

    #END subsample function


def delta_t (first_key_data): #find the difference between first two t_values (frame length, assume its the same for all)  
 #first we get the first entry, third column/index 2 is t_values
    t_values = first_key_data[2]

    # Check if there are at least two elements in t_values before calculating the difference
    if len(t_values) >= 2:
        difference = t_values[1] - t_values[0]
    else:
      print("There are not enough elements in t_values.")
    return difference

def convert_seconds_to_mmss(seconds):
    minutes, seconds = divmod(seconds, 60)
    return f'{int(minutes):02d}:{int(seconds):02d}'
    

def rose_plot (coordinates_dict):  #This function makes a ROSE PLOT accepts a dictionary as input 
    
    from matplotlib import cm
    # Create a figure of size 
    plt.figure(figsize=(3, 2.5), dpi=300)

    # Plot each track      
    for track_id, (x_values, y_values, _) in coordinates_dict.items():

        #Plot the x and y values with the calculated color
        plt.plot(x_values, y_values, label=f'TRACK_ID {track_id}', color='black', alpha=0.3, linewidth=1.0)

    # Add labels and legend (update title to argument variable someday)
        plt.title(f'{title}', loc='center', fontsize=8)


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



#============================
# MAIN Script
#============================
    

#import data - first argument
df = pd.read_csv(sys.argv[1])

filename = sys.argv[1][:-4]
title = sys.argv[2][:-4]
filename = filename.split("/", 1)[-1]

particles = pd.DataFrame(df)

# Create a dictionary
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


#Convert particles_dict to a DataFrame to make the CSV
#df_particles = pd.DataFrame([(track_id, x, y, t) for track_id, (x, y, t) in particles_dict.items()], columns=['TRACK_ID', 'POSITION_X', 'POSITION_Y', 'POSITION_T'])

# Save the DataFrame to a .csv file
#csv_file_path = "particles_dict.csv"
#df_particles.to_csv(csv_file_path, index=False)

#sampled_data = subsample(particles_dict, 5, 50)

rose_plot(particles_dict)




