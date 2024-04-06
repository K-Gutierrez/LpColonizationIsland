# uses Python3
# By Ren Dodge, Ludington Lab at the Carnegie Institute for Science, 2023
# This script generate coefficients of diffusion .csv output and a histogram
# .CSVs must have appropriate column naming, as from TrackMate 7 output 
# note that extraneous rows 2-4 in Trackmate output must be manually deleted prior to anlysis

###########################################
########### set up script #################
###########################################
#import libraries and functions
import matplotlib
import matplotlib.pyplot as plt
import numpy
import pandas
from scipy.optimize import curve_fit
import sys
from tqdm import tqdm
from pathlib import Path

#####################################################
###>>>>>>>>> GLOBAL VARIABLE DEFINITIONS <<<<<<<<<###
###>>>>>>>>>>>>> CHANGE HERE  <<<<<<<<<<<<<##########
cutoff_length = 12 #minimum length of track
tau = 10 #max tau length (consider tau = cutoff_length - 2?)


####################################################
############## FUNCTION DEFINITIONS ################
####################################################

######## simple math ########
def calculate_sq_displacement(x0, y0, x1, y1):
    '''Calculate squared distance (i.e. Pythagorean distance ** 2)'''
    sq_displacement = (x1-x0)**2 + (y1-y0)**2
    return sq_displacement

def calculate_diffusion_stats(list_of_coefficients):
    '''Calculate mean and stdev (assuming normal distribution i think?)'''
    mean_diffusion = numpy.mean(list_of_coefficients)
    stdev_diffusion = numpy.std(list_of_coefficients)
    return (mean_diffusion, stdev_diffusion)

######## functions that assist in the for-loops #######
def find_timesteps(len_trace, dt=10, limit=True):
    if limit:
        if len_trace >= 10:
            max_dt = dt
        else:
            max_dt = len_trace - 1
    else:
        max_dt = len_trace - 1
    return(max_dt)

def displacement_of_TRACK_ID(df_id, max_dt):
    '''Calculates squared displacement to get MSD later'''
    displacements_dict = {}
    dt_displacements = []
    dt = 1
    while dt < max_dt+1:
        #print("Current dt being calculated:", dt)
        dt_displacements = []
        idx = 0
        while idx < len(df_id.index)-dt:
            x0, y0, x1, y1 = df_id.iloc[idx]["POSITION_X"], df_id.iloc[idx]["POSITION_Y"], df_id.iloc[idx+dt]["POSITION_X"], df_id.iloc[idx+dt]["POSITION_Y"]
            displacement = calculate_sq_displacement(float(x0), float(y0), float(x1), float(y1))
            #print(x0, y0, x1, y1, displacement)
            dt_displacements.append(displacement)
            idx += 1
        displacements_dict[dt] = dt_displacements
        #print("Done with time step", dt, "of TRACK_ID", id)
        dt += 1
    return displacements_dict

def fill_gaps(frames_list):
    ''' Checks if there are holes in frames_list = df_id["FRAME"] '''

    no_gaps_list = list(range(min(frames_list), max(frames_list)+1))
    missing_vals = list(set(no_gaps_list).difference(frames_list))
    missing_vals.sort()
    #print("These are the missing frames:", missing_vals)
    #for missing_vals, add items in dataframe like "null"
    new_list = []
    l = frames_list
    for val in missing_vals:
        missing_index = no_gaps_list.index(val)
        l_until_missing = l[:missing_index].to_list()
        # print("This is the current missing index: ", missing_index)
        # print("This is the list up to the missing index: ", l_until_missing)

        l = l[missing_index:]
        # print("This is the new list after removing up to the missing index: ", l)
        no_gaps_list = no_gaps_list[missing_index+1:]
        # print("This is the new no-gaps list after removing up to the missing index: ",no_gaps_list)

        new_list += l_until_missing
        new_list.append(numpy.nan)
        # print("This is the current new list: ", new_list)

    new_list += no_gaps_list
    #print("This is the final new list: ", new_list)
    return new_list

def put_dictdata_in_dictionary(container_dict, data_dict, output):
    '''Accesses displacements in data_dict, calculates the mean and stdevs, and puts them into container_dict'''
    for key in data_dict.keys():
        if output == "mean":
            mean_displacement = numpy.nanmean(data_dict[key]) #calculate average of mean-squared displacements per time step
            if key in container_dict.keys():
                container_dict[key].append(mean_displacement)
            else:
                container_dict[key] = mean_displacement
        elif output == "stdev":
            stdev_displacement = numpy.nanstd(data_dict[key])
            if key in container_dict.keys():
                container_dict[key].append(stdev_displacement)
            else:
                container_dict[key] = stdev_displacement
        else:
            print("Output data variable ('mean' or 'stdev') not specified!")
            break
    return container_dict

def format_data(all_particle_vals_dict, all_particle_err_dict, TRACK_ID):
    '''
    For every key ("TRACK_ID") in all_particle_dict, extracts keys and values and formats them in a way
    that is useable for scipy.curve_fit(). This entails putting the keys and values in arrays and removing NaNs.
    This should also be used for formatting the dictionary containing the error (i.e. using the second argument)

    '''
    x, y = numpy.array(list(all_particle_vals_dict[TRACK_ID].keys())), numpy.array(list(all_particle_vals_dict[TRACK_ID].values()))
    error = numpy.array(list(all_particle_err_dict[TRACK_ID].values()))
    nan_idxs = numpy.argwhere(numpy.isnan(y))
    x = numpy.delete(x, nan_idxs)
    y = numpy.delete(y, nan_idxs)
    error = numpy.delete(error, nan_idxs)

    return x, y, error


##########################################
######## end function definitions ########
##########################################


#import data from first argument .csv file
df = pandas.read_csv(sys.argv[1])

#get the filename with truncated .csv extension for the output

filename = sys.argv[1][:-4]  
filename = filename.split("/")[-1]

#slice out the columns we care about
particles = df[["TRACK_ID", "POSITION_X", "POSITION_Y", "POSITION_T", "FRAME"]]

#sort by time
particles = particles.sort_values(by=["POSITION_T"])

#calculate the frame length 
times = pandas.Series(particles["POSITION_T"].unique())
t = times[1]-times[0]


##########################################################
############### main for loops ############################
##########################################################

#get list of unique particle ids
particle_ids = set(particles["TRACK_ID"])
particle_ids = list(particle_ids)

#set up dictionary to contain all particle displacements, and another to contain trace lengths for each particle
all_particle_mean_displacement_dict = {}
all_particle_stdev_displacement_dict = {}
particle_trace_len = {}

stop = 0 #for testing code

for id in tqdm(particle_ids): #iterate through TRACK_IDs of the dataframe we sliced out
###how this works is:
###1. check if the TRACK_ID is "None", skip the particle if it is
###2. for each TRACK_ID, extract the rows that correspond to that TRACK_ID and calculate the number of frames for that trace
###3. calculate squared displacement for each time step
###4. for each time step, average (get the mean of) the squared displacement and save into dictionary
### - data in dictionary looks like: {TRACK_ID:[<x^2> for tau=1, <x^2> for tau=2, <x^2> for tau=3, etc]}
    
    if id == "None":
        continue
    #print("Current TRACK_ID:", id)
    displacements_dict = {}

    df_id = particles.loc[particles["TRACK_ID"] == id]
    len_trace = len(df_id.index)

    #print("length of trace:", len_trace, "frames")
    if len_trace < cutoff_length:
        continue
    particle_trace_len[id] = len_trace

    displacements_dict = displacement_of_TRACK_ID(df_id, find_timesteps(len_trace, tau))
    #print(displacements_dict)

    ###second bit of code -
    #set up dictionaries to contain displacements and stdevs for a single particle
    single_particle_disp_dict = {}
    single_particle_stdev_dict = {}
    put_dictdata_in_dictionary(single_particle_disp_dict, displacements_dict, 'mean')
    put_dictdata_in_dictionary(single_particle_stdev_dict, displacements_dict, 'stdev')
    all_particle_mean_displacement_dict[id] = single_particle_disp_dict
    all_particle_stdev_displacement_dict[id] = single_particle_stdev_dict

    #print("Track_ID", id, "done")

    ## uncomment below to test code
    # stop += 1
    #
    # if stop == 21:
    #     break


####################################################
########### graph diffusion coefficient#############
####################################################
def fit_func(x, a):
    # Curve fitting function: msd = 2nD*t
    return 2*2*a*x*t   # b=0 is implied; n=2, a=D, x=timestep, t=frame length

# initialize containers
D_coeffs = []
positive_IDs = []
negative_IDs = []
trace_len = []
r_squared_list = []

###!!!
#matplotlib.use('Agg')  #makes plots non interactive for BATCH processing

# initialize plot
plt.figure(1)
plt.title("Mean-squared displacement per tau")
plt.xlabel("tau")
plt.ylabel("MSD")
plt.close()

for TRACK_ID in all_particle_mean_displacement_dict.keys(): #use numpy to fit a line through the squared displacements at each tau
    x, y, error = format_data(all_particle_mean_displacement_dict, all_particle_stdev_displacement_dict, TRACK_ID)
    plt.plot(x, y)
    #plt.errorbar(x, y, error)
    popt, pcov = curve_fit(f=fit_func, xdata=x, ydata=y, sigma=error)

    # for r^2, followed https://stackoverflow.com/questions/19189362/getting-the-r-squared-value-using-curve-fit
    residuals = y - fit_func(x, *popt)
    ss_res = numpy.sum(residuals**2)
    ss_tot = numpy.sum((y - numpy.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    r_squared_list.append(r_squared)

    # add values to lists used for data analysis output
    D_coeffs.append(float(popt[0]))
    trace_len.append(particle_trace_len[TRACK_ID])

    if float(popt[0]) < 0:
        negative_IDs.append(TRACK_ID)
    else:
        positive_IDs.append(TRACK_ID)

mean_D_coeff = calculate_diffusion_stats(D_coeffs)[0]
stdev_D_coeff = calculate_diffusion_stats(D_coeffs)[1]
print("Mean diffusion coeff:", mean_D_coeff, "Std dev:", stdev_D_coeff)

from matplotlib import cm

# Create a figure of size
plt.figure(figsize=(3, 2.67), dpi=300)  # Adjust the figure size as needed

# Plot the normalized histogram manually
counts, bins, _ = plt.hist(D_coeffs, bins=20, range=(0, 0.3), color="dimgray", alpha=0.7)
plt.clf()

# Manually normalize the heights
normalized_counts = counts / sum(counts)

# Plot the normalized histogram
plt.bar(bins[:-1], normalized_counts, width=(0.3 / 20), color="red", alpha=0.5, edgecolor="black")

# Set y-axis range from 0 to 1
plt.ylim(0, 1)

# Set x-axis label
plt.xlabel("Coefficient of Diffusion (µm²/s)")

# Adjusted title with smaller font size
title = (
    f'{filename}\n dt={round(t, 3)}, N tracks={len(trace_len)}, tau={tau}\n'
    f'MeanDcoeff={round(mean_D_coeff, 3)}, StdDev={round(stdev_D_coeff, 3)}'
)
plt.title(title, fontweight='bold', fontsize=10)  # Adjust the fontsize as needed
plt.tight_layout()

# Create the output directory for PDF if it doesn't exist
output_pdf_dir = Path("./output-pdf")
output_pdf_dir.mkdir(parents=True, exist_ok=True)

# Save the figure as a PDF
plt.savefig(output_pdf_dir / f'{filename}_histogram.pdf')
plt.show()

# Output data
dict_with_ID_and_trace_len = {
    f'{filename} TRACK_ID': positive_IDs,
    f'{filename} D_coeffs': D_coeffs,
    f'{filename} r_squared': r_squared_list,
    f'{filename} trace_len': trace_len
}

# Create the output directory for CSV if it doesn't exist
output_csv_dir = Path("./output-csv")
output_csv_dir.mkdir(parents=True, exist_ok=True)

# Save the DataFrame as a CSV file using the appropriate path
pandas.DataFrame(dict_with_ID_and_trace_len).to_csv(output_csv_dir / f'D_coeffs_{filename}.csv')

