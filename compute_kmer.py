#! /usr/bin/env python3

"""
Created by Dana Rocha 2/24/21
"""

import sys
import argparse
import gzip
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

def main(infile):

    fh_in = open(infile, "r")

    kmer_size = int(input("Please enter a k-mer size (31, 41, 51, 61, 71, 81, 91, or 101): "))

    checked = check_size(kmer_size)

    kmers = calc_kmer(fh_in, checked)

    kmer_occurs = count_occurrence(kmers)

    xvals_array, yvals_array, min_pos, real_min_height, occur_sorted_df = calculate_minimas(kmer_occurs)

    threshold_value, merged_data, x_thresh = find_threshold(min_pos, real_min_height, occur_sorted_df)

    remaining_kmers = find_remaining_kmers(merged_data)

    generate_plot(xvals_array, yvals_array, min_pos, real_min_height, x_thresh, threshold_value, remaining_kmers)

    print(threshold_value)
    print(remaining_kmers)

def check_size(size_input):
    """
    This function checks if the user selected one of the allowed k-mer sizes.
    :param size_input:
    :return:
    """

    try:
        sizes = [31, 41, 51, 61, 71, 81, 91, 101]

        if size_input in sizes:
            return size_input
        else:
            raise ValueError
    except ValueError:
        print("Invalid k-mer size entered.")
        sys.exit(1)


def calc_kmer(list_of_seqs, ksize):

    kmer_dict = {}

    for num, seq in enumerate(list_of_seqs):
        for n in range(len(seq) - (ksize - 1)):
            slice = seq[n:n + ksize]

            if slice not in kmer_dict:
                counter = 1
                kmer_dict[slice] = counter
            else:
                kmer_dict[slice] += 1

    return kmer_dict


def count_occurrence(kmer_container):

    occur = dict(Counter(kmer_container.values()))

    return occur


def calculate_minimas(occurrence_dict_in):

    occur_df = pd.DataFrame(list(occurrence_dict_in.items()), columns=['Coverage', 'Count'])
    occur_sorted_df = occur_df.sort_values('Coverage', ascending=True)

    # Convert values to Numpy arrays because find_peaks only takes in arrays
    x_array = occur_sorted_df["Coverage"].to_numpy()
    y_array = occur_sorted_df["Count"].to_numpy()

    # Find the minimas
    y2 = y_array * -1  # Mirror the y values over the horizontal axis because the find_peaks only finds local maximas
    minima = find_peaks(y2)
    min_pos = x_array[minima[0]]  # List containing the positions of the minima
    min_height = y2[minima[0]]  # List containing the height of the mirrored minima
    real_min_height = min_height * -1  # Need to flip back over the horizontal axis to plot

    return x_array, y_array, min_pos, real_min_height, occur_sorted_df


def find_threshold(min_position, mirrored_min_height, occur_pd_sorted):

    # X and Y Values of the minima in a separate dataframe
    minimas_df = pd.DataFrame({'Coverage': min_position, 'Count': mirrored_min_height})

    # Finding first minima in main dataframe
    # Merge to find the indices of the minimas in the occur_df
    merged_df = pd.merge(occur_pd_sorted, minimas_df, on=['Coverage', 'Count'], how='left', indicator='Exist')
    merged_df['Exist'] = np.where(merged_df.Exist == 'both', True, False)  # Getting the indices of the minimas
    first_minima_x = int(minimas_df.loc[0][0])
    first_minima_y = int(minimas_df.loc[0][1])
    threshold_val = "The threshold value (local minimum) is at the point: {}, {}".format(first_minima_x, first_minima_y)

    return threshold_val, merged_df, first_minima_x


def find_remaining_kmers(df_merged):

    # Calculate how many kmers are above the threshold minima
    instance_index = df_merged.query('Exist == True').index.tolist()
    greater_coverage = len(df_merged) - (int(instance_index[0] + 1))
    remaining_kmers = "There are {} k-mers that have greater coverage than the threshold value.".format(greater_coverage)

    return remaining_kmers


def generate_plot(x_values, y_values, minima_x, minima_y, threshold_xvalue, threshold_message, kmers_remains):

    # Plotting the function
    fig = plt.figure()
    ax = fig.subplots()
    ax.plot(x_values, y_values)

    # Plotting the threshold
    ax.axvline(x=threshold_xvalue, color='r', label='Threshold', linestyle='--')

    # Plotting the local minimas
    ax.scatter(minima_x, minima_y, color='g', s=15, marker='X', label='Minima')
    ax.legend()
    ax.grid()
    plt.title("K-mer Coverage vs. Count")
    plt.xlabel("K-mer Coverage")
    plt.ylabel("Count")
    plt.subplots_adjust(bottom=0.2)
    plt.figtext(0.5, 0.01, threshold_message, ha='center', fontsize=12)
    plt.figtext(0.5, 0.05, kmers_remains, ha='center', fontsize=12)
    plt.show()


def get_fh(file_in, r_w_mode):
    """
    This function opens files and returns a file handle name.
    :param file_in: File to be parsed
    :param r_w_mode: Read or write mode
    :return: File handle name
    """
    try:
        if file_in.endswith('.gz'):
            with gzip.open(file_in, 'rb') as gzipped_f:
                fhandle_open = gzip.GzipFile(fileobj=gzipped_f)

        else:
            with open(file_in, r_w_mode) as fh:
                fhandle_open = fh.read().rstrip()
        return fhandle_open

    except IOError:
        print("Cannot open the file: ", file_in)
        sys.exit(1)
    except ValueError:
        print("Invalid mode entered!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate and plot k-mer occurrences.")
    parser.add_argument('-i', '--infile', dest='infile', help='Name of the file to open', required=True)

    args = parser.parse_args()
    main(args.infile)
