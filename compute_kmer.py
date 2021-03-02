#! /usr/bin/env python3

"""
Created by Dana Rocha 2/24/21
"""

import sys
import argparse
#import gzip
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

def main(infile):

    fh_in = open(infile, "r")

    kmer_size = int(input("Please enter a k-mer size (31, 41, 51, 61, 71, 81, 91, or 101): "))
    checked = check_size(kmer_size)

    cleaned_seq = clean_seq(fh_in)

    kmers = calc_kmer(cleaned_seq, checked)

    kmer_occurs = count_occurrence(kmers)

    generate_plot(kmer_occurs)


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


def clean_seq(handle_in):
    seq_list = []

    for line in handle_in:
        seq_list.append(line.rstrip())

    return seq_list


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


def generate_plot(occur_dict):

    df = pd.DataFrame(list(occur_dict.items()), columns=['Coverage', 'Count'])
    df = df.sort_values('Coverage', ascending=True)

    x_vals = df["Coverage"].to_numpy()
    y_vals = df["Count"].to_numpy()

    # Find peaks (maximas)
    peaks = find_peaks(y_vals, height=2, threshold=1, distance=1)
    height = peaks[1]['peak_heights']  # List containing the height of the peaks
    peak_pos = x_vals[peaks[0]]  # List containing the positions of the peaks

    # Find the minimas
    y2 = y_vals * -1 # Mirror the y values over the horizontal axis because this module only finds peaks
    minima = find_peaks(y2)
    min_pos = x_vals[minima[0]]  # list containing the positions of the minimas
    min_height = y2[minima[0]]  # list containing the height of the mirrored minimai
    real_min_height = min_height * -1 # Need to flip back over the horizontal axis to plot

    # Plotting the function
    fig = plt.figure()
    ax = fig.subplots()
    ax.plot(x_vals, y_vals)

    # Plotting the local maximas
    # ax.scatter(peak_pos, height, color='r', s=15, marker='D', label='Maxima')

    # Plotting the local minimas
    ax.scatter(min_pos, real_min_height, color='g', s=15, marker='X', label='Minima')
    ax.legend()
    ax.grid()
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
            #data = gzip.decompress(response.read(file_in))
            #fhandle_open = data.decode("utf-8")
            #fhandle_open = gzip.open(text, r_w_mode)

            fhandle_open = gzip.open(file_in, r_w_mode)

            #with open(file_in, r_w_mode) as file_in:
            #    gzip_fh = gzip.GzipFile(fileobj=file_in)
            #    fhandle_open = gzip.open(gzip_fh)

            #fhandle_open = gzip.GzipFile(file_in, r_w_mode)

            #with codecs.open(file_in, r_w_mode, encoding='utf-8') as f_open:
            #   fhandle_open = f_open.readlines()

            #file_zipped = gzip.open(file_in, 'rb')
            #fhandle_open = file_zipped.read()
            #file_zipped.close()

            #return fhandle_open

        else:
            fhandle_open = open(file_in, r_w_mode)
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
