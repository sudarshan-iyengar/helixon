import numpy as np
import pandas as pd
import csv
import os


def get_SRIR():
    srir_pos = read_SRIR('11_13_pos12_p_pos.csv', '11_13_pos12_pos.csv')
    srir_neg = read_SRIR('11_13_pos12_p_neg.csv', '11_13_pos12_neg.csv')
    return order_SRIR(srir_pos + srir_neg)


def read_SRIR(file_p, file_doa):
    relative_path = 'csvs/'
    folder_path = os.path.join(os.getcwd(), relative_path)
    # Iterate over pos CSV

    file1_path = os.path.join(folder_path, file_p)
    file2_path = os.path.join(folder_path, file_doa)
    srir = []
    # Open and read both CSV files in parallel
    with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
        reader1 = csv.reader(file1)
        reader2 = csv.reader(file2)

        # Iterate through rows in parallel
        for row1, row2 in zip(reader1, reader2):
            # Append data from each row to respective lists
            srir += [[float(row1[0])] + [float(x) for x in row2]]
    srir = remove_zeros(srir)
    return srir


def order_SRIR(arr):
    lengths = []
    for val in arr:
        lengths += [np.linalg.norm(val[-3:])]

    # Combine the two arrays into pairs
    combined = list(zip(lengths, arr))

    # Sort the combined list based on the values from array1
    sorted_combined = sorted(combined, key=lambda x: x[0])

    # Extract the sorted arrays
    srir = [pair[1] for pair in sorted_combined]
    return srir


'''
    Remove all elements where norm of doa is zero
'''
def remove_zeros(srir):
    return [element for element in srir if np.linalg.norm(element[-3:]) > 0]
