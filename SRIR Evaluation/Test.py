import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import os

from read_csv import get_SRIR, read_SRIR
'''
    Normalize the length of the distance vector for a given point in the SRIR
'''
def normalize_len(arr):
    if len(arr) < 4:
        return arr
    arr[-3:] = arr[-3:] / np.linalg.norm(arr[-3:])
    return arr

'''
    Plot x and y of a set of SRIR points
'''
def plot2D(*arrs):
    # Extract x, y coordinates
    points = np.array(arrs)
    x_coords = points[:, 1]
    y_coords = points[:, 2]
    # Create the plot
    plt.figure()
    # Plot the points
    plt.scatter(x_coords, y_coords, c='r', marker='o')
    # Set labels
    plt.xlabel('X')
    plt.ylabel('Y')
    # Show the plot
    plt.grid(True)
    plt.show()

'''
    Calculate average virtual source
'''
def cluster_arrs(*arrs):
    sum_of_p = 0

    sum_of_x = 0
    sum_of_y = 0
    sum_of_z = 0

    for arr in arrs:
        arr = normalize_len(arr)
        sum_of_p += arr[0]
        sum_of_x += arr[0]*arr[1]
        sum_of_y += arr[0]*arr[2]
        sum_of_z += arr[0]*arr[3]

    #sum_of_x = sum_of_x/sum_of_p

    p_avg = sum_of_p/len(arrs)
    return [p_avg, sum_of_x, sum_of_y, sum_of_z]

pot_srir = get_SRIR()
true_srir = read_SRIR('pressure12.csv', 'doa12.csv')
# print(len(true_srir))
# for k in true_srir:
#     print(np.linalg.norm(k[-3:]))
print(len(true_srir))
    #
    # # Concatenate DataFrames into one long array
    # combined_df = pd.concat(srir, ignore_index=True)
    #
    # # Convert DataFrame to array
    # combined_array = combined_df.to_numpy()









# a1 = [1, 2,4,2]
# a2 = [2, 2,2,0]
# a3 = [1, 2,-2,0]
