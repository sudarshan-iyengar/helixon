import numpy as np
import pandas as pd
import csv
import os


def get_SRIR():
    srir_pos = read_SRIR('11_13_pos12_p_pos.csv', '11_13_pos12_pos.csv')
    srir_neg = read_SRIR('11_13_pos12_p_neg.csv', '11_13_pos12_neg.csv')
    return order_SRIR(srir_pos + srir_neg)


def read_SRIR(file_p, file_doa):
    relative_path = 'SRIR Evaluation/csvs/' # 'csvs/'
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

'''
    Calculate average virtual source
'''
def cluster_arrs(arrs):
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


'''
    Normalize the length of the distance vector for a given point in the SRIR
'''
def normalize_len(arr):
    if len(arr) < 4:
        return arr
    arr[-3:] = arr[-3:] / np.linalg.norm(arr[-3:])
    return arr

def get_error(ground, interp, window_size):
    #trunctate ground to window_size
    virtual_interp = []
    ground =  ground[:window_size]
    
    for i in range (0,len(ground)):
        cur_len = np.linalg.norm(interp[i][-3:])
        upper = 0
        lower = 0
        #get the boundaries (in terms of length)
        if i != 0:
            lower = cur_len - (cur_len - np.linalg.norm(ground[i-1][-3:]))/2
        if i != len(ground)-1:
            upper = cur_len + (np.linalg.norm(ground[i+1][-3:]) - cur_len)/2

        
        #find all values in interp, within the boundary
        cluster = []
        for i in range(0, len(interp)):
            val = np.linalg.norm(interp[i][-3:])
            if val >= lower and val < upper:
                cluster += interp[i]
            
            if val > upper:
                break
        
                
        #we find the combined virtual source from these values -> cluster_arrs
        #calculate rms

        #at some point normalize all doa before we compare -> WORK WITH ANGLES INSTEAD
    for i in range(0, len(ground)):

