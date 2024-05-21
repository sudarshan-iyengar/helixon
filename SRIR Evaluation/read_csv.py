import numpy as np
import pandas as pd
import csv
import os


def get_SRIR():
    srir_pos = read_SRIR('11_13_pos12_p_pos.csv', '11_13_pos12_pos.csv')
    srir_neg = read_SRIR('11_13_pos12_p_neg.csv', '11_13_pos12_neg.csv')
    return order_SRIR(srir_pos + srir_neg)

def read_SRIR(technique, position, interp_from):
    
    relative_path = 'srir\\' #SRIR Evaluation\\
    if "SRIR Eval" not in os.getcwd():
        relative_path = 'SRIR Evaluation\\srir\\'

    if interp_from == "0":
        relative_path = 'ground_truth\\' #SRIR Evaluation\\
        if "SRIR Eval" not in os.getcwd():
            relative_path = 'SRIR Evaluation\\ground_truth\\'
    folder_path = os.path.join(os.getcwd(), relative_path)

    file_name = f"{technique}_{position}_{interp_from}"

    file1_path = os.path.join(folder_path, f"{file_name}_p.csv")
    file2_path = os.path.join(folder_path, f"{file_name}_doa.csv")

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
    #order_SRIR(srir)? --> here or outside of function?
    return srir


def combinePosNeg(file_pos_name, file_neg_name):
    relative_path = 'csvs\\'
    if "SRIR Eval" not in os.getcwd():
        relative_path = 'SRIR Evaluation\\csvs\\'
    folder_path = os.path.join(os.getcwd(), relative_path)

    destination_path = 'srir\\'
    if "SRIR Eval" not in os.getcwd():
        destination_path = 'SRIR Evaluation\\srir\\'
    destination_folder_path = os.path.join(os.getcwd(), destination_path)

    filepos_path = os.path.join(folder_path, file_pos_name)
    fileneg_path = os.path.join(folder_path, file_neg_name)

    if not (file_neg_name.endswith("_neg.csv") and file_pos_name.endswith("_pos.csv")):
            raise ValueError("Input files must follow the naming convention *_neg.csv and *_pos.csv")

    common = file_pos_name.rsplit("_", 1)[0]

    file_combined_name = f"{common}.csv"
    file_combined_path = os.path.join(destination_folder_path, file_combined_name)

    try:

        with open(filepos_path, 'r', newline='') as file_pos, \
             open(fileneg_path, 'r', newline='') as file_neg, \
             open(file_combined_path, 'w', newline='') as file_combined:

            reader1 = csv.reader(file_pos)
            reader2 = csv.reader(file_neg)
            writer = csv.writer(file_combined)

            for row in reader1:
                writer.writerow(row)

            for row in reader2:
                writer.writerow(row)

        print(f"CSV files {filepos_path} and {fileneg_path} have been successfully combined into {file_combined_path}.")
        return file_combined_name

    except Exception as e:
        print(f"An error occurred: {e}")
        return None


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
    if len(arrs) == 0:
        return [0,0,0,0]
    sum_of_p = 0

    w_sum_of_x = 0
    w_sum_of_y = 0
    w_sum_of_z = 0

    sum_of_len = 0

    # we get the average pressure, and the direction weighed by the pressure
    for arr in arrs:
        sum_of_len += np.linalg.norm(arr[-3:])

        arr = normalize_len(arr)
        sum_of_p += arr[0]
        w_sum_of_x += arr[0]*arr[1]
        w_sum_of_y += arr[0]*arr[2]
        w_sum_of_z += arr[0]*arr[3]



    len_avg = sum_of_len/len(arrs)
    len_temp = np.linalg.norm([w_sum_of_x, w_sum_of_y, w_sum_of_z])
    normalized_x = w_sum_of_x/len_temp
    normalized_y = w_sum_of_y/len_temp
    normalized_z = w_sum_of_z/len_temp

    scaled_x = normalized_x * len_avg
    scaled_y = normalized_y * len_avg
    scaled_z = normalized_z * len_avg

    # print("average len: ", len_avg, "; scaled len: ", np.linalg.norm([scaled_x, scaled_y, scaled_z]))

    p_avg = sum_of_p/len(arrs)
    return [p_avg, scaled_x, scaled_y, scaled_z]


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
    interp = order_SRIR(interp)
    virtual_interp = []
    ground = ground[:window_size]
    bottom = 0
    for i in range(0, len(ground)):
        # if i % 5000 == 0:
        #     print(i)
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


        #print(bottom)
        updated_bottom = False
        for j in range(bottom, len(interp)):
            val = np.linalg.norm(interp[j][-3:])
            if val >= lower and val < upper:
                if not updated_bottom:
                    bottom = j
                    updated_bottom = True
                cluster += [interp[j]]

            if val > upper:
                break
        
            
        #we find the combined virtual source from these values -> cluster_arrs
        virtual_src = cluster_arrs(cluster)
        virtual_interp += [virtual_src]

    #print(len(virtual_interp))
    return rms_standard(ground, virtual_interp)
    #return virtual_interp
        #calculate rms
        
def rms_error(true_srir, interp_srir):

    for i in range(0, len(true_srir)):
        true_srir[i] = normalize_len(true_srir[i])

    count = 0
    for i in interp_srir:
        if i[0] == 0:
            count += 1


    return -1

        #at some point normalize all doa before we compare -> WORK WITH ANGLES INSTEAD

###########################

def rms_standard(true_srir, interp_srir):
    p_sqrd = []
    doa_sqrd = []

    normalize_len(true_srir)
    for i in range (0, len(true_srir)):

        #pressure error
        e1 = (true_srir[i][0] - interp_srir[i][0])**2
        e2 = np.sqrt((true_srir[i][1] - interp_srir[i][1])**2 + (true_srir[i][2] - interp_srir[i][2])**2 + (true_srir[i][3] - interp_srir[i][3])**2)

        p_sqrd += [e1]
        doa_sqrd += [e2]


    err_p = np.sqrt(np.sum(p_sqrd)/len(p_sqrd))
    err_doa = np.sum(doa_sqrd)/len(doa_sqrd)

    return [err_p, err_doa]

