import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import os

from read_csv import get_SRIR, read_SRIR, order_SRIR, get_error, rms_standard, combinePosNeg, make_negative

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
        sum_of_x += arr[0] * arr[1]
        sum_of_y += arr[0] * arr[2]
        sum_of_z += arr[0] * arr[3]

    # sum_of_x = sum_of_x/sum_of_p

    p_avg = sum_of_p / len(arrs)
    return [p_avg, sum_of_x, sum_of_y, sum_of_z]

# def get_ground_truth(minPos, maxPos):
#     ground_truth=[]
#     for pos in range(minPos, maxPos):
#         ground_truth += [read_SRIR("ground", pos, "0")]
#     return ground_truth



# def get_all_SRIR(minPos, maxPos):
#     ground_arr = get_ground_truth(minPos, maxPos)
#
#     for technique in techniqueArr:
#         for pos in range(minPos, maxPos):
#
#             ground = ground_arr[pos]
#
#             interp_loc_arr = []
#             # figure out interpolation locations and add to array
#
#             for interp_loc in interp_loc_arr:
#                 interp = read_SRIR(technique, pos, interp_loc)
#                 error = rms_standard(interp, ground)
#
#     return -1


def get_ground_truth(minPos, maxPos):
    ground_truths = []
    for pos in range(minPos, maxPos+1):
        ground_truths += [read_SRIR("ground", pos, "0")]
    return ground_truths


# pot_srir = get_SRIR()
# true_srir = order_SRIR(read_SRIR('pressure12.csv', 'doa12.csv'))

# print(len(true_srir))
# for k in true_srir:
#     print(np.linalg.norm(k[-3:]))

# for i in range(0, len(true_srir)-1):
#     print(np.linalg.norm(true_srir[i+1][-3:])-np.linalg.norm(true_srir[i][-3:]))


# print(len(true_srir))

#
# # Concatenate DataFrames into one long array
# combined_df = pd.concat(srir, ignore_index=True)
#
# # Convert DataFrame to array
# combined_array = combined_df.to_numpy()


#combinePosNeg ("pot_2_1-3_p_pos.csv", "pot_2_1-3_p_neg.csv") #this is ok :)
# print(read_SRIR("pot", "12", "11-13")[:3])

# test = get_error(true_srir, pot_srir, 40000)
# snipped = true_srir[:40000]
# #print(len(test))
# #print(len(true_srir[:40000]))
# print (max(test[0]))
# print (max(snipped[0]))
# print (rms_standard(snipped, test))


# a1 = [1, 2,4,2]
# a2 = [2, 2,2,0]
# a3 = [1, 2,-2,0]


make_negative("potMue003_3_1-5_p_neg.csv")
combinePosNeg("potMue003_3_1-5_p_pos.csv", "potMue003_3_1-5_p_neg.csv")
combinePosNeg("potMue003_3_1-5_doa_pos.csv", "potMue003_3_1-5_doa_neg.csv")

# srir = read_SRIR("pot" , "2" , "1-3")

# print(len(srir))
# print(srir)
# read all ground truths and put in array
# TODO: IMPLEMENT READ ALL GROUND TRUTHS -> STORE IN ground_truth
ground_truth_all = get_ground_truth(1,28)
# print(ground_truth_all[5])
# print(len(ground_truth_all))


evaluation_techniques = ["lin", "pot", "potMue003"]#, "pot"]  # ADD NEW AS NEEDED
NUM_POSITIONS = 28
WINDOW_SIZE = 10000

# for loop iterating over different evaluation techniques


# THIS LOOP WE EXECUTE TWICE. WE DO IT ONCE FOR MIN-MAX RESOLUTION AND THEN ONCE FOR ALL THE CENTRES

# EVALUATE MIN TO MAX -> each time the center value is evaluated
# for this we evaluate the following locations:
# 2_1-3, 3_1-5, 4_1-7, ... , 15_1-30
srirs_to_evaluate = []
for i in range(2, 16):
    srirs_to_evaluate.append(f"{i}_1-{2 * i - 1}")

# for now hardcode list:
srirs_to_evaluate = ['3_1-5']

# print(srirs_to_evaluate)


# for loop iterating over positions
for srir_location in srirs_to_evaluate:
    for technique in evaluation_techniques:


        position, interp_from = srir_location.split('_')

        srir = read_SRIR(technique, position, interp_from)
        ground_truth = ground_truth_all[int(position)]


        # COMMENT NEXT PART FOR NO TIME ALIGNMENT
        max_len_int = 0
        for i in srir:
            if  np.linalg.norm(i[-3:]) > max_len_int:
                max_len_int = np.linalg.norm(i[-3:])
        #print("srir: ", max_len_int)

        max_len_ground = 0
        for i in ground_truth:
            if  np.linalg.norm(i[-3:]) > max_len_ground:
                max_len_ground = np.linalg.norm(i[-3:])
        #print("ground: ", max_len_ground)


        # index_to_remove = []
        # if max_len_int > max_len_ground:
        #     for i in range(0, len(srir)):
        #         if np.linalg.norm(srir[i][-3:]) > max_len_ground:
        #             index_to_remove += [i]
        #     for i in reversed(index_to_remove):
        #         del srir[i]
        # else:
        #     for i in range(0, len(ground_truth)):
        #         if np.linalg.norm(ground_truth[i][-3:]) > max_len_int:
        #             index_to_remove += [i]
        #     for i in reversed(index_to_remove):
        #         del ground_truth[i]


        # for i in range(0, min(len(ground_truth), len(srir))):
           # print(ground_truth[i], srir[i])
           # print(np.linalg.norm(ground_truth[i][-3:]), np.linalg.norm(srir[i][-3:]))
        s = 0
        for i in srir:
            s += abs(i[0])
        #print("average p: ", s/len(srir))
        error = get_error(ground_truth, srir, min(len(ground_truth), len(srir)))

        print(srir_location, " ", technique, ": ", error)

        # do something with the error here


