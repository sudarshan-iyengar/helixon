import glob
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import os
import re


from read_csv import read_SRIR, order_SRIR, get_error, combinePosNeg, make_negative, cluster_arrs

def get_ground_truth(minPos, maxPos):
    ground_truths = []
    for pos in range(minPos, maxPos+1):
        ground_truths += [read_SRIR("ground", pos, "0")]
    return ground_truths

def match_pos_neg(pattern):
    relative_path = 'csvs\\'
    if "SRIR Eval" not in os.getcwd():
        relative_path = 'SRIR Evaluation\\csvs\\'
    folder_path = os.path.join(os.getcwd(), relative_path)
    
    files = os.listdir(folder_path)

    file_pairs = {}


    for file in files:
        match = pattern.match(file)
        if match:
            signifier, polarity = match.groups()
            if signifier not in file_pairs:
                file_pairs[signifier] = {}
            file_pairs[signifier][polarity] = file

    for signifier, pair in file_pairs.items():
        pos_file = pair.get('pos')
        neg_file = pair.get('neg')
        
        if pos_file and neg_file:

            combinePosNeg(pos_file, neg_file)
        else:
            print(f"Missing pair for signifier {signifier}: {pair}")




relative_path = 'csvs\\'
if "SRIR Eval" not in os.getcwd():
    relative_path = 'SRIR Evaluation\\csvs\\'
folder_path = os.path.join(os.getcwd(), relative_path)
pattern = os.path.join(folder_path, '*_p_neg.csv')

matching_files = glob.glob(pattern)

for file_name in matching_files:
    make_negative(file_name)

patterns = [re.compile(r'^(.*)_p_(pos|neg)\.csv$'),re.compile(r'^(.*)_doa_(pos|neg)\.csv$')]

for pattern in patterns:
    match_pos_neg(pattern)




# print(len(srir))
# print(srir)
# read all ground truths and put in array
# TODO: IMPLEMENT READ ALL GROUND TRUTHS -> STORE IN ground_truth
ground_truth_all = get_ground_truth(1,28)
# print(groun 30000d_truth_all[5])
# print(len(ground_truth_all))


#evaluation_techniques = ["lin", "potMat", "potMatMu003", "potMatMu5"]#, "pot"]  # ADD NEW AS NEEDED
evaluation_techniques = ["lin", "potMat", "potMatMu0001", "potMatMu0003", "potMatMu0010",  "potMatMu0030", "potMatMu0100", "potMatMu0300", "potMatMu1000", "potMatMu3000", "potMatMu9000"]
NUM_POSITIONS = 28
WINDOW_SIZE = 5000
GROUND_CLUSTER_SIZE = 500

# for loop iterating over different evaluation techniques


# THIS LOOP WE EXECUTE TWICE. WE DO IT ONCE FOR MIN-MAX RESOLUTION AND THEN ONCE FOR ALL THE CENTRES

# EVALUATE MIN TO MAX -> each time the center value is evaluated
# for this we evaluate the following locations:
# 2_1-3, 3_1-5, 4_1-7, ... , 15_1-30
srirs_to_evaluate = []
for i in range(2, 16):
    srirs_to_evaluate.append(f"{i}_1-{2 * i - 1}")

# for now hardcode list:
srirs_to_evaluate = ['8_1-15']
#srirs_to_evaluate = ['3_1-5']

# print(srirs_to_evaluate)


# for loop iterating over positions
for srir_location in srirs_to_evaluate:

    position, interp_from = srir_location.split('_')
    ground_truth = ground_truth_all[int(position)]
    ground_truth = order_SRIR(ground_truth)

    ground_truth = ground_truth[:WINDOW_SIZE]

    ground_truth_cluster = []

    max_len_ground = 0
    for i in ground_truth:
        if np.linalg.norm(i[-3:]) > max_len_ground:
            max_len_ground = np.linalg.norm(i[-3:])

    # print("max len: ", max_len_ground)
    interval = max_len_ground/GROUND_CLUSTER_SIZE
    bottom = 0
    # print("interval: ", interval)


    for i in range(0, GROUND_CLUSTER_SIZE):
        arr = []
        for p in range(bottom, len(ground_truth)):
            j = ground_truth[p]

            if np.linalg.norm(j[-3:]) < (i+1)*interval:
                arr += [j]
            else:
                bottom = p
                break
            # print(arr)

        ground_truth_cluster += [cluster_arrs(arr)]
    #print (ground_truth_cluster)
    for technique in evaluation_techniques:

        srir = read_SRIR(technique, position, interp_from)
        error = get_error(ground_truth_cluster, srir, min(len(ground_truth_cluster), len(srir)))

        print(srir_location, " ", technique, ": ", error)

        # do something with the error here





        # COMMENT NEXT PART FOR NO TIME ALIGNMENT
        # max_len_int = 0
        # for i in srir:
        #     if  np.linalg.norm(i[-3:]) > max_len_int:
        #         max_len_int = np.linalg.norm(i[-3:])
        # #print("srir: ", max_len_int)
        #
        # max_len_ground = 0
        # for i in ground_truth:
        #     if  np.linalg.norm(i[-3:]) > max_len_ground:
        #         max_len_ground = np.linalg.norm(i[-3:])
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