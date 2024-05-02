import pandas as pd
import csv
import os

def get_SRIR():
    srir_pos = []
    srir_neg = []
    relative_path = 'SRIR Evaluation/csvs/'
    folder_path = os.path.join(os.getcwd(), relative_path)

    # Iterate over pos CSV
    for filename in os.listdir(folder_path):
        file1_path = os.path.join(folder_path, "11_13_pos12_p_pos.csv")
        file2_path = os.path.join(folder_path, "11_13_pos12_pos.csv")

        # Open and read both CSV files in parallel
        with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
            reader1 = csv.reader(file1)
            reader2 = csv.reader(file2)

            # Iterate through rows in parallel
            for row1, row2 in zip(reader1, reader2):
                # Append data from each row to respective lists
                srir_pos += [row1 + row2]
    # Iterate over pos CSV
    for filename in os.listdir(folder_path):
        file1_path = os.path.join(folder_path, "11_13_pos12_p_neg.csv")
        file2_path = os.path.join(folder_path, "11_13_pos12_neg.csv")

        # Open and read both CSV files in parallel
        with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
            reader1 = csv.reader(file1)
            reader2 = csv.reader(file2)

            # Iterate through rows in parallel
            for row1, row2 in zip(reader1, reader2):
                # Append data from each row to respective lists
                srir_neg += [row1 + row2]

    # print("pos len: " + str(len(srir_pos)) )
    # print("neg len: " + str(len(srir_neg)) )
    return srir_pos + srir_neg