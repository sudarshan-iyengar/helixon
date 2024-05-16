import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#from ot.partial import partial_wasserstein
from POT import POT  # Assuming the VSOT class is in VSOT.py


# Function to split positive and negative pressures
def split_pos_neg(p):
    p_pos_full = np.zeros_like(p)#, dtype=np.float16)
    p_neg_full = np.zeros_like(p)#, dtype=np.float16)

    p_pos_indices = np.where(p > 0)[0]
    p_neg_indices = np.where(p < 0)[0]

    p_pos_full[p_pos_indices] = p[p_pos_indices]
    p_neg_full[p_neg_indices] = p[p_neg_indices]

    p_neg_full = np.abs(p_neg_full)

    return p_pos_full, p_neg_full


# Load pressure and DOA data for locations 2 and 4
pressure_files = ['helixon/csv/P_1.csv', 'helixon/csv/P_3.csv']
doa_files = ['helixon/csv/doa_1.csv', 'helixon/csv/doa_3.csv']

p_w_pos = []
p_w_neg = []
doa_data = []

for i in range(len(pressure_files)):
    df_pressure = pd.read_csv(pressure_files[i], header=None).values.flatten()#.astype(np.float16)
    df_pressure = df_pressure[400:30000]#.astype(np.float16)
    print("PRESSURE SHAPE: ")
    print(df_pressure.shape)
    print(df_pressure[0] )
    print()


    df_doa = pd.read_csv(doa_files[i], header=None).values#.astype(np.float16)
    df_doa = df_doa[400:30000]#.astype(np.float16)
    print("DOA SHAPE: ")
    print(df_doa.shape)
    print(df_doa[0])
    print()

    combined_df = pd.DataFrame({'pressure': df_pressure, 'doa_x': df_doa[:,0], 'doa_y': df_doa[:,1], 'doa_z': df_doa[:,2]})

    # Drop rows with any NaN values and get corresponding indices
    combined_df = combined_df.dropna()

    # Separate pressure data back
    df_pressure = combined_df['pressure'].values
    df_doa = combined_df[['doa_x', 'doa_y', 'doa_z']].values

    # Splitting positive and negative pressure
    p_pos, p_neg = split_pos_neg(df_pressure)
    p_w_pos.append(p_pos)
    p_w_neg.append(p_neg)
    doa_data.append(df_doa)

doa_data[0] = doa_data[0][0:20000]
doa_data[1] = doa_data[1][0:20000]
p_w_pos[0] = p_w_pos[0][0:20000]
p_w_pos[1] = p_w_pos[1][0:20000]
p_w_neg[0] = p_w_neg[0][0:20000]
p_w_neg[1] = p_w_neg[1][0:20000]

# Create point clouds for positive and negative pressures
PC2_pos = {'pos': doa_data[0].reshape(-1,3), 'mass': p_w_pos[0], 'n': len(doa_data[0])}
PC4_pos = {'pos': doa_data[1].reshape(-1,3), 'mass': p_w_pos[1], 'n': len(doa_data[1])}

PC2_neg = {'pos': doa_data[0].reshape(-1,3), 'mass': p_w_neg[0], 'n': len(doa_data[0])}
PC4_neg = {'pos': doa_data[1].reshape(-1,3), 'mass': p_w_neg[1], 'n': len(doa_data[1])}
print("CHECK IF DIMENSIONS ARE RIGHT: ")
print(PC2_neg['pos'].shape)
print(PC2_pos['mass'].shape)
print(PC2_pos['n'])
print()

# Set expected distance and tolerance
distEx = np.nanmean(np.abs(PC4_pos['pos'] - PC2_pos['pos']))#.astype(np.float16)  # IDK WHAT TO PUT SO I DID THIS
distTol = 0.2 #np.float16(0.2)
print("DIST EX: ")
print(distEx)
print()





### --- Optimal Transport for Positive Pressures --- ###

# Calculate cost matrix
# C_pos = dist(PC2_pos['pos'], PC4_pos['pos'], metric='sqeuclidean').astype(np.float16)
# Instantiate POT object and perform partial OT
pot_pos = POT(PC2_pos, PC4_pos, distEx, distTol)
print("REAACHEEEDDDDDDDDDDDDD end of pot pos")
print(np.mat(pot_pos.T).shape)
print()

# Optimal sRatio is automatically calculated within the POT object

# --- Optimal Transport for Negative Pressures ---
# Calculate cost matrix (put in class code better)
#C_neg = dist(PC2_pos['pos'], PC4_pos['pos'], metric='sqeuclidean').astype(np.float16)

# Instantiate POT object and perform partial OT
pot_neg = POT(PC2_neg, PC4_neg, distEx, distTol)
print("REAACHEEEDDDDDDDDDDDDD end of pot neg")

# --- Interpolation ---
# Choose interpolation parameter k
k = 0.5

# Get interpolated point clouds for positive and negative pressures
PCk_pos = pot_pos.interpPC(k)
PCk_neg = pot_neg.interpPC(k)

# --- You can now access the interpolated point clouds: PCk_pos, PCk_neg ---
print("Interpolated positive pressures:")
print(len(PCk_pos['mass']))

print("\nInterpolated negative pressures:")
print(len(PCk_neg['mass']))

