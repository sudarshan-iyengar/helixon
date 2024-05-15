import numpy as np
import pandas as pd
#from ot.partial import partial_wasserstein
from POT import POT  # Assuming the VSOT class is in VSOT.py


# Function to split positive and negative pressures
def split_pos_neg(p):
    p_pos_full = np.zeros_like(p, dtype=np.float16)
    p_neg_full = np.zeros_like(p, dtype=np.float16)

    p_pos_indices = np.where(p > 0)[0]
    p_neg_indices = np.where(p < 0)[0]

    p_pos_full[p_pos_indices] = p[p_pos_indices]
    p_neg_full[p_neg_indices] = p[p_neg_indices]

    p_neg_full = np.abs(p_neg_full)

    return p_pos_full, p_neg_full


# Load pressure and DOA data for locations 2 and 4
pressure_files = ['csv/pressure2.csv', 'csv/pressure4.csv']
doa_files = ['csv/doa2.csv', 'csv/doa4.csv']

p_w_pos = []
p_w_neg = []
doa_data = []

for i in range(len(pressure_files)):
    df_pressure = pd.read_csv(pressure_files[i], header=None).values.flatten().astype(np.float16)
    df_pressure = df_pressure[0:40000].astype(np.float16)
    df_doa = pd.read_csv(doa_files[i], header=None).values.flatten().astype(np.float16)
    df_doa = df_doa[0:40000].astype(np.float16)

    p_pos, p_neg = split_pos_neg(df_pressure)
    p_w_pos.append(p_pos)
    p_w_neg.append(p_neg)
    doa_data.append(df_doa)

# Create point clouds for positive and negative pressures
PC2_pos = {'pos': doa_data[0][:, None], 'mass': p_w_pos[0], 'n': len(doa_data[0])}
PC4_pos = {'pos': doa_data[1][:, None], 'mass': p_w_pos[1], 'n': len(doa_data[1])}

PC2_neg = {'pos': doa_data[0][:, None], 'mass': p_w_neg[0], 'n': len(doa_data[0])}
PC4_neg = {'pos': doa_data[1][:, None], 'mass': p_w_neg[1], 'n': len(doa_data[1])}

# Set expected distance and tolerance
distEx = np.mean(np.abs(PC4_pos['pos'] - PC2_pos['pos'])).astype(np.float16)  # IDK WHAT TO PUT SO I DID THIS
distTol = np.float16(0.2)

### --- Optimal Transport for Positive Pressures --- ###

# Calculate cost matrix
#C_pos = dist(PC2_pos['pos'], PC4_pos['pos'], metric='sqeuclidean').astype(np.float16)
# Instantiate POT object and perform partial OT
pot_pos = POT(PC2_pos, PC4_pos, distEx, distTol)
# Optimal sRatio is automatically calculated within the POT object

# --- Optimal Transport for Negative Pressures ---
# Calculate cost matrix
#C_neg = dist(PC2_pos['pos'], PC4_pos['pos'], metric='sqeuclidean').astype(np.float16)

# Instantiate POT object and perform partial OT
pot_neg = POT(PC2_neg, PC4_neg, distEx, distTol)

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
