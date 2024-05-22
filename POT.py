import numpy as np
from ot.partial import partial_wasserstein
from ot import dist


class POT:
    """
    Partial Optimal Transport

    Args:
        PC1 (dict): Point cloud 1 (with keys 'pos', 'mass', 'n')
        PC2 (dict): Point cloud 2 (with keys 'pos', 'mass', 'n')
        distEx (float): Expected distance (m) of transported image sources
        distTol (float): Tolerance for distance deviation from distEx, as a ratio
    """


    def __init__(self, PC1, PC2, distEx, distTol):

        # Copy input data
        self.PC1 = PC1
        self.PC2 = PC2
        self.distEx = distEx
        self.distTol = distTol
        self.T = 0

        # Prepare scost matrix
        C = dist(PC1['pos'], PC2['pos'], metric='sqeuclidean')#.astype(np.float16)
        #C=self.calc_cost(PC1['pos'], PC2['pos'], PC1['mass'], PC2['mass'], 0)
        print("COST MATRIX SHAPE: ")
        print(C.shape)



        # Determine dummy cost
        self.dumCost = (self.distEx * self.distTol) ** 2
        print("DUMMY COST: " + str(self.dumCost))

        # sRatio optimizer
        sVec = np.arange(0.98,1.00,0.01)
        sLen = len(sVec)
        costCoarse = np.zeros(sLen)

        last_successful_Tx = None
        last_successful_sRatOpt = None
        last_successful_sOpt = None

        for sInd in range(sLen):
            try:
                # Determine max amount of mass (s)
                maxS = min(np.sum(np.abs(PC1['mass'])), np.sum(np.abs(PC2['mass'])))
                s = sVec[sInd] * maxS
                print("mass transported: " + str(s) + " index: " + str(sInd))

                # Optimization using partial_wasserstein
                Tx, log = partial_wasserstein(a=PC1['mass'], b=PC2['mass'], M=C, m=s, log=True)
                print("DONE WITH OBTAINING TX ONCE")

                if np.isnan(Tx).any():
                    print("Warning: Transport matrix contains NaN elements.")

                print("Number of non-zeros in Tx: " + str(np.count_nonzero(Tx)))
                costCoarse[sInd] = log['cost']
                print("Cost coarse at index " + str(sInd))

                if (sInd > 0) and (costCoarse[sInd] < costCoarse[sInd - 1]):
                    # Stopping condition, don't update optimal values
                    print(f"Found optimal sRatio = {self.sRatOpt}")
                    self.Tx = last_successful_Tx  # Revert to the last successful value
                    self.T = last_successful_Tx
                    self.sRatOpt = last_successful_sRatOpt
                    self.sOpt = last_successful_sOpt
                    break
                else:
                    # Improved, save data
                    last_successful_Tx = Tx  # Store the successful value
                    last_successful_sRatOpt = sVec[sInd]
                    last_successful_sOpt = sVec[sInd] * maxS

                    self.Tx = Tx
                    self.T = Tx
                    self.sRatOpt = sVec[sInd]
                    self.sOpt = sVec[sInd] * maxS

            except Exception as e:
                # Handle the error: print the error message and revert to the last successful values
                print(f"Error at index {sInd}: {e}")
                if last_successful_Tx is not None:
                    print("Reverting to last successful transport matrix and values.")
                    self.Tx = last_successful_Tx
                    self.T = last_successful_Tx
                    self.sRatOpt = last_successful_sRatOpt
                    self.sOpt = last_successful_sOpt
                break  # Optionally, you can break the loop or continue to try the next index


    def calc_cost(self, XA, XB, PA, PB, mue):
        """
        Compute the Euclidean distance between each pair of the two collections of inputs.

        Parameters
        ----------
        XA : array_like
            An m_A by n array of m_A original observations in an n-dimensional space.
        XB : array_like
            An m_B by n array of m_B original observations in an n-dimensional space.

        Returns
        -------
        C : ndarray
            An m_A by m_B distance matrix is returned. For each i and j, the Euclidean
            distance between XA[i] and XB[j] is computed and stored in the ij-th entry.
            Based on the absolute difference in pressure between the pressure values,
            a penalty is added only pressures of opposite signs.
        """
        C = dist(XA,XB, metric='sqeuclidean')
        sign_mask = (PA[:,np.newaxis]*PB[np.newaxis,:])<0
        pressure_diff = np.abs(PA[:,np.newaxis]-PB[np.newaxis,:])
        penalty = mue * pressure_diff * sign_mask
        C+=penalty

        return C    

    def interpPC(self, k):
        """
        Interpolates the partial OT problem at k

        Args:
            k (float): Interpolation parameter (between 0 and 1)

        Returns:
            dict: Interpolated point cloud (with keys 'pos', 'mass', 'n')
        """

        assert 0 <= k <= 1, "k must be between 0 and 1."

        # 1. Transported Mass (TM)
        Ttm = self.T

        # 2. Moving Mass (MM) - Not applicable in partial OT

        # Calculate transported & growing masses at interp value k
        Ttmk = Ttm
        I, J = np.where(Ttmk)
        TM = Ttmk[I, J]
        numTm = len(TM)
        k_array = np.array([k])

        # Calculate transported mass positions at interp value k
        posTm = (1 - k_array)[:, None] * self.PC1['pos'][I] + k_array[:, None] * self.PC2['pos'][J]

        # 3. Stationary Mass (SM) - Handled internally by partial_wasserstein

        # Concatenate masses and form point cloud struct
        PCk = {
            'pos': posTm,
            'mass': TM,
            'n': len(TM)
        }

        # Object remembers the last interpolation (for plotting)
        self.k = k
        self.PCk = PCk

        return PCk
    

    


    def interpolatePC(self, k):
        """
        Interpolates the partial OT problem at k, handling dummy nodes explicitly.

        Args:
            k (float): Interpolation parameter (between 0 and 1)

        Returns:
            dict: Interpolated point cloud (dictionary with 'pressure' and 'doa' keys)
        """

        assert 0 <= k <= 1, "k must be between 0 and 1."

        T_prime = np.copy(self.Tx)
        rk_list = []  # List to store interpolated pressures
        zk_list = []  # List to store interpolated DOA values

        # Handling Vanishing Points
        for i, ui in enumerate(T_prime.sum(axis=1)[:-1]):
            if ui > 0:
                if np.linalg.norm(T_prime[i, :-1], ord=1) == 0:
                    rk = (1 - k) * ui
                    zk = self.PC1['pos'][i]
                    rk_list.append(rk)
                    zk_list.append(zk)
                else:
                    T_prime[i, :-1] += (1 - k) * ui * T_prime[i, :-1] / np.linalg.norm(T_prime[i, :-1], ord=1)

        # Handling Appearing Points
        for j, vj in enumerate(T_prime.sum(axis=0)[:-1]):
            if vj > 0:
                if np.linalg.norm(T_prime[:-1, j], ord=1) == 0:
                    rk = k * vj
                    zk = self.PC2['pos'][j]
                    rk_list.append(rk)
                    zk_list.append(zk)
                else:
                    T_prime[:-1, j] += k * vj * T_prime[:-1, j] / np.linalg.norm(T_prime[:-1, j], ord=1)

        # Handling Moving Points
        for i, j in zip(*np.where(T_prime[:-1, :-1] > 0)):
            rk = T_prime[i, j]
            zk = (1 - k) * self.PC1['pos'][i] + k * self.PC2['pos'][j]
            rk_list.append(rk)
            zk_list.append(zk)

        # Create the dictionary with 'pressure' and 'doa' keys
        PCk = {
            'pressure': np.array(rk_list),
            'doa': np.array(zk_list)
        }

        return PCk


    