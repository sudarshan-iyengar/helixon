import numpy as np
from ot.partial import partial_wasserstein
from ot import dist
import torch
from geomloss import SamplesLoss


class POT:
    """
    Partial Optimal Transport

    Args:
        PC1 (dict): Point cloud 1 (with keys 'pos', 'mass', 'n')
        PC2 (dict): Point cloud 2 (with keys 'pos', 'mass', 'n')
        distEx (float): Expected distance (m) of transported image sources
        distTol (float): Tolerance for distance deviation from distEx, as a ratio
    """

    def __init__(self, PC1, PC2, distEx, distTol, blur=0.05, scaling=0.8, reach=None):

        # Copy input data
        self.PC1 = PC1
        self.PC2 = PC2
        self.distEx = distEx
        self.distTol = distTol
        self.T = 0

        # Prepare scost matrix
        C = dist(PC1['pos'], PC2['pos'], metric='sqeuclidean')#.astype(np.float16)
        print("COST MATRIX SHAPE: ")
        print(C.shape)



        # Determine dummy cost
        self.dumCost = (self.distEx * self.distTol) ** 2
        print("DUMMY COST: " + str(self.dumCost))

        # Define geomloss criterion
        self.loss = SamplesLoss("sinkhorn", p=2, blur=blur, scaling=scaling, reach=reach)


        # sRatio optimizer
        sVec = np.arange(0.01,1.01,0.01)
        sLen = len(sVec)
        costCoarse = np.zeros(sLen)

        for sInd in range(sLen):
            # Determine max amount of mass (s)
            maxS = min(np.sum(np.abs(PC1['mass'])), np.sum(np.abs(PC2['mass'])))
            s = sVec[sInd] * maxS
            print("mass transported: "+ str(s) + " index: " + str(sInd))

            # Optimization using partial_wasserstein
            Tx = self.loss(PC1['mass'], PC2['mass'], C, m=s)
            print("DONE WITH OBTAINING TX ONCE")

            if torch.isnan(Tx).any():
                print("Warning: Transport matrix contains NaN elements.")

            print("Number of non-zeros in Tx: " + str(np.count_nonzero(Tx)))
            costCoarse[sInd] = self.loss.cost.detach().item()
            print("Cost coarse at index " + str(sInd))
            
            if (sInd > 0) and (costCoarse[sInd] < costCoarse[sInd - 1]):
                # Stopping condition, don't update optimal values
                print(f"Found optimal sRatio = {self.sRatOpt}")
                self.Tx = Tx
                self.T = Tx
                self.sRatOpt=sVec[sInd]
                self.sOpt=self.sRatOpt * maxS
                break
            else:
                # Improved, save data
                self.Tx = Tx
                self.T = Tx
                self.sRatOpt = sVec[sInd]
                self.sOpt = self.sRatOpt * maxS

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
