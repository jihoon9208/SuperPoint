import torch
import numpy as np
import math
from torch_scatter import scatter_add, scatter_mean
import torch_geometric.transforms as T

def range_value (x,y,z):
    return sqrt(x*x + y*y + z*z)




def comput_geof(xyz, target, k_nn):
    nver = len(xyz)
    geof = np.zeros(shape=(4,0))
    target_data = target
    xyz_data = xyz
    for iver in range(0, nver):
        position = np.zeros(shape=(k_nn, 3))
        iedg = k_nn *iver
        ind_nei = 0
        tmp_xyz = np.zeros(shape=(0,3))
        tmp_xyz[0][0] = xyz_data[3 * i_ver]
        tmp_xyz[0][1] = xyz_data[3 * i_ver + 1]
        tmp_xyz[0][2] = xyz_data[3 * i_ver + 2]
        for inei in ragne(0, k_nn):
            ind_nei = target_data[iedg]
            xi_range = range_value(tmp_xyz[0][0],tmp_xyz[0][1], tmp_xyz[0][2])
            xk_range = range_value(xyz_data[3 * ind_nei],xyz_data[3 * ind_nei+1] ,xyz_data[3 * ind_nei+2])
            weight = math.exp(-0.2*np.fabs(xk_range-xi_range))

            position[inei][0] = weight*(xyz_data[3 * ind_nei] - tmp_xyz[i_ver])
            position[inei][1] = weight*(xyz_data[3 * ind_nei + 1] - tmp_xyz[i_ver + 1])
            position[inei][2] = weight*(xyz_data[3 * ind_nei + 2] - tmp_xyz[i_ver + 2])
            
            iedg +=1
            inei +=1

        position_T = np.transpose(position)
        norm = 
        argmin_position = 
        for 
        iver +=1