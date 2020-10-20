import torch
import numpy as np
import math
from torch_scatter import scatter_add, scatter_mean
import torch_geometric.transforms as T

def range_value (x,y,z):
    return sqrt(x*x + y*y + z*z)

def point_position(xyz_data, target_data , k_nn) :
    
    position = np.zeros(shape=(k_nn, 3))

    for i in range(0, k_nn):
        ind_nei = target_data[iedg]

        position[i][0] = xyz_data[3 * ind_nei]
        position[i][1] = xyz_data[3 * ind_nei + 1]
        position[i][2] = xyz_data[3 * ind_nei + 2]

        iedg +=1

    return position 

def compute_weight(xi_range, xk_range):
    return math.exp(-0.2*np.fabs(xk_range-xi_range))

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
        position = point_position(xyz_data, target_data, k_nn)

        xi_range = range_value(tmp_xyz[0][0],tmp_xyz[0][1], tmp_xyz[0][2])
        xk_range = range_value(position[k_nn-1][0],position[k_nn-1][1], position[k_nn-1][2] )
        weightk = compute_weight(xi_range, xk_range)
        for inei in ragne(0, k_nn):
            ind_nei = target_data[iedg]

            position[inei][0] = weightk*(xyz_data[3 * ind_nei] - tmp_xyz[i_ver])
            position[inei][1] = weightk*(xyz_data[3 * ind_nei + 1] - tmp_xyz[i_ver + 1])
            position[inei][2] = weightk*(xyz_data[3 * ind_nei + 2] - tmp_xyz[i_ver + 2])
            
            iedg +=1
            inei +=1

        neighboor_point = np.zeros(shape=(4,3))
        for i in range(0, 4):
            neighboor_point[i][0] = position[i][0]
            neighboor_point[i][1] = position[i][1]
            neighboor_point[i][2] = position[i][2]


        for i in range(0, 3):
            xj_range = range_value(neighboor_point[j][0],neighboor_point[j][1] ,neighboor_point[j][2])
            weightj = compute_weight(xi_range, xj_range )
            end_neighboor = weightk * ( neighboor_point[k_nn-1] - tmp_xyz[i] )
            next_neghboor = weightj * ( neighboor_point[i] - tmp_xyz[i] )
            end_neighboor = np.transpose(end_neighboor)
            four_normal += np.dot(end_neighboor,next_neghboor)


        norm =
        position_T = np.transpose(position)
         
        argmin_position = 
        for 
        iver +=1