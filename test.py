import torch
import numpy as np
import math
from torch_scatter import scatter_add, scatter_mean
import torch_geometric.transforms as T

def range_value (x,y,z):
    return sqrt(x*x + y*y + z*z)

<<<<<<< HEAD
def point_position(xyz_data, target_data , k_nn) :
    
    position = np.zeros(shape=(k_nn, 3))

    for i in range(0, k_nn):
        ind_nei = target_data[iedg]

        position[i][0] = xyz_data[3 * ind_nei]
        position[i][1] = xyz_data[3 * ind_nei + 1]
        position[i][2] = xyz_data[3 * ind_nei + 2]

        iedg +=1

=======
    return torch.sqrt(x*x + y*y + z*z)

def point_position(xyz, target , k_nn, iver):
    target_data = target
    xyz_data = xyz
    iedg = k_nn * iver
    iedg_data = iedg + 44
    ind_nei = 0
    position = torch.zeros(3, dtype=float)
    ind_nei = target_data[iedg_data]
    position[0] = xyz_data[ind_nei][0]
    position[1] = xyz_data[ind_nei][1]
    position[2] = xyz_data[ind_nei][2]    
  
>>>>>>> 8711a7a61ae817f77eb5d098569d42aa18dec902
    return position 

def neighboor_position(xyz, target, k_nn, iver):
    target_data = target
    xyz_data = xyz
    iedg = k_nn * iver
    iedg_data = iedg
    position = torch.zeros(4,3, dtype=float)
    for i in range(0,4):
        ind_nei = target_data[iedg_data]
        position[i][0] = xyz_data[ind_nei][0]
        position[i][1] = xyz_data[ind_nei][1]
        position[i][2] = xyz_data[ind_nei][2]
    
        iedg_data += 1
    return position

def compute_weight(xi_range, xk_range):
<<<<<<< HEAD
    return math.exp(-0.2*np.fabs(xk_range-xi_range))
=======
    return torch.exp(-0.2*torch.abs(xk_range-xi_range))
>>>>>>> 8711a7a61ae817f77eb5d098569d42aa18dec902

def comput_geof(xyz, target, k_nn):
    nver = len(xyz)
<<<<<<< HEAD
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
=======
    geof = np.zeros(shape=(nver, 4), dtype=float)
    xyz_data = torch.from_numpy(xyz)
    target = target.astype('int64')
    target_data = torch.from_numpy(target)

    sver = 0
    for iver in range(0, nver):
        position = torch.zeros(k_nn, 3, dtype=float)
        end_position = torch.zeros(3,dtype=float)
        iedg = k_nn * iver
        ind_nei = 0
        tmp_xyz = torch.zeros(3, dtype=float)
        tmp_xyz[0] = xyz_data[iver][0]
        tmp_xyz[1] = xyz_data[iver][1]
        tmp_xyz[2] = xyz_data[iver][2]
        end_position = point_position(xyz_data, target_data, k_nn, iver)
        
        # each point range value
        xi_range = range_value(tmp_xyz[0], tmp_xyz[1], tmp_xyz[2])
        xk_range = range_value(end_position[0],end_position[1],end_position[2] )

>>>>>>> 8711a7a61ae817f77eb5d098569d42aa18dec902
        weightk = compute_weight(xi_range, xk_range)
        for inei in ragne(0, k_nn):
            ind_nei = target_data[iedg]
<<<<<<< HEAD

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
=======
        
            position[inei][0] = weightk*(xyz_data[ind_nei][0] - tmp_xyz[0])
            position[inei][1] = weightk*(xyz_data[ind_nei][1] - tmp_xyz[1])
            position[inei][2] = weightk*(xyz_data[ind_nei][2] - tmp_xyz[2])
            
            iedg +=1
            inei +=1
        neighboor_point = torch.zeros(4,3, dtype=float)
        neighboor_point = neighboor_position(xyz_data, target_data, k_nn, iver)
        four_normal = 0
        for i in range(0, 4):
            xj_range = range_value(neighboor_point[i][0],neighboor_point[i][1] ,neighboor_point[i][2])
            weightj = compute_weight(xi_range, xj_range )
            end_neighboor = weightk * ( position[k_nn-1] - tmp_xyz )
            next_neghboor = weightj * ( neighboor_point[i] - tmp_xyz)
            #end_neighboor = np.transpose(end_neighboor)
            four_normal += torch.dot(end_neighboor,next_neghboor)
            
        position_norm = position*four_normal
        position_norm = position_norm.numpy()
        # compute the eigen values and vectors
        center_point = position_norm - position_norm.mean(axis=0)
        cp_cov= np.dot(center_point.T, center_point) /44
        w, v = np.linalg.eig(cp_cov)

        eginlambda = torch.from_numpy(w)
        v = torch.from_numpy(v)
        v1 = torch.FloatTensor([v[0][0], v[0][1],v[0][2]])
        v2 = torch.FloatTensor([v[1][0], v[1][1],v[1][2]])
        v3 = torch.FloatTensor([v[2][0], v[2][1],v[2][2]])

        # compute the dimensionality features
        linearity = (torch.sqrt(eginlambda[0]) - torch.sqrt(eginlambda[1])) / torch.sqrt(eginlambda[0]) 
        planarity = (torch.sqrt(eginlambda[1]) - torch.sqrt(eginlambda[2])) / torch.sqrt(eginlambda[0]) 
        scattering = torch.sqrt(eginlambda[2]) / torch.sqrt(eginlambda[0])

        unary_vector = torch.FloatTensor([
            [eginlambda[0] * torch.abs(v1[0]) + eginlambda[1] * torch.abs(v2[0]) + eginlambda[2] * torch.abs(v3[0])],
            [eginlambda[0] * torch.abs(v1[1]) + eginlambda[1] * torch.abs(v2[1]) + eginlambda[2] * torch.abs(v3[1])],
            [eginlambda[0] * torch.abs(v1[2]) + eginlambda[1] * torch.abs(v2[2]) + eginlambda[2] * torch.abs(v3[2])]
            ])

        norm = torch.sqrt(unary_vector[0]*unary_vector[0] + unary_vector[1]*unary_vector[1] 
                            + unary_vector[2]*unary_vector[2])    

        verticality = unary_vector[2] / norm
>>>>>>> 8711a7a61ae817f77eb5d098569d42aa18dec902


<<<<<<< HEAD
        norm =
        position_T = np.transpose(position)
         
        argmin_position = 
        for 
        iver +=1
=======
        sver +=1
        if (sver % 10000 == 0):
            print("%0.f%% done" % np.ceil(sver*100/nver))
        iver +=1
    return geof
>>>>>>> 8711a7a61ae817f77eb5d098569d42aa18dec902
