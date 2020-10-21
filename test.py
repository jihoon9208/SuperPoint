import torch
import numpy as np

def range_value (x,y,z):

    return np.sqrt(x*x + y*y + z*z)

def point_position(xyz, target , k_nn):
    nver = len(xyz)
    target_data = target
    xyz_data = xyz
    for iver in range(0, nver):
        position = torch.zeros(k_nn, 3)
        iedg = k_nn *iver
        for inei in range(0, k_nn):
            ind_nei = target_data[iedg]
            print(ind_nei)
            position[inei][0] = xyz_data[3*ind_nei]
            position[inei][1] = xyz_data[3*ind_nei + 1]
            position[inei][2] = xyz_data[3*ind_nei + 2]
            
            iedg += 1
            inei += 1
    print(position)
    return position 

def compute_weight(xi_range, xk_range):
    return np.exp(-0.2*np.fabs(xk_range-xi_range))

def compute_geof(xyz, target, k_nn):
    nver = len(xyz)
    geof = torch.zeros(nver, 4, dtype=float)
    target = np.int64(target)
    target_data = torch.from_numpy(target)
    np.add(target,1 , out=target)
    xyz_data = torch.from_numpy(xyz)
    np.add(xyz, 1, out=xyz)
    
    for iver in range(0, nver):
        position = torch.zeros(k_nn, 3)
        iedg = k_nn *iver
        ind_nei = 0
        tmp_xyz = torch.zeros(3)
        tmp_xyz[0] = xyz_data[3*iver][0]
        tmp_xyz[1] = xyz_data[3*iver][1]
        tmp_xyz[2] = xyz_data[3*iver][2]

        position = point_position(xyz_data, target_data, k_nn)

        # each point range value
        xi_range = range_value(tmp_xyz[0].item(), tmp_xyz[1].item(), tmp_xyz[2].item())
        xk_range = range_value(position[k_nn-1][0],position[k_nn-1][1], position[k_nn-1][2] )
        
        weightk = compute_weight(xi_range, xk_range)
        for inei in range(0, k_nn):
            ind_nei = target_data[iedg]
            
            position[inei][0] = weightk*(xyz_data[3 * ind_nei][0] - tmp_xyz[iver][0])
            position[inei][1] = weightk*(xyz_data[3 * ind_nei][1] - tmp_xyz[iver][1])
            position[inei][2] = weightk*(xyz_data[3 * ind_nei][2] - tmp_xyz[iver][2])
            
            iedg +=1
            inei +=1
        
        neighboor_point = torch.zeros(4,3)
        for i in range(0, 3):
            neighboor_point[i][0] = position[i][0]
            neighboor_point[i][1] = position[i][1]
            neighboor_point[i][2] = position[i][2]
        four_normal = 0
        for i in range(0, 3):
            xj_range = range_value(neighboor_point[i][0],neighboor_point[i][1] ,neighboor_point[i][2])
            weightj = compute_weight(xi_range, xj_range )
            end_neighboor = weightk * ( position[k_nn-1] - tmp_xyz )
            next_neghboor = weightj * ( neighboor_point[i] - tmp_xyz)
            #end_neighboor = np.transpose(end_neighboor)
            four_normal += torch.dot(end_neighboor,next_neghboor)
            

        position_norm = position*four_normal

        # compute the eigen values and vectors
        centred_position = np.transpose(position_norm) - position_norm.mean(axis = 0)
        cov = torch.dot(centred_position.T, centred_position) / 45
        w, v = np.linalg.eig(cov)

        eginlambda = w
        v1 = torch.FloatTensor([v[0][0], v[0][1],v[0][2]])
        v2 = torch.FloatTensor([v[1][0], v[1][1],v[1][2]])
        v3 = torch.FloatTensor([v[2][0], v[2][1],v[2][2]])

        # compute the dimensionality features
        linearity = (np.sqrt(eginlambda[0]) - np.sqrt(eginlambda[1])) / np.sqrt(eginlambda[0]) 
        planarity = (np.sqrt(eginlambda[1]) - np.sqrt(eginlambda[2])) / np.sqrt(eginlambda[0]) 
        scattering = np.sqrt(eginlambda[2]) / np.sqrt(eginlambda[0])

        unary_vector = torch.FloatTensor([
            [eginlambda[0] * torch.abs(v1[0]) + eginlambda[1] * torch.abs(v2[0]) + eginlambda[2] * torch.abs(v3[0])],
            [eginlambda[0] * torch.abs(v1[1]) + eginlambda[1] * torch.abs(v2[1]) + eginlambda[2] * torch.abs(v3[1])],
            [eginlambda[0] * torch.abs(v1[2]) + eginlambda[1] * torch.abs(v2[2]) + eginlambda[2] * torch.abs(v3[2])]
            ])

        norm = torch.sqrt(unary_vector[0]*unary_vector[0] + unary_vector[1]*unary_vector[1] 
                            + unary_vector[2]*unary_vector[2])    

        verticality = unary_vector[2] / norm

        geof[iver][0] = linearity
        geof[iver][1] = planarity
        geof[iver][2] = scattering
        geof[iver][3] = verticality

        sver +=1
        if (sver % 10000 == 0):
            print("%0.f%% done" % sver)
            print("%0.f%% done" % torch.ceil(sver*100/nver))
        iver +=1
    return 