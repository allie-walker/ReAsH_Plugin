from math import *
import numpy as np

ideal_dists = [];
#ideal_dists.append([0, 3.3, 5.2, 7.4])
#ideal_dists.append([3.3, 0, 7.4, 8.4])
#ideal_dists.append([5.2, 7.4, 0, 3.3])
#ideal_dists.append([7.4, 8.4, 3.3, 0])
ideal_dists.append([0, 3.6, 5.1, 7.2])
ideal_dists.append([3.6, 0, 7.2, 7.6])
ideal_dists.append([5.1, 7.2, 0, 3.6])
ideal_dists.append([7.2, 7.6, 3.6, 0])

distance_dict = {}

def dist(v1, v2):
    return sqrt((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1])+(v1[2]-v2[2])*(v1[2]-v2[2]))

def MatrixVecMult(matrix, translation_vec, vector):
    new_vector = []
    for i in range(0, len(matrix)):
        value= 0;
        for j in range(0, len(matrix[i])):
            value += vector[j]*matrix[i][j]
        new_vector.append(value)
    new_vector += translation_vec;
    #for i in range(0, len(new_vector)):
     #   new_vector[i] += translation_vec[i]
    return new_vector;

def MultMatrices(m1, m2):
    new_matrix = [];
    new_matrix.append([0.0,0.0,0.0])
    new_matrix.append([0.0,0.0,0.0])
    new_matrix.append([0.0,0.0,0.0])
    for i1 in range(0, len(m1)):
        for j1 in range(0, len(m1[i1])):
            for i2 in range(0, len(m2)):
                for j2 in range(0, len(m2[i2])):
                    new_matrix[i1][j2] += m1[i1][j1]*m2[i2][j2]
    return new_matrix;

def findInverse(matrix):
    det = 0.0;
    det += matrix[0][0] *(matrix[1][1]*matrix[2][2] - matrix[1][2]*matrix[2][1])
    det -= matrix[0][1] *(matrix[1][0]*matrix[2][2] - matrix[1][2]*matrix[2][0])
    det += matrix[0][2] *(matrix[1][0]*matrix[2][1] - matrix[1][1]*matrix[2][0])

    inverse = [];
    inverse.append([matrix[0][0]/det, matrix[0][1]/det, matrix[0][2]/det])
    inverse.append([matrix[1][0]/det, matrix[1][1]/det, matrix[1][2]/det])
    inverse.append([matrix[2][0]/det, matrix[2][1]/det, matrix[2][2]/det])

    return inverse;

def findDet(matrix):
    det = 0.0;
    det += matrix[0][0] *(matrix[1][1]*matrix[2][2] - matrix[1][2]*matrix[2][1])
    det -= matrix[0][1] *(matrix[1][0]*matrix[2][2] - matrix[1][2]*matrix[2][0])
    det += matrix[0][2] *(matrix[1][0]*matrix[2][1] - matrix[1][1]*matrix[2][0])
    return det;

def findSulfurCoords(coords):
    S1 = np.matrix.transpose(np.array([0.1430,    2.5514,    1.1983]))
    S2 = np.matrix.transpose(np.array([ 2.3520,    0.9126,    1.2657]))
    S3 = np.matrix.transpose(np.array([0.0700,    0.0645,    2.8205]))

    trans_vec =[-1* float(coords[2][0]), -1* float(coords[2][1]), -1* float(coords[2][2])]
    trans_coords = []
    trans_coords.append([0.0, 0.0, 0.0])
    trans_coords.append([0.0, 0.0, 0.0])
    trans_coords.append([0.0, 0.0, 0.0])
    for i in range(0, len(coords)):
        for j in range(0, len(coords)):
            trans_coords[i][j] = float(coords[i][j]) + trans_vec[j]

    CA_N_len = sqrt(trans_coords[0][0]*trans_coords[0][0] + trans_coords[0][1]*trans_coords[0][1] + trans_coords[0][2]*trans_coords[0][2])
    CA_C_len = sqrt(trans_coords[1][0]*trans_coords[1][0] + trans_coords[1][1]*trans_coords[1][1] + trans_coords[1][2]*trans_coords[1][2])
    trans_coords_np = np.array(trans_coords);
    #print trans_coords_np

    #rotate in y direction so that z=0
    if trans_coords[0][2] <= 0 and  trans_coords[0][0] > 0:
        theta_y = 270.0*pi/180.0 + acos(abs(trans_coords[0][2])/sqrt(trans_coords[0][0]*trans_coords[0][0]+ trans_coords[0][2]*trans_coords[0][2]));
    elif  trans_coords[0][2] > 0 and  trans_coords[0][0] > 0:
        theta_y = asin(abs(trans_coords[0][2])/sqrt(trans_coords[0][0]*trans_coords[0][0]+ trans_coords[0][2]*trans_coords[0][2]));
    elif  trans_coords[0][2] <= 0 and  trans_coords[0][0] <= 0:
        theta_y = -1*(90.0*pi/180.0 +acos(abs(trans_coords[0][2])/sqrt(trans_coords[0][0]*trans_coords[0][0]+ trans_coords[0][2]*trans_coords[0][2])));
    else: #z pos x neg
        theta_y = 90.0*pi/180.0 + acos(abs(trans_coords[0][2])/sqrt(trans_coords[0][0]*trans_coords[0][0]+ trans_coords[0][2]*trans_coords[0][2]));

    #print theta_y *180/3.14
    R_y = [[cos(theta_y), 0, sin(theta_y)],[0, 1,0], [-1*sin(theta_y), 0, cos(theta_y)]]
    R_y_np = np.array(R_y)
    y_rot =  np.dot(R_y_np, np.matrix.transpose(trans_coords_np))
    #print y_rot
    y_rot_trans = np.matrix.transpose(y_rot)

    #rot in z direction so y = 0 for N
    if y_rot_trans[0][1] <= 0 and y_rot_trans[0][0] > 0:
        #theta_z = 270.0*pi/180.0 + acos(abs(y_rot_trans[0][1])/sqrt(y_rot_trans[0][0]*y_rot_trans[0][0]+ y_rot_trans[0][1]*y_rot_trans[0][1]));
        theta_z = 90.0*pi/180.0 + -1* acos(abs(y_rot_trans[0][1])/sqrt(y_rot_trans[0][0]*y_rot_trans[0][0]+ y_rot_trans[0][1]*y_rot_trans[0][1]));
    elif y_rot_trans[0][1] > 0 and  y_rot_trans[0][0] > 0:
        theta_z = -1*asin(abs(y_rot_trans[0][1])/sqrt(y_rot_trans[0][0]*y_rot_trans[0][0]+ y_rot_trans[0][1]*y_rot_trans[0][1]));
        #print 180/3.14*theta_z
    elif  y_rot_trans[0][1] <= 0 and  y_rot_trans[0][0] <= 0:
        theta_z = -1*(90.0*pi/180.0 + acos(abs(y_rot_trans[0][1])/sqrt(y_rot_trans[0][0]*y_rot_trans[0][0]+y_rot_trans[0][1]*y_rot_trans[0][1])));
    else: #y pos x neg
        theta_z = 270.0*pi/180.0 + acos(abs(y_rot_trans[0][1])/sqrt(y_rot_trans[0][0]*y_rot_trans[0][0]+ y_rot_trans[0][1]*y_rot_trans[0][1]));

    R_z = [[cos(theta_z), -1* sin(theta_z), 0], [sin(theta_z), cos(theta_z), 0], [0, 0, 1]]
    R_z_np = np.array(R_z)
    z_rot =  np.dot(R_z_np, y_rot)
    #print z_rot
    z_rot_trans = np.matrix.transpose(z_rot)


    #rot in x direction so that C z = 0
    if z_rot_trans[1][2] <= 0 and z_rot_trans[1][1] > 0:
        theta_x = 90.0*pi/180.0 + -1*acos(abs(z_rot_trans[1][2])/sqrt(z_rot_trans[1][1]*z_rot_trans[1][1]+ z_rot_trans[1][2]*z_rot_trans[1][2]));
    elif z_rot_trans[1][2] > 0 and  z_rot_trans[1][1] > 0:
        theta_x = -1*asin(abs(z_rot_trans[1][2])/sqrt(z_rot_trans[1][1]*z_rot_trans[1][1]+ z_rot_trans[1][2]*z_rot_trans[1][2]));
    elif  z_rot_trans[1][2] <= 0 and  z_rot_trans[1][1] <= 0:
        #print 180/3.14* acos(abs(z_rot_trans[1][2])/sqrt(z_rot_trans[1][1]*z_rot_trans[1][1]+z_rot_trans[1][2]*z_rot_trans[1][2]))
        theta_x = (90.0*pi/180.0 + acos(abs(z_rot_trans[1][2])/sqrt(z_rot_trans[1][1]*z_rot_trans[1][1]+z_rot_trans[1][2]*z_rot_trans[1][2])));
        #theta_x = -1*(acos(abs(z_rot_trans[1][2])/sqrt(z_rot_trans[1][1]*z_rot_trans[1][1]+z_rot_trans[1][2]*z_rot_trans[1][2])));
    else: #y pos x neg
        theta_x = 270.0*pi/180.0 - acos(abs(z_rot_trans[1][2])/sqrt(z_rot_trans[1][1]*z_rot_trans[1][1]+ z_rot_trans[1][2]*z_rot_trans[1][2]));

    #print theta_x*180/3.14
    R_x = [[1, 0, 0], [0, cos(theta_x), -1* sin(theta_x)], [0, sin(theta_x), cos(theta_x)]]
    R_x_np = np.array(R_x)
    x_rot =  np.dot(R_x_np, z_rot)
    #print x_rot
    flip = 1
    if x_rot[0][0] >0:
        flip = -1;

    S1 = np.array([flip*S1[0], flip*S1[1], S1[2]])
    S2 = np.array([flip*S2[0], flip*S2[1], S2[2]])
    S3 = np.array([flip*S3[0], flip*S3[1], S3[2]])
    R_x_inv = np.linalg.inv(R_x_np)
    R_z_inv = np.linalg.inv(R_z_np)
    R_y_inv = np.linalg.inv(R_y_np)
    S1 = np.dot(R_y_inv, np.dot(R_z_inv, np.dot(R_x_inv, S1)))
    S2 = np.dot(R_y_inv, np.dot(R_z_inv, np.dot(R_x_inv, S2)))
    S3 = np.dot(R_y_inv, np.dot(R_z_inv, np.dot(R_x_inv, S3)))
    S1_trans= [0.0, 0.0, 0.0]
    S2_trans= [0.0, 0.0, 0.0]
    S3_trans= [0.0, 0.0, 0.0]
    for j in range(0, 3):
        S1_trans[j] = S1[j] - trans_vec[j]
        S2_trans[j] = S2[j] - trans_vec[j]
        S3_trans[j] = S3[j] - trans_vec[j]
    #print S1_trans

    s_coords = [S1_trans, S2_trans, S3_trans]
    return s_coords

def findRotMatrix(coords):
    ref_mat = []
    ref_mat.append([-1.458, 0.0, 0.0])
    ref_mat.append([.551, -1.198, -.766])
    ref_mat.append([0.0, 0.0, 0.0])
    ref_mat_np = np.array(ref_mat)
    ref_centroid = np.sum(ref_mat_np, 0);
    ref_centroid = ref_centroid/3.0
    trans_ref_mat = []
    trans_ref_mat.append([0.0, 0.0, 0.0])
    trans_ref_mat.append([0.0, 0.0, 0.0])
    trans_ref_mat.append([0.0, 0.0, 0.0])
    for i in range(0, len(ref_mat)):
        for j in range(0, len(ref_mat[i])):
            trans_ref_mat[i][j] = ref_mat[i][j] - ref_centroid[j]
    trans_ref_mat_np = np.array(trans_ref_mat);
    #print trans_ref_mat_np
    #print ""


    for i in range(0, len(coords)):
        for j in range(0, len(coords[i])):
            coords[i][j] = float(coords[i][j])
    coords_np = np.array(coords);
    coords_centroid = np.sum(coords_np, 0)
    coords_centroid = coords_centroid/3.0
    trans_vec =  coords_centroid

    trans_coords = []
    trans_coords.append([0.0, 0.0, 0.0])
    trans_coords.append([0.0, 0.0, 0.0])
    trans_coords.append([0.0, 0.0, 0.0])
    for i in range(0, len(coords)):
        for j in range(0, len(coords[i])):
            trans_coords[i][j] = float(coords[i][j]) - float(trans_vec[j])
    trans_coords_np = np.array(trans_coords);
    #print trans_coords_np
    #print ""

    A_matrix = np.dot(np.matrix.transpose(trans_coords_np), trans_ref_mat_np)
    #A_matrix = np.dot(np.matrix.transpose(trans_ref_mat_np), trans_coords_np)

    (V, S, W_T) = np.linalg.svd(A_matrix)
    det = findDet(np.dot(np.matrix.transpose(W_T), np.matrix.transpose(V)))

    if det > 0:
        d = 1;
    else:
        d = -1;
    d_matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, d]])
    U = np.dot(np.dot(np.matrix.transpose(W_T), d_matrix), np.matrix.transpose(V))
    #print np.dot(U,np.matrix.transpose(trans_coords_np))
    #print (np.dot(U,trans_coords_np) - trans_ref_mat_np) * (np.dot(U,trans_coords_np) - trans_ref_mat_np)
    #print sqrt(sum(sum((np.dot(U,trans_coords_np) - trans_ref_mat_np) * (np.dot(U,trans_coords_np) - trans_ref_mat_np))))
    exit()
    return (U, trans_vec)


def find_dist(c1, c2):
    return sqrt((c1[0]-c2[0])*(c1[0]-c2[0]) + (c1[1]-c2[1])*(c1[1]-c2[1]) + (c1[2]-c2[2])*(c1[2]-c2[2]))

def calc_score(ideal_dist, actual_dist):
    score = 0;
    for i in range(0, len(ideal_dist)):
        for j in range(0, len(ideal_dist[0])):
            score += abs(ideal_dist[i][j] - actual_dist[i][j])
            if score > score_threshold:
                break
    return score;


def make_dist_mat(c1s, c2s, c3s, c4s, i1, i2, i3, i4):
    best_rotamer_index = [];
    best_score = 1000000000;
    best_actual_mat = [];
    #calc distances first? or in a different method all together and reuse??

    i = 0;
    j = 0;
    k =0;
    l = 0;
    dist1_2 = distance_dict[i1][i2]
    dist1_3 = distance_dict[i1][i3]
    dist1_4 = distance_dict[i1][i4]
    score = 0
    for c1 in c1s:
        rot_dist = [[0, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 0]];
        j=0
        for c2 in c2s:
            rot_dist[0][1] =dist1_2[i][j]
            #rot_dist[1][0] = rot_dist[0][1]
            score_1 = 2*abs(ideal_dists[0][1] - rot_dist[0][1])
            if score_1 > score_threshold:
                j+=1
                continue
            k=0
            for c3 in c3s:
                score_2 = score_1;
                rot_dist[0][2] =dist1_3[i][k]
                rot_dist[2][0] = rot_dist[0][2]
                rot_dist[1][2] = distance_dict[i2][i3][j][k]
                rot_dist[2][1] = rot_dist[1][2]
                score_2 += 2*abs(ideal_dists[0][2] - rot_dist[0][2])
                score_2 += 2*abs(ideal_dists[1][2] - rot_dist[1][2])
                if score_2 > score_threshold:
                    k+=1
                    continue
                l=0;
                for c4 in c4s:
                    score_3 = score_2
                    rot_dist[0][3] = dist1_4[i][l]
                    rot_dist[3][0] = rot_dist[0][3]
                    rot_dist[1][3] = distance_dict[i2][i4][j][l]
                    rot_dist[3][1] = rot_dist[1][3]
                    rot_dist[2][3] = distance_dict[i3][i4][k][l]
                    rot_dist[3][2] = rot_dist[3][2]
                    score_3 += 2*abs(ideal_dists[0][3] - rot_dist[0][3])
                    score_3 += 2*abs(ideal_dists[1][3] - rot_dist[1][3])
                    score_3 += 2*abs(ideal_dists[2][3] - rot_dist[2][3])
                    if score_3 > score_threshold:
                        l+=1
                        continue
                    #rot_dist = []
                   # rot_dist.append([0, distance_dict[i1][i2][i][j], distance_dict[i1][i3][i][k], distance_dict[i1][i4][i][l]])
                    #rot_dist.append([distance_dict[i1][i2][i][j], 0, distance_dict[i2][i3][j][k], distance_dict[i2][i4][j][l]])
                    #rot_dist.append([distance_dict[i1][i3][i][k], distance_dict[i2][i3][j][k], 0, distance_dict[i3][i4][k][l]])
                    #rot_dist.append([distance_dict[i1][i4][i][l], distance_dict[i2][i4][j][l], distance_dict[i3][i4][k][l], 0])

                    #score = calc_score(rot_dist, ideal_dists)
                    if score_3 < best_score:
                        rotamer_index =[i, j, k, l]
                        best_score = score_3;
                        best_actual_mat = rot_dist;
                        best_rotamer_index = rotamer_index;
                    l+=1;
                k+=1;
            j+=1;
        i+=1;
    return (best_score, best_actual_mat, best_rotamer_index)


