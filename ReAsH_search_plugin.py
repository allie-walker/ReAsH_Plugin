#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Allison
#
# Created:     08/10/2015
# Copyright:   (c) Allison 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import Tkinter
import tkFileDialog
import pymol
import time
import numpy as np
import ReAsH_search_tools as RST
import operator
from math import *
results_list = []
results_dict = {}
aa_alphabet = ['GLY', 'ALA', 'VAL', 'LEU','ILE', 'MET', 'CYS', 'PRO', 'TYR', 'PHE', 'TRP', 'HIS', 'ARG', 'LYS', 'GLN', 'GLU', 'ASP', 'ASN', 'SER', 'THR']
distance_dict = {}
distances={}
allowed_pairs = {}
score_threshold = 30



def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command',
                             'ReAsH Search',
                             label = 'ReAsH Search',
                             command = lambda s=self : ReAsHSearchPlugin(s))

def filterByTypeSite(results_list, site_type):
    return results_list

#Convert a reash binding site search result to a pymol selection command
def resultToSelection(result, aa_numbers, chains):
    residues = result[1]
    print residues
    selection_string = "("
    for i in range(0, 4):
        if i != 0:
            selection_string += "+ "

        selection_string += "chain " + str(chains[int(residues[i])]) + " and resi " + str(aa_numbers[int(residues[i])])

    selection_string +=")"
    return selection_string

#filter amino acids in the pymol object to all have solvent accessible surface area (SASA) above the specified minimum
def filterSasa(min_sasa, target_obj, aa_nums,aa_names,chains):
    print "filteirng SASA"
    sasa_too_low = []
    filtered_aa = []
    filtered_chains = []
    pymol.cmd.set('dot_solvent', 1)
    pymol.cmd.set('dot_density', 3)
    for i in range (0, len(aa_nums)):
        area = 0
        area = float(pymol.cmd.get_area(target_obj + " and chain " + chains[i] + " and resi " + aa_nums[i] ))
        if area < min_sasa:
            sasa_too_low.append(i)
        else:
            filtered_aa.append(aa_nums[i])
            filtered_chains.append(chains[i])
    print sasa_too_low
    return (sasa_too_low, filtered_aa, filtered_chains)

#filter amino acids in the pymol object to all have b-factor above the specified minimum
def filterBFactor(min_b, target_obj, aa_nums, aa_names,chains):
    b_factor_too_low = []
    filtered_aa = []
    filtered_chains = []
    for i in range(0, len(aa_nums)):
        sele_string = target_obj + " and resi " + aa_nums[i] + " and name CA"
        m = pymol.cmd.get_model(sele_string)
        if m.atom ==[]:
            b_factor = 0
            continue
        b_factor = float(m.atom[0].b)
        if b_factor < min_b:
            b_factor_too_low.append(i)
        else:
            filtered_aa.append(aa_nums[i])
            filtered_chains.append(chains[i])
    return (b_factor_too_low, filtered_aa, filtered_chains)

#calculate the possible sulfur coordinates when a residue is mutated to a cysteine
def calcSulfurCoords(coords):
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

def sulfurCoords(bb_x,bb_y,bb_z,bb_name,b_too_low,sasa_too_low):
    sulfur_coords = [];

    i =0
    print "finding possible sulfur coordinates"
    ineligible_residues = 0
    while i < len(bb_name) -2:
        if i/3 in b_too_low or i/3 in sasa_too_low:
            i+=3
            ineligible_residues += 1
            continue
        print "num ineligible residues: " + str(ineligible_residues)

        N_x = float(bb_x[i])
        N_y = float(bb_y[i])
        N_z = float(bb_z[i])

        CA_x = float(bb_x[i+1])
        CA_y = float(bb_y[i+1])
        CA_z = float(bb_z[i+1])

        C_x = float(bb_x[i+2])
        C_y = float(bb_y[i+2])
        C_z = float(bb_z[i+2])

        coord_mat = [];
        coord_mat.append([N_x,N_y, N_z])
        coord_mat.append([C_x,C_y, C_z])
        coord_mat.append([CA_x,CA_y, CA_z])
        #print coord_mat

        s_coords = calcSulfurCoords(coord_mat);
        #s_coords = []

        i+=3

        sulfur_coords.append(s_coords)
    return sulfur_coords

def makeDistDict(sulfur_coords):
    i = 0
    for c1 in sulfur_coords:
        j = 0
        allowed_list = []
        if i not in distances:
            distances[i] = {}
        if i not in distance_dict:
            distance_dict[i] = {}
        for  j in range(i+1, len(sulfur_coords)):
            c2 = sulfur_coords[j]
            if RST.dist(c1[0],c2[0]) <= 15:
                allowed_list.append(j)
                distance_dict[i][j] = {}
                if j not in distance_dict:
                    distance_dict[j] = {}
                distance_dict[j][i] = {}
                distance_dict[j][i][0] = []
                distance_dict[i][j][0] = []
                distance_dict[j][i][1] = []
                distance_dict[i][j][1] = []
                distance_dict[j][i][2] = []
                distance_dict[i][j][2] = []
                for k in range(0, 3):
                    for l in range(0, 3):
                        distance_dict[i][j][k].append(RST.dist(c1[k],c2[l]))
                        distance_dict[j][i][k].append(RST.dist(c1[l],c2[k]))
            distances[i][j] =RST.dist(c1[0],c2[0])
            if j not in distances:
                distances[j] = {}
            distances[j][i] = RST.dist(c1[0],c2[0])
            j+=1;
        allowed_pairs[i] = allowed_list
        i +=1
    #return (allowed_pairs, distances)

def divideList(coord_list, n):
    for i in xrange(0, len(coord_list), n):
        yield (coord_list[i:i+n], i)

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
    dist2_3 = distance_dict[i2][i3]
    dist2_4 = distance_dict[i2][i4]
    dist3_4 = distance_dict[i3][i4]
    score = 0
    for c1 in c1s:
        j=0
        dist1_2_i = dist1_2[i]
        dist1_3_i = dist1_3[i]
        dist1_4_i = dist1_4[i]
        for c2 in c2s:
            score_1 = 2*abs(3.6 - dist1_2_i[j])
            if score_1 > score_threshold or score_1 > best_score:
                j+=1
                continue
            k=0
            dist2_3_j = dist2_3[j]
            dist2_4_j = dist2_4[j]
            for c3 in c3s:
                score_2 = score_1 + 2*abs(5.1 - dist1_3_i[k]) + 2*abs(7.2 - dist2_3_j[k]);
                if score_2 > score_threshold or score_2 > best_score:
                    k+=1
                    continue
                l=0;
                for c4 in c4s:
                    score_3 = score_2 + 2*abs(7.2 - dist1_4_i[l]) + 2*abs(7.6 -  dist2_4_j[l]) + 2*abs( 3.6 -dist3_4[k][l])
                    if score_3 > score_threshold:
                        l+=1
                        continue
                    if score_3 < best_score:
                        rotamer_index =[i, j, k, l]
                        best_score = score_3;
                        best_rotamer_index = rotamer_index;
                    l+=1;
                k+=1;
            j+=1;
        i+=1;
    return (best_score, best_rotamer_index)

def siteSearch(all_coords, allowed_pairs):
        i = 0;
        j = 0;
        k = 0;
        l = 0;

        results = []


        for c1 in all_coords:
            print i
            for j in allowed_pairs[i]:
                c2 = all_coords[j]
                for k in allowed_pairs[i]:
                    if k == j:
                        continue
                    if distances[j][k] >= 15:
                        continue;
                    c3 = all_coords[k]
                    for l in allowed_pairs[i]:
                        if l ==k or l ==j:
                            continue
                        if distances[j][l] >15 or distances[k][l] >15:
                            continue;
                        c4 = all_coords[l]
                        (dist_matrix_score, rotamer_indices) = make_dist_mat(c1, c2, c3, c4, i, j, k, l)

                        if dist_matrix_score > 30:
                            continue
                        index = [i, j, k, l]
                        results.append((dist_matrix_score, index, rotamer_indices))
            i+=1;
        return results

def isNumber(input_str):
    try:
        float(input_str)
        return True
    except ValueError:
        return False

def newSearch(target_obj, search_name, min_sasa, min_b,site_type):
    bb_x = [];
    bb_y = [];
    bb_z = [];
    bb_bscore = [];
    bb_name = [];
    aa_numbers = []
    aa_names = []
    chains = []

    #get list of all bb atoms in the selection
    pymol.stored.list = []
    pymol.cmd.iterate("("+ target_obj+" and (name ca + name c + name n))","stored.list.append((chain,resi,resn,name))")
    atom_list = pymol.stored.list
    print atom_list
    pymol.stored.list = []

    #print atom_list
    previous_atom = ""

    for atom in atom_list:
        if atom[2] not in aa_alphabet:
            print atom[2]
            continue

        if atom[3] == previous_atom:
            #print "HERE"
            continue
        pymol.stored.list = []
        sele_string ="("+ target_obj+" and chain "+ atom[0] +" and resi "+atom[1] + " and name " + atom[3] + ")"
        pymol.cmd.iterate_state(1, sele_string, 'stored.list.append((x,y,z))')
        bb_x.append(pymol.stored.list[0][0])
        bb_y.append(pymol.stored.list[0][1])
        bb_z.append(pymol.stored.list[0][2])

        bb_name.append(atom[3])

        #print atom[3]
        #print atom[1]
        if atom[3] == 'N':
            chains.append(atom[0])
            aa_numbers.append(atom[1])
            aa_names.append(atom[2])
        previous_atom = atom[3]

    (sasa_too_low, filtered_aa1, filtered_chains1) = filterSasa(min_sasa,target_obj,aa_numbers,aa_names,chains)
    (b_too_low, filtered_aa2, filtered_chains2) = filterBFactor(min_sasa,target_obj,aa_numbers,aa_names,chains)
    filtered_aa = []
    filtered_chains = []
    i = 0
    j = 0


    while i < len(filtered_aa1) and j < len(filtered_aa2):
        #print filtered_aa1[i]
        #print filtered_aa1[j]

        if int(filtered_aa1[i]) == int(filtered_aa2[j]):
            filtered_aa.append(filtered_aa1[i])
            filtered_chains.append(filtered_chains1[i])
            j+=1
            i+=1
            continue
        elif int(filtered_aa1[i]) < int(filtered_aa2[j]):
            i +=1
            continue
        else:
            j+=1
            continue

    sulfur_coords = sulfurCoords(bb_x,bb_y,bb_z,bb_name,b_too_low,sasa_too_low)
    #print sulfur_coords

    i = 0;
    j = 0;
    k = 0;
    l = 0;

    print "building distance dictionary"
    makeDistDict(sulfur_coords)

    i = 0;
    j = 0;
    k = 0;
    l = 0;
    #partial_siteSearch = partial(siteSearch, all_coords=sulfur_coords, allowed_pairs=allowed_pairs)

    print "starting search"
    results = siteSearch(sulfur_coords,allowed_pairs)
    sorted_results = sorted(results, key=lambda score:score[0])

    print "done"
    #still need to filter out for correct site type (intramolecular, etc)!!
    sorted_results = filterByTypeSite(sorted_results, site_type)
    results_obj = SearchResultList(search_name, target_obj, sorted_results,filtered_aa, filtered_chains)
    return results_obj

def getAverageSasa(residue_list, sasa_dict):
    avg_sasa = 0.0
    #return avg_sasa
    for res in residue_list:
        #print res
        if res in sasa_dict:
            avg_sasa += sasa_dict[res]
        else:
            print "here"
            avg_sasa += 0
    return avg_sasa/4.0

def getAverageBFactor(residue_list, b_dict):
    avg_b =0.0
    for res in residue_list:
        #sele_string = selected_obj + " and resi " + " and chain " + res.split(",")[0] + " and resi " + res.split(",")[1] + " and name CA"
        #m = pymol.cmd.get_model(sele_string)
        if res in b_dict:
            avg_b += b_dict[res]
        else:
            avg_b += 0
    return avg_b/4

#return dictionary of all b factors
#keys are chain letter, str(aa_index)
def getAllBFactors(aa_list, chain_list, selected_obj):
    b_dict = {}
    i = 0
    for a in aa_list:
        sele_string = selected_obj + " and chain " + chain_list[i] + " and resi " + str(a) + " and name CA"
        #print sele_string
        #return b_dict
        m = pymol.cmd.get_model(sele_string)
        if m.atom ==[]:
            b_dict[chain_list[i] + "," + str(a)] = 0
            continue
        b_dict[chain_list[i] + "," + str(a)] = float(m.atom[0].b)
        i+=1
    return b_dict

def getAllSASA(aa_list, chain_list, selected_obj):
    sasa_dict = {}
    i = 0
    for a in aa_list:
        #what to do about Gly?
        sele_string = selected_obj + " and chain " + chain_list[i] + " and resi " + str(a) + " and (name CA + name CB)"
        sasa_dict[chain_list[i] + "," + str(a)] = float(pymol.cmd.get_area(sele_string))
        i+=1
    return sasa_dict


class SearchResultList:
    def __init__(self, name, pymol_ojbect, results, aa_list, chain_list):
        self.name = name
        self.pymol_ojbect = pymol_ojbect #name of the object/selection in pymol search was done on
        self.results = results
        self.aa_list= aa_list
        self.chain_list = chain_list
    def getPymolObject(self):
        return self.pymol_ojbect
    def getName(self):
        return self.name
    def getAAList(self):
        return self.aa_list
    def getResultsList(self):
        return self.results
    def getChainsList(self):
        return self.chain_list

class SearchErrorDialog:
    def __init__(self, parent, message):
        top = self.top = Tkinter.Toplevel(parent)
        Tkinter.Label(top,text=message).pack(pady=2)
        b=Tkinter.Button(top, text= "OK", command=self.ok)
        b.pack(pady=2)
    def ok(self):
        self.top.destroy()

class WaitingDialog:
    def __init__(self,parent):
        top=self.top = Tkinter.Toplevel(parent)
        message = "Waiting for calculations to finish"
        Tkinter.Label(top,text=message).pack(pady=2)


class ReAsHSearchPlugin:
         #   min_sasa_float = float(min_sasa_str)
        #print results_list


    def getPymolObjects(self):
        #selection_list =pymol.cmd.get_names('selections')
        object_list =  pymol.cmd.get_names('objects')
        #sele_and_obj_list = []
        #sele_and_obj_list.extend(object_list)
        #sele_and_obj_list.extend(selection_list)
        return object_list


    def nameDropDownEvents(self,event):
        #need to update list everytime it is clicked
        sele_and_obj_list = self.getPymolObjects()
        print sele_and_obj_list
        self.obj_name_options['menu'].delete(0, Tkinter.END)
        if len(sele_and_obj_list)> 0:
            self.selected_object.set(sele_and_obj_list[0])
            for obj in sele_and_obj_list:
                print self.obj_name_options['menu'].add_command(label=obj,command=lambda v=self.selected_object,l=obj:v.set(l))
        else:
            self.selected_object.set("")



    def startNewSearch(self):
        #check that an object has been selected
        if self.selected_object == "":
            d = SearchErrorDialog(self.window,"Please select an object to run the search on")

        #check that object is valid?
        elif self.selected_object.get() not in self.getPymolObjects():
            d = SearchErrorDialog(self.window, "Please choosing and existing object or selection")
        #check that a name has been entered into the name text
        elif len(self.search_name_text.get()) == 0 or len(self.search_name_text.get().replace(" ", "")) == 0:
            d = SearchErrorDialog(self.window, "Please enter a name for your search")
        #check that the name is not the same as any existing searches
        elif self.search_name_text in results_list:
            d = SearchErrorDialog(self.window, "This name has already been used")
        elif not isNumber(self.SASA_text_box.get()) or not isNumber(self.B_factor_label_text_box.get()):
            d = SearchErrorDialog(self.window, "Please enter a number for min SASA and B-factor")

        #open the serach results screen
        else:
            #Open a progress dialog box, wait for the search to be done, then close the progress box and open the results box
            #d = WaitingDialog(self.window)
            new_results = newSearch(self.selected_object.get(), self.search_name_text.get(), float(self.SASA_text_box.get()), float(self.B_factor_label_text_box.get()),self.site_type_var.get())
            #print results_list
            #print new_results
            results_list.append(new_results.getName())
            results_dict[new_results.getName()] = new_results
            b_dict = getAllBFactors(new_results.getAAList(), new_results.getChainsList(), new_results.getPymolObject())
            sasa_dict = getAllSASA(new_results.getAAList(), new_results.getChainsList(), new_results.getPymolObject())
            results = ResultsDisplay(self.app, results_list, results_dict, b_dict, sasa_dict)

    def viewCurrentResults(self):
        results = ResultsDisplay(self.app, results_list, results_dict, {}, {})

    def __init__(self,app):
        self.app = app
        self.parent = app.root
        window =Tkinter.Tk()
        self.window= window
        window.title("ReAsH Search")

        welcome_banner_frame = Tkinter.Frame(window)
        welcome_banner_frame.pack(fill=Tkinter.BOTH, expand=1)
        welcome_banner_frame.config(background="red")
        welcome_banner_frame.config(width=70)
        welcome_banner = Tkinter.Label( welcome_banner_frame, text="ReAsH Search Program \n Developed by Allison Walker \n Schepartz Lab \n Yale University", font=("Helvetica", 16),justify=Tkinter.LEFT)
        welcome_banner.configure(background="red")
        #welcome_banner.config(width=50)
        welcome_banner.pack(side=Tkinter.LEFT)

        self.name_frame = Tkinter.Frame(window)
        self.name_frame.pack(fill=Tkinter.BOTH, expand=1,pady=2)
        obj_name_lbl = Tkinter.Label(self.name_frame, text="Run search on: ").pack(side=Tkinter.LEFT)

        self.selected_object = Tkinter.StringVar(window)
        self.obj_name_options = Tkinter.OptionMenu(self.name_frame, self.selected_object, "")  #change this to get the objects in the pymol file!!
        sele_and_obj_list = self.getPymolObjects()
        if len(sele_and_obj_list)>0:
            self.selected_object.set(sele_and_obj_list[0])
            self.obj_name_options = Tkinter.OptionMenu(self.name_frame, self.selected_object, *sele_and_obj_list)
            #self.obj_name_options = apply(OptionMenu, (self.name_frame, self.selected_object) + tuple(obj_name_options))
        self.obj_name_options.config(width=20)
        self.obj_name_options.bind("<Button-1>", self.nameDropDownEvents)
        self.obj_name_options.pack(side=Tkinter.LEFT)

        search_name_frame = Tkinter.Frame(window)
        search_name_frame.pack(fill=Tkinter.BOTH,expand=1)
        search_name_label = Tkinter.Label(search_name_frame,text="Name of Search: ")
        search_name_label.pack(side=Tkinter.LEFT)
        self.search_name_text = Tkinter.Entry(search_name_frame)
        self.search_name_text.config(width = 20)
        self.search_name_text.pack(side=Tkinter.LEFT)

        options_frame =Tkinter.Frame(window)
        options_frame.pack(fill=Tkinter.BOTH, expand =1)
        options_label = Tkinter.Label(options_frame, text= "Options:")
        options_label.pack(side=Tkinter.LEFT)

        option_values_frame = Tkinter.Frame(window)
        option_values_frame.pack(fill=Tkinter.BOTH, expand =1,padx=5)
        SASA_label = Tkinter.Label(option_values_frame, text ="    Minimum SASA:  ")
        SASA_label.pack(side=Tkinter.LEFT)

        self.SASA_text_box =Tkinter.Entry(option_values_frame)
        self.SASA_text_box.insert(0, "0")
        self.SASA_text_box.config(width =15)
        self.SASA_text_box.pack(side=Tkinter.LEFT)

        B_factor_label = Tkinter.Label(option_values_frame, text = "   Minimum B factor:   ")
        B_factor_label.pack(side = Tkinter.LEFT)

        self.B_factor_label_text_box = Tkinter.Entry(option_values_frame)
        self.B_factor_label_text_box.insert(0, "0")
        self.B_factor_label_text_box.config(width = 15)
        self.B_factor_label_text_box.pack(side=Tkinter.LEFT)

        site_type_frame = Tkinter.Frame(window)
        site_type_frame.pack(fill=Tkinter.BOTH, expand =1)
        site_type_label = Tkinter.Label(site_type_frame, text="Site Type:")
        site_type_label.pack(side=Tkinter.LEFT)

        self.site_type_var = Tkinter.IntVar()
        site_type_radio_frame_1 = Tkinter.Frame(window)
        site_type_radio_frame_1.pack(fill=Tkinter.BOTH, expand=1)
        rad1 = Tkinter.Radiobutton(site_type_radio_frame_1, text="Intramolecular Site", variable=self.site_type_var, value=1)
        rad1.pack(side=Tkinter.LEFT)
        site_type_radio_frame_2 = Tkinter.Frame(window)
        site_type_radio_frame_2.pack(fill=Tkinter.BOTH, expand=1)
        rad2 = Tkinter.Radiobutton(site_type_radio_frame_2, text="Intermolecular Site", variable=self.site_type_var, value=2)
        rad2.pack(side=Tkinter.LEFT)
        site_type_radio_frame_3 = Tkinter.Frame(window)
        site_type_radio_frame_3.pack(fill=Tkinter.BOTH, expand=1)
        rad3 = Tkinter.Radiobutton(site_type_radio_frame_3, text="Either Type of Site", variable=self.site_type_var, value=3)
        rad3.pack(side=Tkinter.LEFT)
        rad3.select()


        self.start_new_search = Tkinter.Button(window, text="Start New Search", command=self.startNewSearch)
        self.start_new_search.pack(pady=2)

        self.view_current_results = Tkinter.Button(window, text="View Current Results", command=self.viewCurrentResults)
        self.view_current_results.pack(pady=2)



class ResultsDisplay():
    def resultsDropDownEvents(self,event):
        #need to update list everytime it is clicked
        print "HERE!!"

    def loadFromFile(self):
         in_file= tkFileDialog.askopenfile(mode='r', **self.open_opt)
         if in_file == None:
            return
         #add the new results and list
         pymol_obj = ""
         search_name = ""
         chain_list = []
         aa_and_chain_list = []
         aa_list =[]
         results_list = []
         index = 0
         result_lines = []

         for line in in_file:
            # print line
             line = line.replace("\n", "")
             if len(line.split(",")) == 2:
                pymol_obj = line.split(",")[0]
                search_name = line.split(",")[1]
                if search_name in results_dict:
                    print "HERE"
                    d = SearchErrorDialog(self.parent,"There is already a search by this name")
                    return
             else:
                result_lines.append(line)
                for i in (0, 2, 4, 6):
                    aa_id = (line.split(",")[i+1], int(line.split(",")[i+2]))
                    if aa_id not in aa_and_chain_list:
                        aa_and_chain_list.append(aa_id)


         aa_and_chain_list = sorted(aa_and_chain_list, key = lambda x : (x[0], x[1]))
         #now make an aa_list and a chains_list
         for a in aa_and_chain_list:
               #print a
                chain_list.append(a[0])
                aa_list.append(a[1])

         #now make the results list by using the indices above
         results_list = []
         for l in result_lines:
                score = float(l.split(",")[0])
                aa_index_list = []
                for i in (0, 2, 4, 6):
                    aa_id = (l.split(",")[i+1], int(l.split(",")[i+2]))
                    #print aa_id
                    #print aa_and_chain_list[0]
                    #print aa_and_chain_list.index(aa_id)
                    aa_index_list.append(aa_and_chain_list.index(aa_id))
                rotamers = []
                for i in range(9, len(l.split(","))):
                    rotamers.append(int(l.split(",")[i]))
                #print rotamers
                results_list.append([score, aa_index_list, rotamers])  #make the proper things into ints or numbers!!!!!!!!!!!

         #make the object
         search_results = SearchResultList(search_name, pymol_obj, results_list, aa_list, chain_list)

         #add everything it to the dictionary
         results_dict[search_name] = search_results
         self.search_results_names.append(search_name)
         #return
         #refresh the drop down menu?
         #self.display_results_dropdown = Tkinter.OptionMenu(result_choice_frame, self.selected_results, "") #add a command to this so when a different one is selected results show
         menu = self.display_results_dropdown['menu']
         menu.delete(0, 'end')
         for name in self.search_results_names:
            menu.add_command(label=name, command=lambda search=name: self.selected_results.set(name))
         menu2 = self.list1_dropdown['menu']
         menu2.delete(0, 'end')
         for name in self.search_results_names:
             menu2.add_command(label=name, command=lambda search=name: self.list1_var.set(name))

         menu3 = self.list2_dropdown['menu']
         menu3.delete(0, 'end')
         for name in self.search_results_names:
            menu3.add_command(label=name, command=lambda search=name:self.list2_var.set(name))

         if len(self.search_results_names) == 1:
            self.selected_results.set(self.search_results_names[0])
            #get b dict and sasa dict!
            b_dict = getAllBFactors(search_results.getAAList(), search_results.getChainsList(), search_results.getPymolObject())
            sasa_dict = getAllSASA(search_results.getAAList(), search_results.getChainsList(), search_results.getPymolObject())
            self.fillListBox(results_dict, b_dict, sasa_dict)


    def saveToFile(self):
         #HANDLE CASE WHERE THERE IS NO SEARCH OBJECT SELECTED!!!!!!!!!!!  DOES THIS WORK OK?
         if len(results_dict) == 0:
            return
         outfile = tkFileDialog.asksaveasfile(mode='w', **self.save_opt)
         if outfile == None:
            return

         #print filename
         selected_search = self.selected_results.get()

         #file format will be csv formatted file with
         #first line is object name, search name
         object_name = results_dict[selected_search].getPymolObject()
         search_name = results_dict[selected_search].getName()
         outfile.write(object_name + "," + search_name + "\n")

         #all other lines:
         #score, chain 1, resi 1, chain 2, resi 2, chain 3, resi 3, chain 4, resi 4, rotamer 1, rotamer 2, rotamer 3, rotamer 4
         aa_list = results_dict[selected_search].getAAList()
         chain_list = results_dict[selected_search].getChainsList()


         for r in  self.working_results_list:
            #format of working_results_list: score, aa_index, rotamer index
            print r
            indices = r[1]
            print r[0]
            outfile.write(str(r[0]))#write the score

            #write the chain and resi nums
            for i in indices:
             outfile.write("," + str(chain_list[i]) + "," + str(aa_list[i]))

            #now write the rotamers
            rotamers = r[2]
            for rot in rotamers:
                outfile.write("," + str(rot))
            outfile.write("\n")
         outfile.close()

         #how should I deal with people saving over previously saved list that might be open currently in editor?

    def restoreListToDefualt(self):
        print "blah"


    def filterAromatic(self,event):
        if self.require_aromatic.get() == 0:
            print "require aromatics"
        else:
            print "don't requrie aromatics"

    def filterTrp(self, event):
        if self.disallow_trp.get() == 0:
            print "don't allow trp"
        else:
            print "allow trp"

    def filterPos(self, event):
        if self.require_pos.get() == 0:
            print "require nearby positive residue"
        else:
            print "don't require positive residue"

    def filterNeg(self, event):
        if self.disallow_neg.get() == 0:
            print "disallow nearby negative residue"
        else:
            print "don't disallow negative residue"

    def filterRotamer(self, event):
        if self.require_allowed_rots.get() == 0:
            print "require allowed rotamers"
        else:
            print "don't require allowed rotamers"

    def changeSort(self):
        print self.sort_by_var.get()
        if self.sort_by_var.get() == 1:
            print "sort by score"
            self.fillListBox(results_dict, self.b_dict, self.sasa_dict)
        elif self.sort_by_var.get() == 2:
            print "sort by sasa"
            self.fillListBox(results_dict, self.b_dict, self.sasa_dict)
        else:
            print "sort by b factor"
            self.fillListBox(results_dict, self.b_dict, self.sasa_dict)

    def displayAromatic(self):
        if self.aromatic_display.get() == 1:
            print "display aromatic residues"
            pymol.cmd.delete("near_phe")
            line_text = self.results_listbox.get(Tkinter.ANCHOR)
            if len(line_text) == 0:
                return
            if "Residues" in line_text:
                return

            selected_search = self.selected_results.get()
            selection_text = ""
            for i in range(0, 4):
                selection_text += results_dict[selected_search].getPymolObject()  + " and "
                selection_text += "chain " + line_text.split(",")[i].split()[0]
                selection_text += " and resi " + line_text.split(",")[i].split()[1]
                if i != 3:
                    selection_text += " + "

            pymol.cmd.select("site", selection_text)
            pymol.cmd.select("near", "site around 5")
            pymol.cmd.select("near_phe_atom", "near and (resn Phe)")
            pymol.cmd.select("near_phe","br. near_phe_atom")
            pymol.cmd.show(representation="sticks", selection="near_phe")
            #delete site and near slections
            pymol.cmd.delete("site")
            pymol.cmd.delete("near")
            pymol.cmd.delete("near_phe_atom")
            pymol.cmd.show(representation="sticks", selection="near_phe")
        else:
            if "near_phe" in pymol.cmd.get_names():
                pymol.cmd.hide(representation="sticks", selection="near_phe")
            pymol.cmd.delete("near_phe")
            print "don't dipslay aromatic residues"

    def displayTrp(self):
        if self.trp_display.get() == 1:
            pymol.cmd.delete("near_trp_tyr")
            print "dipslay Trp"
            line_text = self.results_listbox.get(Tkinter.ANCHOR)
            if len(line_text) == 0:
                return
            if "Residues" in line_text:
                return

            selected_search = self.selected_results.get()
            selection_text = ""
            for i in range(0, 4):
                selection_text += results_dict[selected_search].getPymolObject()  + " and "
                selection_text += "chain " + line_text.split(",")[i].split()[0]
                selection_text += " and resi " + line_text.split(",")[i].split()[1]
                if i != 3:
                    selection_text += " + "

            pymol.cmd.select("site", selection_text)
            pymol.cmd.select("near", "site around 5")
            pymol.cmd.select("near_trp_atom", "near and (resn Trp + resn Tyr)")
            pymol.cmd.select("near_trp_tyr","br. near_trp_atom")
            pymol.cmd.show(representation="sticks", selection="near_trp_tyr")
            #delete site and near slections
            pymol.cmd.delete("site")
            pymol.cmd.delete("near")
            pymol.cmd.delete("near_trp_atom")
            pymol.cmd.show(representation="sticks", selection="near_trp_tyr")
        else:
            if "near_trp_tyr" in pymol.cmd.get_names():
                pymol.cmd.hide(representation="sticks", selection="near_trp_tyr")
            print "don't display Trp"
            pymol.cmd.delete("near_trp_tyr")

    def displayPos(self):
        if self.pos_display.get() == 1:
            pymol.cmd.delete("near_pos")
            print "display positive"
            line_text = self.results_listbox.get(Tkinter.ANCHOR)
            if len(line_text) == 0:
                return
            if "Residues" in line_text:
                return

            selected_search = self.selected_results.get()
            selection_text = ""
            for i in range(0, 4):
                selection_text += results_dict[selected_search].getPymolObject()  + " and "
                selection_text += "chain " + line_text.split(",")[i].split()[0]
                selection_text += " and resi " + line_text.split(",")[i].split()[1]
                if i != 3:
                    selection_text += " + "

            pymol.cmd.select("site", selection_text)
            pymol.cmd.select("near", "site around 5")
            pymol.cmd.select("near_pos_atom", "near and (resn Lys + resn Arg)")
            pymol.cmd.select("near_pos","br. near_pos_atom")
            pymol.cmd.show(representation="sticks", selection="near_pos")
            #delete site and near slections
            pymol.cmd.delete("site")
            pymol.cmd.delete("near")
            pymol.cmd.delete("near_pos_atom")
            pymol.cmd.show(representation="sticks", selection="near_pos")
        else:
            if "near_pos" in pymol.cmd.get_names():
                pymol.cmd.hide(representation="sticks", selection="near_pos")
            print "don't display positive"
            pymol.cmd.delete("near_pos")

    def displayNeg(self):
        if self.neg_display.get() == 1:
            print "display negative"
            pymol.cmd.delete("near_neg")
            line_text = self.results_listbox.get(Tkinter.ANCHOR)
            if len(line_text) == 0:
                return
            if "Residues" in line_text:
                return

            selected_search = self.selected_results.get()
            selection_text = ""
            for i in range(0, 4):
                selection_text += results_dict[selected_search].getPymolObject()  + " and "
                selection_text += "chain " + line_text.split(",")[i].split()[0]
                selection_text += " and resi " + line_text.split(",")[i].split()[1]
                if i != 3:
                    selection_text += " + "

            pymol.cmd.select("site", selection_text)
            pymol.cmd.select("near", "site around 5")
            pymol.cmd.select("near_neg_atom", "near and (resn Glu + resn Asp)")
            pymol.cmd.select("near_neg","br. near_neg_atom")
            pymol.cmd.show(representation="sticks", selection="near_neg")
            #delete site and near slections
            pymol.cmd.delete("site")
            pymol.cmd.delete("near")
            pymol.cmd.delete("near_neg_atom")
            pymol.cmd.show(representation="sticks", selection="near_neg")
        else:
            if "near_neg" in pymol.cmd.get_names():
                pymol.cmd.hide(representation="sticks", selection="near_neg")
            print "don't display negative"
            pymol.cmd.delete("near_neg")

    #display all polar residues within 5 angstroms (or 10?)
    def displayPolar(self):
        if self.polar_display.get() == 1:
            pymol.cmd.delete("near_polar")
            print "display polar"
            line_text = self.results_listbox.get(Tkinter.ANCHOR)
            if len(line_text) == 0:
                return
            if "Residues" in line_text:
                return

            selected_search = self.selected_results.get()
            selection_text = ""
            for i in range(0, 4):
                selection_text += results_dict[selected_search].getPymolObject()  + " and "
                selection_text += "chain " + line_text.split(",")[i].split()[0]
                selection_text += " and resi " + line_text.split(",")[i].split()[1]
                if i != 3:
                    selection_text += " + "

            pymol.cmd.select("site", selection_text)
            pymol.cmd.select("near", "site around 5")
            pymol.cmd.select("near_polar_atom", "near and (resn Tyr + resn Ser + resn Thr + resn His + resn Glu + resn Gln + resn Asp + resn Asn + resn Arg + resn Lys + resn Trp)")
            pymol.cmd.select("near_polar","br. near_polar_atom")
            pymol.cmd.show(representation="sticks", selection="near_polar")
            pymol.cmd.delete("site")
            pymol.cmd.delete("near")
            pymol.cmd.delete("near_polar_atom")
            #delete site and near slections
            pymol.cmd.show(representation="sticks", selection="near_polar")
        else:
            if "near_polar" in pymol.cmd.get_names():
                pymol.cmd.hide(representation="sticks", selection="near_polar")
            print "don't display polar"
            pymol.cmd.delete("near_polar")

    def displayReAsH(self, event):
        if self.ReAsH_display.get() == 0:
            print "display ReAsH"
        else:
            print "don't display ReAsH"

    def displayCys(self):
        if self.display_side_chains.get() == 1:
            print "display side chains"
            line_text = self.results_listbox.get(Tkinter.ANCHOR)
            if len(line_text) == 0:
                return
            if "Residues" in line_text:
                return
            selected_search = self.selected_results.get()
            selection_text = ""
            for i in range(0, 4):
                selection_text += results_dict[selected_search].getPymolObject()  + " and "
                if i == 0:
                    selection_text += "chain " + line_text.split(",")[i].split()[1]
                    selection_text += " and resi " + line_text.split(",")[i].split()[2]
                else:
                    selection_text += "chain " + line_text.split(",")[i].split()[0]
                    selection_text += " and resi " + line_text.split(",")[i].split()[1]
                if i != 3:
                    selection_text += " + "
            print selection_text
            pymol.cmd.create("Tetracysteine_model", selection_text, 0, 0)
            pymol.cmd.cmd.show_as(representation="sticks", selection="Tetracysteine_model")
            pymol.cmd.wizard("mutagenesis")
            selected_line= self.results_listbox.curselection()
            selected_line_zero = selected_line[0]
            #print int(selected_line_zero)
            #print self.working_results_list[int(selected_line_zero)]
            rotamers = self.working_results_list[int(selected_line[0])-1][2]
            print rotamers
            print "HERE!!"

            for i in range(0, 4):
                if i == 0:
                    selection_text = "Tetracysteine_model and chain " + line_text.split(",")[i].split()[1] +" and resi " + line_text.split(",")[i].split()[2]
                else:
                    selection_text = "Tetracysteine_model and chain " + line_text.split(",")[i].split()[0] +" and resi " + line_text.split(",")[i].split()[1]
                #selection_text = "Tetracysteine_model and resi " + line_text.split(",")[i].split()[1]
                print selection_text
                rotamer = int(rotamers[i]) +2
                print rotamer
                pymol.cmd.refresh_wizard()
                pymol.cmd.get_wizard().set_dep("ind")
                pymol.cmd.select("sele", selection_text)
                pymol.cmd.get_wizard().do_select('''sele''')

                pymol.cmd.frame(1)
                pymol.cmd.get_wizard().set_mode('CYS')
                print rotamer
                for i in range(2, rotamer):
                    pymol.cmd.get_wizard().do_state(rotamer)
                    pymol.cmd.forward()

                pymol.cmd.get_wizard().do_state(rotamer)
                #pymol.cmd.frame(1)
                pymol.cmd.get_wizard().apply()
                pymol.cmd.select('sele','None')
                #get the rotamer for this result
                #cmd.set_wizard()
            pymol.cmd.set_wizard()
            pymol.cmd.zoom("Tetracysteine_model")
            #now do the mutagensis of each resi in the tetracysteine model
        else:
            #TODO: get rid of the tetracysteine model object
            #check if tetracysteine model exists and delete it if it does
            object_list =  pymol.cmd.get_names('objects')
            if "Tetracysteine_model" in object_list:
                pymol.cmd.delete("Tetracysteine_model")
            print "don't display side chains"

    def displayMeasurements(self):
        if self.display_measurements.get() == 1:
            print "display measurements"
            #delete all previous measurements!!
            #only delete measurements that this program has made
            line_text = self.results_listbox.get(Tkinter.ANCHOR)
            if len(line_text) == 0:
                return
            if "Residues" in line_text:
                return
            pymol.cmd.delete("cys-cys")
            selected_search = self.selected_results.get()
            pymol.cmd.distance("cys-cys", "Tetracysteine_model and chain " + line_text.split(",")[0].split()[1] +" and resi " + line_text.split(",")[0].split()[2] + " and name SG", "Tetracysteine_model and chain " + line_text.split(",")[1].split()[0] +" and resi " + line_text.split(",")[1].split()[1]  + " and name SG")
            pymol.cmd.distance("cys-cys", "Tetracysteine_model and chain " + line_text.split(",")[0].split()[1] +" and resi " + line_text.split(",")[0].split()[2] + " and name SG", "Tetracysteine_model and chain " + line_text.split(",")[2].split()[0] +" and resi " + line_text.split(",")[2].split()[1]  + " and name SG")
            pymol.cmd.distance("cys-cys", "Tetracysteine_model and chain " + line_text.split(",")[1].split()[0] +" and resi " + line_text.split(",")[1].split()[1] + " and name SG", "Tetracysteine_model and chain " + line_text.split(",")[3].split()[0] +" and resi " + line_text.split(",")[3].split()[1]  + " and name SG")
            pymol.cmd.distance("cys-cys", "Tetracysteine_model and chain " + line_text.split(",")[2].split()[0] +" and resi " + line_text.split(",")[2].split()[1] + " and name SG", "Tetracysteine_model and chain " + line_text.split(",")[3].split()[0] +" and resi " + line_text.split(",")[3].split()[1]  + " and name SG")
        else:
            pymol.cmd.delete("cys-cys")
            print "don't display measurements"

    def selectResultLine(self,event):
        selected_line= self.results_listbox.curselection()
        if selected_line[0] != 0:
            #selected a result
            selected_search = self.selected_results.get()
            selected_result = results_dict[selected_search]
            line_text = self.results_listbox.get(Tkinter.ANCHOR)
            selection_text = ""
            for i in range(0, 4):
                selection_text += results_dict[selected_search].getPymolObject()  + " and "
                if i == 0:
                    selection_text += "chain " + line_text.split(",")[i].split()[1]
                    selection_text += " and resi " + line_text.split(",")[i].split()[2]
                else:
                    selection_text += "chain " + line_text.split(",")[i].split()[0]
                    selection_text += " and resi " + line_text.split(",")[i].split()[1]
                if i != 3:
                    selection_text += " + "
            print selection_text
            pymol.cmd.select(selection_text)
            pymol.cmd.zoom(selection_text)

            #display nearby sidechains
            self.displayCys()
            self.displayMeasurements()
            self.displayPolar()
            self.displayNeg()
            self.displayPos()
            self.displayTrp()
            self.displayAromatic()


    def deleteSelection(self):
        if "Residues" not in self.results_listbox.get(Tkinter.ANCHOR):
            self.results_listbox.delete(Tkinter.ANCHOR)

    def compareLists(self):
        if self.list1_var.get() == "" or self.list2_var.get() == "":
            d = SearchErrorDialog(self.parent,"Please select lists to compare")
            self.window.wait_window(d.top)
        elif self.list1_var.get() == self.list2_var.get():
             d = SearchErrorDialog(self.parent,"Please select two different lists to compare")
             self.window.wait_window(d.top)
        else:
            compare = ListCompare(self.parent, self.list1_var.get(), self.list2_var.get())

    #fills list box will call other fillListBox methods if alternative sorting is being used
    def fillListBox(self, results_dict, b_dict, sasa_dict):
        #get the currently selected search
        selected_search = self.selected_results.get()
        if selected_search == "":
            return

        aa_list = results_dict[selected_search].getAAList()
        results_list = results_dict[selected_search].getResultsList()
        self.working_results_list = results_list
        chain_list = results_dict[selected_search].getChainsList()

        if self.sort_by_var.get() == 2:
            print "sort by sasa"
            self.fillListBoxSASAsort(results_dict, b_dict, sasa_dict)
            return ""
        elif self.sort_by_var.get() == 3:
            self.fillListBoxBsort(results_dict, b_dict, sasa_dict)
            print "sort by b factor"
            return ""

        #results = dist_matrix_score, index, rotamer_indices
        #convert results amino acids to chain and AA number
        print "sorting by score"
        i = 1
        print str(i)
        for result in results_list:
            text_line = ""
            if len(str(i)) == 1:
                text_line = str(i) + "               "
            elif len(str(i)) == 2:
                text_line = str(i) + "              "
            elif len(str(i)) == 3:
                text_line = str(i) + "             "
            elif len(str(i)) == 4:
                text_line = str(i) + "            "
            else:
                text_line = str(i) + "           "
            i += 1
            #print result
            score = result[0]
            aa_indices = result[1]
            site_aa_numbers =  []
            for index in aa_indices:
                chain_num = chain_list[index]
                aa_num = aa_list[index]
                if len(str(aa_num)) == 3:
                    text_line = text_line + str(chain_num) + " " + str(aa_num) + ", "
                elif len(str(aa_num)) == 2:
                    text_line = text_line + str(chain_num) + "   " + str(aa_num) + ", "
                else:
                    text_line = text_line + str(chain_num) + "     " + str(aa_num) + ", "
                aa_item =str(chain_num) + "," + str(aa_num)
                site_aa_numbers.append(str(aa_item))
            text_line = text_line + "\t            " + (str(score))[0:5]
            selected_obj = results_dict[selected_search].getPymolObject()
            text_line = text_line + ",       " + str(getAverageSasa(site_aa_numbers, sasa_dict))
            text_line =  text_line + ",    " + str(getAverageBFactor(site_aa_numbers, b_dict))
            self.results_listbox.insert(Tkinter.END, text_line)

        return ""

    #fill the listbox and sort by SASA
    def fillListBoxSASAsort(self, results_dict, b_dict, sasa_dict):
        sasa_and_results = []
        selected_search = self.selected_results.get()
        aa_list = results_dict[selected_search].getAAList()
        results_list = results_dict[selected_search].getResultsList()
        self.working_results_list = results_list
        chain_list = results_dict[selected_search].getChainsList()
        for result in results_list:
             text_line = ""
             #print result
             score = result[0]
             aa_indices = result[1]
             site_aa_numbers =  []
             for index in aa_indices:
                chain_num = chain_list[index]
                aa_num = aa_list[index]
                aa_item =str(chain_num) + "," + str(aa_num)
                site_aa_numbers.append(str(aa_item))
             sasa = getAverageSasa(site_aa_numbers, sasa_dict)
             sasa_and_results.append((sasa, result))

        #sort the it by sasa
        sorted_sasa_and_results = sorted(sasa_and_results, key=operator.itemgetter(0), reverse=True)
        self.results_listbox.delete(1, Tkinter.END)
        for sasa_result in sorted_sasa_and_results:
            result = sasa_result[1]
            text_line = ""
            score = result[0]
            aa_indices = result[1]
            site_aa_numbers =  []
            for index in aa_indices:
                chain_num = chain_list[index]
                aa_num = aa_list[index]
                if len(str(aa_num)) == 3:
                    text_line = text_line + str(chain_num) + " " + str(aa_num) + ", "
                elif len(str(aa_num)) == 2:
                    text_line = text_line + str(chain_num) + "   " + str(aa_num) + ", "
                else:
                    text_line = text_line + str(chain_num) + "     " + str(aa_num) + ", "
                aa_item =str(chain_num) + "," + str(aa_num)
                site_aa_numbers.append(str(aa_item))
            text_line = text_line + "\t            " + (str(score))[0:5]
            selected_obj = results_dict[selected_search].getPymolObject()
            text_line = text_line + ",       " + str(sasa_result[0])
            text_line =  text_line + ",    " + str(getAverageBFactor(site_aa_numbers, b_dict))
            self.results_listbox.insert(Tkinter.END, text_line)
        return ""


    #fill the listbox and sort by B factor
    def fillListBoxBsort(self, results_dict, b_dict, sasa_dict):
        b_and_results = []
        selected_search = self.selected_results.get()
        aa_list = results_dict[selected_search].getAAList()
        results_list = results_dict[selected_search].getResultsList()
        self.working_results_list = results_list
        chain_list = results_dict[selected_search].getChainsList()
        for result in results_list:
             score = result[0]
             aa_indices = result[1]
             site_aa_numbers =  []
             for index in aa_indices:
                chain_num = chain_list[index]
                aa_num = aa_list[index]
                aa_item =str(chain_num) + "," + str(aa_num)
                site_aa_numbers.append(str(aa_item))
             b_fact = getAverageBFactor(site_aa_numbers, b_dict)
             b_and_results.append((b_fact, result))
        sorted_b_and_result = sorted(b_and_results, key=operator.itemgetter(0), reverse=True)
        self.results_listbox.delete(1, Tkinter.END)
        for b_and_result in sorted_b_and_result:
            restult = b_and_result[1]
            text_line = ""
            score = result[0]
            aa_indices = result[1]
            site_aa_numbers =  []
            for index in aa_indices:
                chain_num = chain_list[index]
                aa_num = aa_list[index]
                if len(str(aa_num)) == 3:
                    text_line = text_line + str(chain_num) + " " + str(aa_num) + ", "
                elif len(str(aa_num)) == 2:
                    text_line = text_line + str(chain_num) + "   " + str(aa_num) + ", "
                else:
                    text_line = text_line + str(chain_num) + "     " + str(aa_num) + ", "
                aa_item =str(chain_num) + "," + str(aa_num)
                site_aa_numbers.append(str(aa_item))
            text_line = text_line + "\t            " + (str(score))[0:5]
            selected_obj = results_dict[selected_search].getPymolObject()
            text_line = text_line + ",       " + str(getAverageSasa(site_aa_numbers, sasa_dict))
            text_line =  text_line + ",    " + str(b_and_result[0])
            self.results_listbox.insert(Tkinter.END, text_line)
        return ""


    def __init__(self,app,search_results_names, search_results, b_dict, sasa_dict):
        self.parent = app.root
        self.search_results_names = search_results_names
        window =Tkinter.Tk()
        self.window = window
        self.open_opt = options = {}
        options['defaultextension'] = '.txt'
        options['filetypes'] = [('text files', '.txt')]
        options['initialdir'] = 'C:\\'
        options['parent'] = self.window
        options['title'] = 'Open Search Results File'
        self.save_opt = s_options = {}
        s_options['defaultextension'] = '.txt'
        s_options['filetypes'] = [('text files', '.txt')]
        s_options['initialdir'] = 'C:\\'
        s_options['parent'] = self.window
        s_options['title'] = 'Save Search Results File'
        window.title("ReAsH Search Results")


        result_options_frame = Tkinter.Frame(window,relief=Tkinter.RAISED, borderwidth=1)
        result_options_frame.pack(side=Tkinter.LEFT,padx=1,pady=1)

        results_title_frame = Tkinter.Frame(result_options_frame)
        results_title_frame.pack(fill=Tkinter.BOTH, expand=1,pady=4)
        results_title_label = Tkinter.Label(results_title_frame, text="Search Results")
        results_title_label.pack(side=Tkinter.LEFT)

        result_choice_frame = Tkinter.Frame(result_options_frame)
        result_choice_frame.pack(fill=Tkinter.BOTH, expand=1)
        display_results_label = Tkinter.Label(result_choice_frame, text="Display results for: ")
        display_results_label.pack(side=Tkinter.LEFT)

        self.selected_results = Tkinter.StringVar(window)
        self.display_results_dropdown = Tkinter.OptionMenu(result_choice_frame, self.selected_results, "") #add a command to this so when a different one is selected results show
        if len(search_results_names) > 0:
            self.selected_results.set(search_results_names[0])
            self.display_results_dropdown = apply(Tkinter.OptionMenu, (result_choice_frame, self.selected_results) + tuple(search_results_names))
        self.display_results_dropdown.bind("<Button-1>", self.resultsDropDownEvents)
        self.display_results_dropdown.config(width=10)
        self.display_results_dropdown.pack(side=Tkinter.LEFT)
        self.open_results_from_file = Tkinter.Button(result_choice_frame, text="Load from file",command=self.loadFromFile)
        self.open_results_from_file.pack(side=Tkinter.LEFT)

        spacer_1_frame = Tkinter.Frame(result_options_frame)
        spacer_1_frame.pack(fill=Tkinter.BOTH, expand=1,pady=7)

        filter_results_title_frame = Tkinter.Frame(result_options_frame)
        filter_results_title_frame.pack(fill=Tkinter.BOTH, expand=1,pady=2)
        filter_results_title = Tkinter.Label(filter_results_title_frame, text="Filter Results By:")
        filter_results_title.pack(side=Tkinter.LEFT)

        require_aromatic_frame = Tkinter.Frame(result_options_frame)
        require_aromatic_frame.pack(fill=Tkinter.BOTH, expand=1)
        self.require_aromatic = Tkinter.IntVar(window)
        self.require_aromatic_check = Tkinter.Checkbutton(require_aromatic_frame, text="Require Nearby Phe          ", variable=self.require_aromatic)
        self.require_aromatic_check.bind("<ButtonRelease-1>", self.filterAromatic)
        self.require_aromatic_check.pack(side=Tkinter.LEFT,padx=20)

        self.disallow_trp = Tkinter.IntVar(window)
        self.disallow_trp_check = Tkinter.Checkbutton(require_aromatic_frame,text="Don't Allow Nearby Trp or Tyr", variable=self.disallow_trp)
        self.disallow_trp_check.bind("<ButtonRelease-1>", self.filterTrp)
        self.disallow_trp_check.pack(side=Tkinter.LEFT,padx=20)

        require_pos_frame = Tkinter.Frame(result_options_frame)
        require_pos_frame.pack(fill=Tkinter.BOTH, expand=1)
        self.require_pos = Tkinter.IntVar(window)
        self.require_pos_check = Tkinter.Checkbutton(require_pos_frame, text="Require Nearby Positive Residue", variable=self.require_pos)
        self.require_pos_check.bind("<ButtonRelease-1>", self.filterPos)
        self.require_pos_check.pack(side=Tkinter.LEFT,padx=20)

        self.disallow_neg = Tkinter.IntVar(window)
        self.disallow_neg_check = Tkinter.Checkbutton(require_pos_frame, text="Don't Allow Nearby Negative Residue", variable=self.disallow_neg)
        self.disallow_neg_check.bind("<ButtonRelease-1>", self.filterNeg)
        self.disallow_neg_check.pack(side=Tkinter.LEFT,padx=20)

        require_all_allowed_rots_frame = Tkinter.Frame(result_options_frame)
        require_all_allowed_rots_frame.pack(fill=Tkinter.BOTH,expand=1)
        self.require_allowed_rots = Tkinter.IntVar(window)
        self.require_allowed_rots_check = Tkinter.Checkbutton(require_all_allowed_rots_frame, text="Require Allowed Rotamers", variable=self.require_allowed_rots)
        self.require_allowed_rots_check.bind("<ButtonRelease-1>", self.filterRotamer)
        self.require_allowed_rots_check.pack(side=Tkinter.LEFT,padx=20)

        spacer_2_frame = Tkinter.Frame(result_options_frame)
        spacer_2_frame.pack(fill=Tkinter.BOTH, expand=1,pady=2)

        sorting_opt_label_frame = Tkinter.Frame(result_options_frame)
        sorting_opt_label_frame.pack(fill=Tkinter.BOTH,expand=1,pady=2)
        sorting_opt_label =Tkinter.Label(sorting_opt_label_frame,text="Sort Results By: ")
        sorting_opt_label.pack(side=Tkinter.LEFT)

        self.sort_by_var = Tkinter.IntVar(result_options_frame)
        self.sort_by_var.set(1)
        sort_score_frame = Tkinter.Frame(result_options_frame)
        sort_score_frame.pack(fill=Tkinter.BOTH,expand=1)
        self.sort_score_option = Tkinter.Radiobutton(sort_score_frame, text="Sort By Score", variable=self.sort_by_var, value=1,command=self.changeSort)
       # self.sort_score_option = Tkinter.Radiobutton(sort_score_frame, text="Sort By Score", variable=self.sort_by_var, value=1)
        self.sort_score_option.pack(side=Tkinter.LEFT, padx=20)
        self.sort_score_option.select()

        sort_SASA_frame = Tkinter.Frame(result_options_frame)
        sort_SASA_frame.pack(fill=Tkinter.BOTH,expand=1)
        self.sort_SASA_option = Tkinter.Radiobutton(sort_SASA_frame, text="Sort By SASA", variable=self.sort_by_var,value=2,command=self.changeSort)
        self.sort_SASA_option.pack(side=Tkinter.LEFT, padx=20)

        sort_B_frame = Tkinter.Frame(result_options_frame)
        sort_B_frame.pack(fill=Tkinter.BOTH,expand=1)
        self.sort_B_option = Tkinter.Radiobutton(sort_B_frame, text ="Sort By B-factor", variable=self.sort_by_var, value=3,command=self.changeSort)
        self.sort_B_option.pack(side=Tkinter.LEFT,padx=20)



        spacer_3_frame = Tkinter.Frame(result_options_frame)
        spacer_3_frame.pack(fill=Tkinter.BOTH, expand=1,pady=2)

        display_options_label_frame = Tkinter.Frame(result_options_frame)
        display_options_label_frame.pack(fill=Tkinter.BOTH, expand=1,pady=2)
        display_options_label = Tkinter.Label(display_options_label_frame,text="Site Display Options:")
        display_options_label.pack(side=Tkinter.LEFT)

        display_options_row1_frame = Tkinter.Frame(result_options_frame)
        display_options_row1_frame.pack(fill=Tkinter.BOTH,expand=1)
        self.aromatic_display = Tkinter.IntVar(result_options_frame)
        self.aromatic_display_option = Tkinter.Checkbutton(display_options_row1_frame, text="Display Nearby Phe",variable=self.aromatic_display, command=self.displayAromatic)
        self.aromatic_display_option.pack(side=Tkinter.LEFT,padx=20)
        self.trp_display = Tkinter.IntVar(result_options_frame)
        self.trp_display_option = Tkinter.Checkbutton(display_options_row1_frame,text="Display Nearby Trp and Tyr",variable=self.trp_display, command=self.displayTrp)
        self.trp_display_option.pack(side=Tkinter.LEFT,padx=37)

        display_options_row2_frame = Tkinter.Frame(result_options_frame)
        display_options_row2_frame.pack(fill=Tkinter.BOTH,expand=1)
        self.pos_display = Tkinter.IntVar(result_options_frame)
        self.pos_display_option = Tkinter.Checkbutton(display_options_row2_frame,text="Display Nearby Positive Residues", variable=self.pos_display, command=self.displayPos)
        self.pos_display_option.pack(side=Tkinter.LEFT,padx=20)
        self.neg_display = Tkinter.IntVar(result_options_frame)
        self.neg_display_option = Tkinter.Checkbutton(display_options_row2_frame,text="Display Nearby Negative Residues", variable=self.neg_display, command=self.displayNeg)
        self.neg_display_option.pack(side=Tkinter.LEFT,padx=10)
        #self.neg_display_option.bind("<ButtonRelease-1>",self.displayNeg)

        display_options_row3_frame = Tkinter.Frame(result_options_frame)
        display_options_row3_frame.pack(fill=Tkinter.BOTH,expand=1)
        self.polar_display =Tkinter.IntVar(result_options_frame)
        self.polar_display_option = Tkinter.Checkbutton(display_options_row3_frame,text="Display Nearby Polar Residues",variable=self.polar_display, command=self.displayPolar)
        self.polar_display_option.pack(side=Tkinter.LEFT, padx=20)
        #self.polar_display_option.bind("<ButtonRelease-1>",self.displayPolar)
        self.ReAsH_display = Tkinter.IntVar()
        self.ReAsH_display_option = Tkinter.Checkbutton(display_options_row3_frame,text="Display ReAsH",variable=self.ReAsH_display)
        self.ReAsH_display_option.pack(side=Tkinter.LEFT,padx=23)
        self.ReAsH_display_option.bind("<ButtonRelease-1>",self.displayReAsH)

        display_options_row4_frame = Tkinter.Frame(result_options_frame)
        display_options_row4_frame.pack(fill=Tkinter.BOTH,expand=1)
        self.display_side_chains = Tkinter.IntVar(result_options_frame)
        self.display_side_chains_option = Tkinter.Checkbutton(display_options_row4_frame, text="Display Cys Sidechains",variable=self.display_side_chains,command=self.displayCys)
        self.display_side_chains_option.pack(side=Tkinter.LEFT,padx=20)
        #self.display_side_chains_option.bind("<ButtonRelease-1>",self.displayCys)
        self.display_measurements = Tkinter.IntVar(result_options_frame)
        self.display_measurements_options = Tkinter.Checkbutton(display_options_row4_frame,text="Display Sulfur-Sulfur Distances",variable=self.display_measurements,command=self.displayMeasurements)
        self.display_measurements_options.pack(side=Tkinter.LEFT,padx=61)
        #self.display_measurements_options.bind("<ButtonRelease-1>",self.displayMeasurements)

        results_and_compare_frame =Tkinter.Frame(window,relief=Tkinter.RAISED, borderwidth=1)
        results_and_compare_frame.pack(side=Tkinter.LEFT,padx=5,pady=1)

        results_frame = Tkinter.Frame(results_and_compare_frame)
        results_frame.pack(side=Tkinter.TOP)

        results_label_frame =Tkinter.Frame(results_frame)
        results_label_frame.pack(fill=Tkinter.BOTH,expand=1)
        results_label = Tkinter.Label(results_label_frame,text="List of Results:")
        results_label.pack(side=Tkinter.LEFT,padx=10)

        results_list_frame = Tkinter.Frame(results_frame)
        results_list_frame.pack(fill=Tkinter.BOTH,expand=1)
        results_scrollbar = Tkinter.Scrollbar(results_list_frame,orient=Tkinter.VERTICAL)
        #results_scrollbar.config(height=50)
        self.results_listbox = Tkinter.Listbox(results_list_frame,yscrollcommand=results_scrollbar)
        #self.results_listbox.config(yscrollcommand=results_scrollbar.set)
        results_scrollbar.config(command=self.results_listbox.yview)
        self.results_listbox.bind("<ButtonRelease-1>", self.selectResultLine)
        results_scrollbar.pack(side=Tkinter.RIGHT,padx=2)
        self.results_listbox.pack(side=Tkinter.LEFT)
        self.results_listbox.insert(Tkinter.END, "Result Number   Residues                                         \t Score      \t SASA                    \t B-score")
       # self.results_listbox.insert(Tkinter.END, "blah blah blah")
        self.results_listbox.config(width =80,height=15)

        results_buttons_frame = Tkinter.Frame(results_frame)
        results_buttons_frame.pack(fill=Tkinter.BOTH,expand =1,pady=10)
        self.delete_selection_button = Tkinter.Button(results_buttons_frame,text="Delete Selection",command=self.deleteSelection)
        self.delete_selection_button.pack(side=Tkinter.LEFT,padx=30)
        self.restore_default_button = Tkinter.Button(results_buttons_frame,text="Restore List to Original",command=self.restoreListToDefualt)
        self.restore_default_button.pack(side=Tkinter.LEFT,padx=30)
        self.save_button =Tkinter.Button(results_buttons_frame,text="Save List",command=self.saveToFile)
        self.save_button.pack(side=Tkinter.LEFT,padx=30)


        compare_frame = Tkinter.Frame(results_and_compare_frame)
        compare_frame.pack(fill=Tkinter.BOTH,expand=1,pady=10)
        compare_label_frame = Tkinter.Frame(compare_frame)
        compare_label_frame.pack(fill=Tkinter.BOTH,expand=1)
        compare_label = Tkinter.Label(compare_label_frame, text="Compare Two Search Lists:")
        compare_label.pack(side=Tkinter.LEFT)

        compare_dropdown_frame = Tkinter.Frame(compare_frame)
        compare_dropdown_frame.pack(fill=Tkinter.BOTH,expand=1)
        list1_label = Tkinter.Label(compare_dropdown_frame,text="List 1")
        list1_label.pack(side=Tkinter.LEFT,padx=10)
        self.list1_var = Tkinter.StringVar(window)
        self.list1_dropdown = Tkinter.OptionMenu(compare_dropdown_frame, self.list1_var, "")
        if len(search_results_names) > 0:
            self.list1_var.set(search_results_names[0])
            self.list1_dropdown = apply(Tkinter.OptionMenu, (compare_dropdown_frame, self.list1_var) + tuple(search_results_names))
        self.list1_dropdown.config(width=15)
        self.list1_dropdown.pack(side=Tkinter.LEFT)
        list2_label = Tkinter.Label(compare_dropdown_frame,text="List 2")
        list2_label.pack(side=Tkinter.LEFT,padx=10)
        self.list2_var = Tkinter.StringVar(window)
        self.list2_dropdown = Tkinter.OptionMenu(compare_dropdown_frame, self.list2_var, "")
        if len(search_results_names) > 0:
            self.list2_var.set(search_results_names[0])
            self.list2_dropdown = apply(Tkinter.OptionMenu, (compare_dropdown_frame, self.list2_var) + tuple(search_results_names))
        self.list2_dropdown.config(width=15)
        self.list2_dropdown.pack(side=Tkinter.LEFT)

        compare_button_frame = Tkinter.Frame(compare_frame)
        compare_button_frame.pack(fill=Tkinter.BOTH,expand=1)
        compare_button = Tkinter.Button(compare_button_frame, text="Compare These Lists",command=self.compareLists)
        compare_button.pack()
        self.working_results_list = []
        self.search_results = search_results
        self.b_dict = b_dict
        self.sasa_dict = sasa_dict
        self.fillListBox(search_results, b_dict, sasa_dict)

class ListCompare:
    def __init__(self,  app, list1_var, list2_var):
        print list1_var
        print list2_var
        self.parent = app.root
        window =Tkinter.Tk()
        window.title("Compare Two Searches")

        list_frames = Tkinter.Frame(window)
        list_frames.pack(fill=Tkinter.BOTH,expand=1)

        list_1_frame = Tkinter.Frame(list_frames,relief=Tkinter.RAISED, borderwidth=1)
        list_1_frame.pack(side=Tkinter.LEFT,padx=2,pady=2)

        list_1_dropdown_frame = Tkinter.Frame(list_1_frame)
        list_1_dropdown_frame.pack(fill=Tkinter.BOTH,expand=1)
        list_1_dropdown_label = Tkinter.Label(list_1_dropdown_frame, text="List 1: ")
        list_1_dropdown_label.pack(side=Tkinter.LEFT)
        self.list_1 = Tkinter.StringVar(window)
        self.list_1_option = Tkinter.OptionMenu(list_1_dropdown_frame, self.list_1, "")
        self.list_1_option.config(width=20)
        self.list_1_option.pack(side=Tkinter.LEFT)

        list_1_box_frame = Tkinter.Frame(list_1_frame)
        list_1_box_frame.pack(fill=Tkinter.BOTH,expand=1)
        list1_scrollbar = Tkinter.Scrollbar(list_1_box_frame,orient=Tkinter.VERTICAL)
        self.listbox1 = Tkinter.Listbox(list_1_box_frame,yscrollcommand=list1_scrollbar)
        list1_scrollbar.pack(side=Tkinter.RIGHT,padx=2)
        self.listbox1.pack(side=Tkinter.LEFT)
        self.listbox1.insert(Tkinter.END, "Result Number   Residues                  \t Score      \t SASA                    \t B-score")
        self.listbox1.config(width =80,height=15)

        list_1_sort_label_frame = Tkinter.Frame(list_1_frame)
        list_1_sort_label_frame.pack(fill=Tkinter.BOTH,expand=1)
        list_1_sort_label = Tkinter.Label(list_1_sort_label_frame,text="Sort Results By: ")
        list_1_sort_label.pack(side=Tkinter.LEFT)
        self.list1_sort_val = Tkinter.IntVar()
        list_1_sort_score_frame = Tkinter.Frame(list_1_frame)
        list_1_sort_score_frame.pack(fill=Tkinter.BOTH, expand=1)
        self.list_1_sort_score=Tkinter.Radiobutton(list_1_sort_score_frame, text="Sort By Score", variable=self.list1_sort_val, value=1)
        self.list_1_sort_score.select()
        self.list_1_sort_score.pack(side=Tkinter.LEFT, padx=10)
        list_1_sort_sasa_frame = Tkinter.Frame(list_1_frame)
        list_1_sort_sasa_frame.pack(fill=Tkinter.BOTH,expand=1)
        self.list_1_sort_sasa =Tkinter.Radiobutton(list_1_sort_sasa_frame, text = "Sort By SASA", variable=self.list1_sort_val,value=2)
        self.list_1_sort_sasa.pack(side=Tkinter.LEFT,padx=10)
        list_1_sort_b_factor_frame = Tkinter.Frame(list_1_frame)
        list_1_sort_b_factor_frame.pack(fill=Tkinter.BOTH,expand=1)
        self.list_1_sort_b_factor = Tkinter.Radiobutton(list_1_sort_b_factor_frame, text="Sort By B-factor",variable=self.list1_sort_val,value=3)
        self.list_1_sort_b_factor.pack(side=Tkinter.LEFT,padx=10)

        list_2_frame = Tkinter.Frame(list_frames,relief=Tkinter.RAISED, borderwidth=1)
        list_2_frame.pack(fill=Tkinter.BOTH,padx=2,pady=2)

        list_2_dropdown_frame = Tkinter.Frame(list_2_frame)
        list_2_dropdown_frame.pack(fill=Tkinter.BOTH,expand=1)
        list_2_dropdown_label = Tkinter.Label(list_2_dropdown_frame, text="List 2: ")
        list_2_dropdown_label.pack(side=Tkinter.LEFT)
        self.list_2 = Tkinter.StringVar(window)
        self.list_2_option = Tkinter.OptionMenu(list_2_dropdown_frame, self.list_2, "")
        self.list_2_option.config(width=20)
        self.list_2_option.pack(side=Tkinter.LEFT)

        list_2_box_frame = Tkinter.Frame(list_2_frame)
        list_2_box_frame.pack(fill=Tkinter.BOTH,expand=1)
        list2_scrollbar = Tkinter.Scrollbar(list_2_box_frame,orient=Tkinter.VERTICAL)
        self.listbox2 = Tkinter.Listbox(list_2_box_frame,yscrollcommand=list2_scrollbar)
        list2_scrollbar.pack(side=Tkinter.RIGHT,padx=2)
        self.listbox2.pack(side=Tkinter.LEFT)
        self.listbox2.insert(Tkinter.END, "Result Number   Residues                  \t Score      \t SASA                    \t B-score")
        self.listbox2.config(width =80,height=15)

        list_2_sort_label_frame = Tkinter.Frame(list_2_frame)
        list_2_sort_label_frame.pack(fill=Tkinter.BOTH,expand=1)
        list_2_sort_label = Tkinter.Label(list_2_sort_label_frame,text="Sort Results By: ")
        list_2_sort_label.pack(side=Tkinter.LEFT)
        self.list2_sort_val = Tkinter.IntVar()
        list_2_sort_score_frame = Tkinter.Frame(list_2_frame)
        list_2_sort_score_frame.pack(fill=Tkinter.BOTH, expand=1)
        self.list_2_sort_score=Tkinter.Radiobutton(list_2_sort_score_frame, text="Sort By Score", variable=self.list2_sort_val, value=1)
        self.list_2_sort_score.select()
        self.list_2_sort_score.pack(side=Tkinter.LEFT, padx=10)
        list_2_sort_sasa_frame = Tkinter.Frame(list_2_frame)
        list_2_sort_sasa_frame.pack(fill=Tkinter.BOTH,expand=1)
        self.list_2_sort_sasa =Tkinter.Radiobutton(list_2_sort_sasa_frame, text = "Sort By SASA", variable=self.list2_sort_val,value=2)
        self.list_2_sort_sasa.pack(side=Tkinter.LEFT,padx=10)
        list_2_sort_b_factor_frame = Tkinter.Frame(list_2_frame)
        list_2_sort_b_factor_frame.pack(fill=Tkinter.BOTH,expand=1)
        self.list_2_sort_b_factor = Tkinter.Radiobutton(list_2_sort_b_factor_frame, text="Sort By B-factor",variable=self.list2_sort_val,value=3)
        self.list_2_sort_b_factor.pack(side=Tkinter.LEFT,padx=10)

        display_frame = Tkinter.Frame(window,relief=Tkinter.RAISED,borderwidth=1)
        display_frame.pack(fill=Tkinter.BOTH,expand=1,padx=2,pady=2)
        display_label_frame = Tkinter.Frame(display_frame)
        display_label_frame.pack(fill=Tkinter.BOTH,expand=1)
        display_label = Tkinter.Label(display_label_frame, text="Display Options:")
        display_label.pack(side=Tkinter.LEFT)

        display_options_line1_frame = Tkinter.Frame(display_frame)
        display_options_line1_frame.pack(fill=Tkinter.BOTH,expand=1)
        self.display_aromatic = Tkinter.IntVar()
        self.display_aromatic_option = Tkinter.Checkbutton(display_options_line1_frame, text="Display Nearby Phe and Trp", variable=self.display_aromatic)
        self.display_aromatic_option.pack(side=Tkinter.LEFT,padx=10)
        self.display_trp = Tkinter.IntVar()
        self.display_trp_option = Tkinter.Checkbutton(display_options_line1_frame,text="Display Nearby Trp   ", variable=self.display_trp)
        self.display_trp_option.pack(side=Tkinter.LEFT,padx=10)
        self.display_pos = Tkinter.IntVar()
        self.display_pos_option = Tkinter.Checkbutton(display_options_line1_frame,text="Display Nearby Positive Residues",variable=self.display_pos)
        self.display_pos_option.pack(side=Tkinter.LEFT,padx=10)
        self.display_neg = Tkinter.IntVar()
        self.display_neg_option = Tkinter.Checkbutton(display_options_line1_frame,text="Display Nearby Negative Residues",variable=self.display_neg)
        self.display_neg_option.pack(side=Tkinter.LEFT,padx=10)
        self.display_polar = Tkinter.IntVar()
        self.display_polar_option = Tkinter.Checkbutton(display_options_line1_frame,text="Display Nearby Polar Residues",variable=self.display_polar)
        self.display_polar_option.pack(side=Tkinter.LEFT, padx=10)

        display_options_line2_frame = Tkinter.Frame(display_frame)
        display_options_line2_frame.pack(fill=Tkinter.BOTH, expand =1)
        self.display_reash = Tkinter.IntVar()
        self.display_reash_option = Tkinter.Checkbutton(display_options_line2_frame,text="Display ReAsH                       ",variable=self.display_reash)
        self.display_reash_option.pack(side=Tkinter.LEFT,padx=10)
        self.display_sidechains = Tkinter.IntVar()
        self.display_sidechains_option = Tkinter.Checkbutton(display_options_line2_frame,text="Display Cys Sidechains",variable=self.display_sidechains)
        self.display_sidechains_option.pack(side=Tkinter.LEFT,padx=10)
        self.display_measure = Tkinter.IntVar()
        self.display_measure_option = Tkinter.Checkbutton(display_options_line2_frame, text="Display Sulfur-Sulfur Distances    ",variable=self.display_measure)
        self.display_measure_option.pack(side=Tkinter.LEFT,padx=0)
        self.display_comparison = Tkinter.IntVar()
        self.display_comparison_option = Tkinter.Checkbutton(display_options_line2_frame,text="Display Site in Other Object")
        self.display_comparison_option.pack(side=Tkinter.LEFT,padx=20)

        self.compareList(list1_var, list2_var)

    def compareLists(self, list1_var, list2_var):
        aa_list1 = results_dict[list1_var].getAAList()
        aa_list2 = results_dict[list2_var].getAAList()
        results_list1 = results_dict[list1_var].getResultsList()
        results_list2 = results_dict[list2_var].getResultsList()
        chain_list1 = results_dict[list1_var].getChainsList()
        chain_list2 = results_dict[list2_var].getChainsList()

        r1_indices_list = []
        r2_indices_list = []
        r1_unique_indices = []
        r2_unique_indices = []
        for i in results_list1[1]:
            r1_indices_list.append(i)
        for i in results_list2[1]:
            r2_indices_list.append(i)
        for i in r1_indices_list:
            if i not in r2_indices_list:
                r1_indices_list.append(i)
        for i in r2_indices_list:
            if i not in r1_indices_list:
                r2_unique_indices.append(i)



        return






#if __name__ == '__main__':
 #   main()

