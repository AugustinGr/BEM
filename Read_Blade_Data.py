# -*- coding: utf-8 -*-

"""
Created on Thu May 1 13:58:00 2018

@author: augus
"""

import numpy as np
import csv, math, re
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from Read_Thickness_Data import get_thickness_data

#WT='NREL_5MW'
#WT='DTU_10MW'

########################################################################################################################
def Extract_to_CharList(WT):
    """
    Read, and extract Blade_Data.vda with a list of lines
    :param WT: names of wind turbine blade
    :return: List of Lines
    """
    fichier = 'WT_data/'+str(WT)+'/Blade_Data.vda'
    f = open(fichier, "rb")
    fichier_csv = csv.reader(f, delimiter=" ")
    tab = list(fichier_csv)
    f.close()
    return tab

"""
L=Extract_to_CharList(WT)
print L
"""

########################################################################################################################
def Look_for_Table_thickness(L):
    """
    Extract the data from the list of lines, the second table (wich concern the thickness)
    :param L: List of lines from the data files .vda
    :return: List of lines the second table related to the thickness
    """
    check=0
    new_L=[]
    for l in L:
        for i in range(len(l)):
            if 't/c' in l[i]:
                #return l
                check=1
            if 'Rwt.pro' in l[i]:
                check=0
        if check==1:
            new_L.append(l)
    return new_L

"""
L=Look_for_Table_thickness(L)
#print L
"""

########################################################################################################################
def create_table_of_floatting_values(L):
    """
    Convert the char value to floatting value
    :param L:
    :return: [[radius, chord, beta, t/c, yac/c],...,]
    """
    new_L=[]
    for l in L:
        new_l=[]
        for c in l:
            x = re.findall(r"[-+]?\d*\.\d+|\d+", c)
            if len(x) == 1:
                new_l = new_l + [float(x[0])]
        new_L=new_L+[new_l]
    return new_L

"""
L=create_table_of_floatting_values(L)
#print L
"""

########################################################################################################################
def List_to_Matrix(L):
    """
    Convert the list to array
    """
    lines=len(L)
    rows=len(L[0])
    M=np.zeros((lines,rows))
    for i in range(lines):
        for k in range(rows):
            M[i,k]=L[i][k]
    #add hub radius
    m=np.zeros((1,rows))
    for k in range(1,rows):
        m[0,k]=M[0,k]
    M=np.concatenate((m,M))

    return M

"""
M= List_to_Matrix(L)
#print M
"""

########################################################################################################################
def normalized_radius(M):
    """
    Normalized the radius r to r/R
    """
    Blade_radius=M[-1,0]
    M[:,0]=M[:,0]/Blade_radius
    return M

"""
M=normalized_radius(M)
#print M
"""

########################################################################################################################
def interpolate_thickness_function_of_radius(M):
    """
    1D-Interpolation giving: thickness = f(radius)
    """
    x=M[:,0]
    y=M[:,3]
    f = interp1d(x,y,kind='linear')
    return f

def interpolate_radius_function_of_thickness(M):
    """
    1D-Interpolation giving: thickness = f(radius)
    """
    y=M[:,0]
    x=M[:,3]
    f = interp1d(x,y,kind='linear')
    return f

"""
f= interpolate_thickness_function_of_radius(M)

X=np.arange(0.,1,0.001)
plt.plot(X,f(X))
plt.show()
"""

########################################################################################################################
def interpolate_chord_function_of_radius(M):
    """
    1D-Interpolation giving: chord = f(radius)
    """
    x=M[:,0]
    y=M[:,1]
    f = interp1d(x,y,kind='linear')
    return f

"""
f = interpolate_chord_function_of_radius(M)

X = np.arange(0., 1, 0.001)
plt.plot(X, f(X))
plt.show()
"""

def interpolate_beta_fuction_of_radius(M):
    """
    1D-Interpolation giving: beta = f(radius)
    """
    x = M[:, 0]
    y = M[:, 2]
    f = interp1d(x, y, kind='linear')
    return f
########################################################################################################################
########################################################################################################################
########################################################################################################################

def Interpolate_Blade_Thickness_Chord_function_of_radius(WT):
    """
    Main Function
    :param WT: Kind of Blade
    :return:
    """

    L = Extract_to_CharList(WT)
    L = Look_for_Table_thickness(L)
    L = create_table_of_floatting_values(L)
    M = List_to_Matrix(L)
    #M = normalized_radius(M)

    f_thickness = interpolate_thickness_function_of_radius(M)
    f_chord = interpolate_chord_function_of_radius(M)
    f_beta = interpolate_beta_fuction_of_radius(M)
    f_radius = interpolate_radius_function_of_thickness(M)

    return [f_thickness,f_chord, f_beta, f_radius]

"""
WT='NREL_5MW'
#WT='DTU_10MW'

F=Interpolate_Blade_Thickness_Chord_function_of_radius(WT)
f_thick=F[0]
f_chord=F[1]

X = np.arange(0., 1, 0.001)
plt.plot(X, f_thick(X))
plt.show()

plt.plot(X,f_chord(X))
plt.show()
"""

def extract_hub_radius_and_total_radius(WT):
    L = Extract_to_CharList(WT)
    L = Look_for_Table_thickness(L)
    L = create_table_of_floatting_values(L)
    R_root=L[0][0]
    R_total=L[-1][0]
    Thickness=get_thickness_data(WT)
    f=Interpolate_Blade_Thickness_Chord_function_of_radius(WT)[3]
    R_end_cylinder=f(Thickness[-2])
    return [R_root,R_total,R_end_cylinder]


"""
L=extract_hub_radius_and_total_radius('NREL_5MW')
print L
"""