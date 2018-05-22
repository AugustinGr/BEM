# -*- coding: utf-8 -*-

"""
Created on Thu Apr 26 15:46:00 2018

@author: augus
"""

import csv, math, numpy, re
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


#WT='NREL_5MW'
#WT='DTU_10MW'

#We are looking for the data stored on 2nd lines (0th, 1st, 2nd)

def get_thickness_data(WT):
    fichier = 'WT_data/' + WT + '/Rwt.pro'
    f = open(fichier, "rb")
    fichier_csv = csv.reader(f, delimiter=" ")
    tab = list(fichier_csv)
    #print tab[2]
    L=[]
    for c in tab[2]:
        x = re.findall(r"[-+]?\d*\.\d+|\d+", c)

        if len(x)==1:
            x = float(x[0])
            L.append(x)
    #print L

    Thickness=list(L)
    f.close()
    return Thickness

#print get_thickness_data(WT)