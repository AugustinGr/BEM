# -*- coding: utf-8 -*-

"""
Created on Wed Apr 25 15:32:00 2018

@author: augus
"""
# Lecture des fichiers texte en CSV

# importer les modules
import csv, math, numpy, re
from scipy import interpolate
from scipy.interpolate import interp1d, interp2d, griddata, Rbf
#from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
from Read_Thickness_Data import get_thickness_data
from mpl_toolkits.mplot3d.axes3d import Axes3D

#WT='NREL_5MW'
#WT='DTU_10MW'

""" Function to Extract Aerofoil data """
def Make_AirFoils_Profile(tab,WT):
    AirFoils_Profile=[]
    for i in range(len(tab)):

        if 'DTU_10MW' in WT:
            if 'F' in tab[i][0] or 'Cyl' in tab[i][0]:
                #print tab[i][0]
                AirFoils_Profile.append(tab[i][0])
                #print AirFoils_Profile
        if 'NREL_5MW' in WT:
            if 'NACA' in tab[i][0] or 'Round' in tab[i][0] or 'DU' in tab[i][0] :
                #print tab[i][0]
                if len(tab[i][0])>1:
                    AirFoils_Profile.append(tab[i][0])
                #print AirFoils_Profile
    #print AirFoils_Profile
    #print len(AirFoils_Profile)
    return AirFoils_Profile

""" Function to Extract Aerofoil data """
def get_aerofoils_data(WT):
    """ Function to extract AeroCoeff data
    Input:
        WT: Name of the choosen Wind Turbine
        AirFoil_Profile: list of char containing names for all profiles for this Wind Turbine
    Output:
        Profiles_DATA: list of list of aerocoeff for each profiles
    """
    fichier = 'WT_data/'+WT+'/Rwt.pro'
    f = open(fichier, "rb")
    fichier_csv = csv.reader(f, delimiter=" ")
    tab = list(fichier_csv)
    ###########
    AirFoil_Profile=Make_AirFoils_Profile(tab,WT)
    ###########
    Profiles_DATA=[]
    k=0
    index=-1
    profil_data=[]

    for i in range(len(tab)):
        #print k
        for profil in AirFoil_Profile:
            if profil in tab[i][0] and k!=0:
                #print str(i) + "  out"
                #print tab[i]+['out']

                k=0
                Profiles_DATA.append(profil_data)#=Profiles_DATA+[profil_data]
        if k==1:
            #print 'i= '+str(i)
            L=[]
            for c in tab[i]:
                x = re.findall(r"[-+]?\d*\.\d+|\d+", c)
                #print x
                if x != []:
                    L = L + x
            profil_data=profil_data+[L]

        for profil in AirFoil_Profile:
            if profil in tab[i][0]:
                #print str(i) + "  in"
                #print tab[i]
                k=1
                index=index+1
                profil_data=[]
    Profiles_DATA.append(profil_data) #Add the last Profile Data

    #print Profiles_DATA[0][:][0]
    #print len(Profiles_DATA[:])
    f.close()
    return Profiles_DATA[:]

#########TEST######
#AirFoil_Profile=['NACA64','DU21','DU25','DU30','DU35','DU40','Round']
#WT='NREL_5MW'
#print len(get_aerofoils_data(WT,AirFoil_Profile))
########################################################################################################################
""" Function to convert char data from get_aerofoils_data to floatting values"""
def CharToFloat(AeroChar):
    """
    Function to convert char data from get_aerofoils_data to floatting values
    :param AeroChar: containing char
    :return: AeroChar: containing float
    """

    for i0 in range(len(AeroChar)):
        for i1 in range(len(AeroChar[i0])):
            for i2 in range(len(AeroChar[i0][i1])):
                AeroChar[i0][i1][i2]=float(AeroChar[i0][i1][i2])
    return AeroChar

##############TEST#######################
#AirFoil_Profile=['NACA64','DU21','DU25','DU30','DU35','DU40','Round']
#WT='NREL_5MW'
#print CharToFloat(get_aerofoils_data(WT,AirFoil_Profile))

########################################################################################################################
""" Function to organize the data"""
def Organize_data_V2(Aerofoils_Data,Thickness):
    """
    Create a aero matrix (number_of_profile*number_of_data_per_profil,4), (Thickness ,Alpha, Cl, Cd)

    :param Aerofoils_Data: Data under list format
    :param Thickness: list format
    :return: Matrix Aero
    """
    #determins the dim of data:
    L=[]
    for i in range(len(Aerofoils_Data)):
        for k in range(len(Aerofoils_Data[i])):
            Aerofoils_Data[i][k]=[Thickness[i]]+Aerofoils_Data[i][k]
            L.append(Aerofoils_Data[i][k])
    #print L
    #print len(L)
    New_Aerofoils_Data_mat=numpy.zeros((len(L),4))
    for i in range(len(L)):
        New_Aerofoils_Data_mat[i, 0] = L[i][0]   # Thickness
        New_Aerofoils_Data_mat[i, 1] = L[i][1]   # Angular
        New_Aerofoils_Data_mat[i, 2] = L[i][2]   # Cl
        New_Aerofoils_Data_mat[i, 3] = L[i][3]   # Cd
    #print New_Aerofoils_Data_mat
    return New_Aerofoils_Data_mat

def Organize_data_V3(Aerofoils_Data,Thickness):
    """
    Create a aero matrix (number_of_profile*number_of_data_per_profil,4), (Thickness ,Alpha, Cl, Cd)
    +
    Include the thickness = 0 data, i.e: Cl=Cd=0
    :param Aerofoils_Data: Data under list format
    :param Thickness: list format
    :return: Matrix Aero
    """
    L=[]
    for i in range(len(Aerofoils_Data)):
        for k in range(len(Aerofoils_Data[i])):
            Aerofoils_Data[i][k]=[Thickness[i]]+Aerofoils_Data[i][k]
            L.append(Aerofoils_Data[i][k])
    #print L
    #print len(L)
    New_Aerofoils_Data_mat=numpy.zeros((len(L),4))
    for i in range(len(L)):
        New_Aerofoils_Data_mat[i, 0] = L[i][0]   # Thickness
        New_Aerofoils_Data_mat[i, 1] = L[i][1]   # Angular
        New_Aerofoils_Data_mat[i, 2] = L[i][2]   # Cl
        New_Aerofoils_Data_mat[i, 3] = L[i][3]   # Cd

    # create a matrix wich represent the data at thickness=0.
    lenght_0_matrix=len(Aerofoils_Data[0])
    zeros_matrix=numpy.zeros((lenght_0_matrix,4))
    for k in range(lenght_0_matrix):
        zeros_matrix[k,1]=L[k][1]

    # Concatenate the two matrix
    New_Aerofoils_Data_mat=numpy.concatenate((zeros_matrix,New_Aerofoils_Data_mat))


    #print New_Aerofoils_Data_mat
    return New_Aerofoils_Data_mat

#############TEST#############
#WT='NREL_5MW'
#result=Organize_data_V2(CharToFloat(get_aerofoils_data(WT)),get_thickness_data(WT))
#print result
#print len(result)
#print len(result[0])

########################################################################################################################

def interpolation2d_aerofoil_data(Aero_Mat):
    """
        Function returning for a value of an angle the list (Cl,Cd) interpolate for this angle
        :param profile_index: Index of Profile Number (Exemple for NREL_5MW):
            NACA64: 0
            DU21:   1
            DU25:   2
            DU30:   3
            DU35:   4
            DU40:   5
            Round:  6
        :param alpha: angle of attack in degree
        :return: list Cl,Cd interpolate for this angle
        """
    ALPHA=Aero_Mat[:,1]
    THICKNESS=Aero_Mat[:,0]
    CL=Aero_Mat[:,2]
    CD=Aero_Mat[:,3]

    xi = THICKNESS
    yi = ALPHA
    zi1 = CL
    zi2 = CD

    f1 = interp2d(xi, yi, zi1, kind='linear', copy='True', bounds_error=False)
    f2 = interp2d(xi, yi, zi2, kind='linear', copy='True', bounds_error=False)

    """ Sans interpolation """
    #"""
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(ALPHA, THICKNESS)
    ax.plot_trisurf(ALPHA,THICKNESS, CD)
    plt.xlabel('Alpha attack angular in degree'), plt.ylabel('Thickness in percent'),  # plt.zlabel('Drag Coeffictient')
    plt.title('NACA64 Interp')
    plt.show()
    #"""



    """ f1 interp2d """
    """
    # Create cube to match aCube dimensions
    xi = THICKNESS
    yi = ALPHA
    zi = CL

    # Interpolate scattered points
    X, Y= numpy.meshgrid(xi, yi)
    bCube = griddata((xi, yi), zi, (X, Y), method='linear')
    print bCube
    print len(bCube)
    f=interp2d(xi,yi,zi)
    """
    """ Bon résultat with gridddata malgré un plateau à 0"""
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xi, yi, zi,c='r')
    ax.plot_surface(X, Y, bCube)
    # plt.xlabel('Alpha attack angular in degree'), plt.ylabel('Thickness in percent'),  # plt.zlabel('Drag Coeffictient')
    plt.title('NACA64 Interp')
    plt.show()
    """
    """ Plot interp2d function of thickness """
    """
    fig=plt.figure()
    plt.plot(xi,f1(xi,12))
    plt.show()
    """
    """ Plot interp2d function of angular """
    """
    fig = plt.figure()
    plt.plot(yi, f1(50, yi))
    plt.show()
    """


    """
    grid_x, grid_y=numpy.meshgrid(THICKNESS, ALPHA,sparse=True, copy=False)

    grid_z1=interp2d(grid_x,grid_y,CL,kind='quintic',bounds_error=True,fill_value=0)
    #grid_z1 = interp2d(grid_x, grid_y, CL, kind='linear')


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(ALPHA, THICKNESS)
    ax.plot_surface(grid_x, grid_y, grid_z1(THICKNESS,ALPHA))
    #plt.xlabel('Alpha attack angular in degree'), plt.ylabel('Thickness in percent'),  # plt.zlabel('Drag Coeffictient')
    plt.title('NACA64 Interp')
    plt.show()
    """

    """

    xArray=ALPHA
    yArray=THICKNESS
    heightArray=CL
    #point=
    #estimatedHeightList=

    numIndexes = 500
    xi = numpy.linspace(numpy.min(xArray), numpy.max(xArray), numIndexes)
    yi = numpy.linspace(numpy.min(yArray), numpy.max(yArray), numIndexes)

    XI, YI = numpy.meshgrid(ALPHA, THICKNESS)
    points = numpy.vstack((xArray, yArray)).T
    values = numpy.asarray(heightArray)
    points = numpy.asarray(points)
    #values = numpy.asarray(estimatedHeightList)
    DEM = griddata(xArray, yArray, heightArray, xi, yi, interp='linear')
    #DEM = interpolate.griddata(points, values, (XI, YI), method='linear')
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(ALPHA, THICKNESS)
    ax.plot_surface(XI, YI, DEM)
    plt.xlabel('Alpha attack angular in degree'), plt.ylabel('Thickness in percent'),  # plt.zlabel('Drag Coeffictient')
    plt.title('NACA64 Interp')
    plt.show()
    """

    return [f1,f2]


def interpolation2d_aerofoil_data_V2(Aero_Mat):
    """
        Function returning for a value of an angle the list (Cl,Cd) interpolate for this angle
        :param profile_index: Index of Profile Number (Exemple for NREL_5MW):
            NACA64: 0
            DU21:   1
            DU25:   2
            DU30:   3
            DU35:   4
            DU40:   5
            Round:  6
        :param alpha: angle of attack in degree
        :return: list Cl,Cd interpolate for this angle
        """
    ALPHA = Aero_Mat[:, 1]
    THICKNESS = Aero_Mat[:, 0]
    CL = Aero_Mat[:, 2]
    CD = Aero_Mat[:, 3]

    x = THICKNESS
    y = ALPHA
    z1 = CL
    z2 = CD

    #f1 = interp2d(xi, yi, zi1, kind='linear', copy='True', bounds_error=False)
    #f2 = interp2d(xi, yi, zi2, kind='linear', copy='True', bounds_error=False)

    f1=rbfi1 = Rbf(x,y,z1,function='linear')
    f2=rbfi2 = Rbf(x,y,z2,function='linear')
    #print rbfi1(50,0)
    #"""
    xi=numpy.linspace(0,100,180)
    yi=numpy.linspace(-180,180,180)

    di1 = rbfi1(xi,yi)
    di2 = rbfi2(xi,yi)
    """
    fig = plt.figure()
    ax = fig.add_subplot(121, projection='3d')
    ax.plot_trisurf(x, y, rbfi1(x, y))
    plt.title('Interpolation of airfoil profile lift coefficient'), plt.xlabel('Thickness in percent'), plt.ylabel(
        'Angular of attack in degree')
    ax = fig.add_subplot(122, projection='3d')
    ax.plot_trisurf(x, y, rbfi2(x, y))
    plt.title('Interpolation of airfoil profile drag coefficient'), plt.xlabel('Thickness in percent'), plt.ylabel(
        'Angular of attack in degree')
    plt.show()
    #"""
    """ Sans interpolation """
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(ALPHA, THICKNESS)
    ax.plot_trisurf(ALPHA, THICKNESS, CL)
    plt.xlabel('Alpha attack angular in degree'), plt.ylabel('Thickness in percent'),  # plt.zlabel('Drag Coeffictient')
    plt.title('NACA64 Interp')
    plt.show()
    #"""

    return [f1,f2]
######TEST#####

#######################################################################################################################
def Interpolate_ClCd_from_Data(WT):
    """
    Main Function wich extract Aerofoils data for a specific Blade from a file.pos
    and return the interpolated functions derived from this data


    :param WT: Name of the specific Blade
    :return: [f_Cl(thickness, alpha), f_Cd(thickness, alpha)]
    """
    """ Extract Data from a .pos file """
    DATA_Char=get_aerofoils_data(WT)
    Thickness=get_thickness_data(WT)
    """ Convert Char value in floatting Value """
    DATA_Float=CharToFloat(DATA_Char)

    """ Organize DATA to be more easily interpolated """
    DATA=Organize_data_V3(DATA_Float,Thickness)
    Aero_Mat=DATA

    """ Interpolate DATA to return a [[Cl(alpha),Cd()(alpha)],...for each profile]"""

    #L = interpolation2d_aerofoil_data(Aero_Mat)
    L = interpolation2d_aerofoil_data_V2(Aero_Mat)
    Cl = L[0]
    Cd = L[1]

    return [Cl,Cd]


########TEST#####
"""
WT='NREL_5MW'
F=Interpolate_ClCd_from_Data(WT)
print F

Cl=F[0]
Cd=F[1]
#"""
######TO PLOT WITH TEST####
"""
f=Cl
t=15
a=numpy.arange(-180,181,1)
plt.plot(a,f(t,a))
plt.show()

t=numpy.arange(0,10,0.5)
for i in t:
    plt.plot(a,f(i,a), label=str(i))
plt.legend()
plt.show()

plt.plot(a,f(35,a))
plt.show()
t=numpy.arange(40,105,5)
for i in t:
    plt.plot(a,f(i,a), label=str(i))
plt.legend()
plt.show()
#"""
