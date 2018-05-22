# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 15:40:47 2018

@author: augus
Reference: Grant Ingram, Wind Turbine Blade Analysis using the Blade
Element Momentum Method.
"""
########################################################################################################################
import numpy as np
import math as m
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from Read_ClCd_Data_V2 import Interpolate_ClCd_from_Data
from Read_Thickness_Data import get_thickness_data
from Read_Blade_Data import Interpolate_Blade_Thickness_Chord_function_of_radius, extract_hub_radius_and_total_radius


def BEM(V_var,WT_char,number_blades,TSR_var,plot_bool, drag_bool, q_bool,N_discre,J,PITCH):
    """

    :param V_var: Wind Speed (m/s)
    :param WT_char: char defining the kind of wind turbine
    :param number_blades: number of blade for the wind turbine
    :param TSR_var: Global Tip Speed Ratio
    :param plot_bool: boolean to know if I want to plot some result along the computation
    :param drag_bool: boolean to know if I want drag effect or not
    :param q_bool: boolean to know if I want Glauert correction or not
    :param N_discre: discretization of the blade
    :param J: number of iterative calculation on the blade
    :param PITCH: degree, Global Pitch Angle for the blade
    :return:
    """
    ########################################################################################################################
    """ Flow parameters """
    V_inf=V_var                                 # velocity of the free stream (m/s)
    rho_inf=1.225                               # flow density  (kg/m^3)

    ########################################################################################################################
    """ Wind Turbine Parameters """
    #WT='NREL_5MW'
    #WT='DTU_10MW'
    WT=WT_char

    B=number_blades                             # Number of Blade
    L_R=extract_hub_radius_and_total_radius(WT)
    R=L_R[1]                                    # Radius of the wind turbine rotor (total radius root included)
    R0=L_R[0]                                   # Root radius
    R_end_cylinder = L_R[2]
    TSR=TSR_var
    lambd=TSR
    Omega=TSR*V_inf/R                           # Rotationnal speed of the rotor calculated with TSR
    lambd_root=Omega*R0/V_inf
    lambd_end_cylinder=Omega*R_end_cylinder/V_inf

    ########################################################################################################################
    print "Rotationnal Speed is "+str(Omega*30/m.pi)+" tr/mn."
    ########################################################################################################################
    """Numerical Parameters"""
    N=N_discre                                              # Dividing the blade into N element. Typically 10 to 20
    #r = np.linspace(R0,R,N)                                # Non-Normalized radius discretization (put in commentary the  function which normalized radius in Read_real_blade
    #r = np.linspace(0.09*R,R,N)
    r = np.linspace(R_end_cylinder,R,N)
    r = r[:]
    #lambd_r=np.linspace(lambd_root,lambd,N)                # Lambda discretization
    lambd_r = np.linspace(lambd_end_cylinder, lambd, N )    # Lambda discretization
    lambd_r = lambd_r[ :]
    ########################################################################################################################
    """ AeroDynamic Factor (VERSION 2): Interpolate Aerofoil data """
    F=Interpolate_ClCd_from_Data(WT)             # Extract Cl(Thickness, Alpha), Cd(Thickness, Alpha) for each airfoil
    f_Cl=F[0]
    f_Cd=F[1]
    ########################################################################################################################
    """ Blade Chord and Thickness Interpolation function of normalized radius """
    F = Interpolate_Blade_Thickness_Chord_function_of_radius(WT)
    f_thick = F[0]
    f_chord = F[1]
    f_gamma = F[2]

    #######PLOT#####
    #"""

    if plot_bool:
        x=r
        plt.subplot(221)
        plt.plot(x,f_thick(x),label='Thickness(r)'), plt.xlabel('r'),plt.ylabel('percent'), plt.legend()#, plt.show()
        plt.subplot(212)
        plt.plot(x,f_chord(x),label='Chord(r)'), plt.xlabel('r'),plt.ylabel('m'), plt.legend()#, plt.show()
        plt.subplot(222)
        plt.plot(x,f_gamma(x) ,label='Turbine Twist (Global pitch = '+str(PITCH)+')'),plt.plot(x,f_gamma(x)+PITCH)
        plt.xlabel('r'),plt.ylabel('degree'), plt.legend()#, plt.show()
        plt.show()
    #"""
    ########################################################################################################################
    """ Shape Parameters Initialization"""

    beta0=[]                                                                           # Relative flow angle onto blades
    c0=f_chord(r)                                                                      # Chord distribution
    sigma0=B*c0/(2*m.pi*r)                                                             # Local Solidity
    i0=[]                                                                              # Incidence Angle
    gamma0=m.pi/2-(f_gamma(r)+PITCH)*m.pi/180.                                         # Aerofoil inlet angle
    Cl0=[]                                                                             # Lift Coefficient
    Cd0=[]                                                                             # Drag Coefficient
    a1=[]                                                                              # Axial Induction Factor
    a2=[]                                                                              # Angular Induction Factor

    ########################################################################################################################
    """ Tip loss Correction """
    Q0=[]

    ########################################################################################################################
    """ Initial shape and pre-iteration """

    for k in range(N):

        """ Geometrical Computation """
        beta0 = beta0 + [m.pi / 2. - 2. / 3. * np.arctan(1. / lambd_r[k])]
        i0.append( gamma0[k] - beta0[k])

        """ Tip Loss Correction """
        if q_bool:
            Q0 = Q0 + [2. * np.arccos(np.exp(-(B / 2. * (R - r[k])) / (r[k] * np.cos(beta0[k])))) / m.pi]
        else:
            Q0 = Q0 +[1.]

        """ Get Cd, Cl """
        Cl0=Cl0+[f_Cl(f_thick(r[k]),i0[k]*180./m.pi)]
        if drag_bool:
            Cd0=Cd0+[f_Cd(f_thick(r[k]),i0[k]*180./m.pi)]
        else:
            Cd0= Cd0 + [0.]

        """ Initial Guess of induction factor """

        a1 = a1 + [ (1. + (4. * Q0[k] * (np.cos(beta0[k])) ** 2) / (sigma0[k] * (Cl0[k] * np.sin(beta0[k]) + Cd0[k] * np.cos(beta0[k])))) ** (-1.)]
        a2 = a2 + [sigma0[k] * ((Cl0[k] * np.cos(beta0[k]) - Cd0[k] * np.sin(beta0[k])) / (4 * Q0[k] * lambd_r[k] * (np.cos(beta0[k])) ** 2)) * (1. - a1[k])]

        #a1=a1+[(1.+(4.*1.*(np.cos(beta0[k]))**2.)/(sigma0[k]*(Cl0[k]*np.sin(beta0[k]))))**(-1.)]
        #a2=a2+[(1.-3.*a1[k])/(4.*a1[k]-1.)]
    beta=list(beta0)
    gamma=list(gamma0)
    i=list(i0)

    Cl=list(Cl0)
    Cd=list(Cd0)

    c=list(c0)
    Q=list(Q0)
    """ Angular in degree"""
    #"""
    Phi=[]
    Gamma=[]
    Alpha=[]

    for k in range(N):
        Phi=Phi+[90-beta[k]*180./m.pi]
        Alpha=Alpha+[i[k]*180./m.pi]
        Gamma=Gamma+[90-gamma[k]*180./m.pi]
    #"""
    ########PLOT#######
    #"""
    if plot_bool:
        plt.plot(r,sigma0,label='sigma'), plt.legend(),plt.show()
        plt.subplot(221)
        plt.plot(r, Phi, label='Phi'), plt.plot(r, Gamma, label='Gamma'), plt.plot(r, Alpha, label='Alpha'), plt.legend()
        plt.subplot(222)
        plt.plot(r,Q,label='Tip Loss Correction'), plt.legend()
        plt.subplot(223)
        plt.plot(r,Cl,label='Cl'), plt.plot(r,Cd,label='Cd'), plt.legend()
        plt.subplot(224)
        plt.plot(r,a1,label='a1'), plt.plot(r,a2,label='a2'), plt.legend()
        plt.show()
    #"""
    ########################################################################################################################
    """ iterations """
    for j in range(J):

        for k in range(N):


            """ Step 1 """
            beta[k] = np.arctan(lambd_r[k] * (1. + a2[k]) / (1. - a1[k]))

            """ Step 2 """
            i[k]=gamma[k]-beta[k]

            """ Tip Loss Correction """
            if q_bool:
                Q[k] = 2. * np.arccos(np.exp(-(B / 2. * (R - r[k])) / (r[k] * np.cos(beta[k])))) / m.pi
            """ Get Cd, Cl """
            Cl[k]=f_Cl(f_thick(r[k]),i[k]*180./m.pi)
            if drag_bool:
                Cd[k]=f_Cd(f_thick(r[k]),i[k]*180./m.pi)
            """ Step 3 """
            a1[k]=(1. + (4.*Q[k]*(np.cos(beta[k]))**2) / (sigma0[k]*(Cl[k]*np.sin(beta[k])+Cd[k]*np.cos(beta[k])) ) )**(-1.)
            a2[k]=sigma0[k]*((Cl[k]*np.cos(beta[k])-Cd[k]*np.sin(beta[k]))/(4*Q[k]*lambd_r[k]*(np.cos(beta[k]))**2))*(1.-a1[k])

    """ Angular in degree"""
    #"""
    # Angle in Søren Reference
    Phi=[]
    Gamma=[]
    Alpha=[]

    for k in range(N):
        Phi=Phi+[90-beta[k]*180./m.pi]
        Alpha=Alpha+[i[k]*180./m.pi]
        Gamma=Gamma+[90-gamma[k]*180./m.pi]

    #"""
    ########PLOT########
    #"""
    if plot_bool:
        plt.subplot(221)
        plt.title('after convergence'), plt.plot(r,Phi,label='Phi'), plt.plot(r,Gamma,label='Gamma'), plt.plot(r,Alpha,label='Alpha'),plt.ylabel('Angle in degree'), plt.legend()

        plt.subplot(222)
        plt.title('after convergence'), plt.plot(r,Cl,label='Cl'), plt.plot(r,Cd,label='Cd'), plt.legend()

        plt.subplot(223)
        plt.title('after convergence'), plt.plot(r,a1,label='a1'), plt.plot(r,a2,label='a2'), plt.legend()

        plt.subplot(224)
        plt.plot(r, Q, label='Tip Loss Correction'), plt.legend()

        plt.show()
    #"""
    ########################################################################################################################
    """ Power Computation """
    Cp=0.


    Cp_plot=[]

    ##### no drag, no tip losses
    def f_x_noDrag_noQ(k, lambd_r, a1, a2):
        return lambd_r[k] ** (3) * a2[k] * (1 - a1[k])
    ###### the good f_x
    def f(k,lambd_r,a1,a2,beta,Cl,Cd,Q):
        return Q[k] * lambd_r[k] ** (3) * a2[k] * (1 - a1[k]) * (1 - np.tan(beta[k])*Cd[k]/Cl[k])

    Cp = 0.
    Cp_plot = []
    #f_x=f_x_noDrag_noQ
    f_x=f
    for k in range(N-2):

        ###### TRAPEZIUM RULES ########
        y_k_0 = f_x(k, lambd_r, a1, a2,beta,Cl,Cd,Q)
        y_k_1 = f_x(k + 1, lambd_r, a1, a2,beta,Cl,Cd,Q)
        Cp_k = (lambd_r[k + 1] - lambd_r[k]) * (y_k_0 + y_k_1) / 2
        Cp = Cp + Cp_k
        Cp_plot=Cp_plot+[Cp_k]

    if plot_bool:

        plt.plot(r[:N-2],Cp_plot,label='Cp(r)')
        plt.title('Cp evolution along iterative computation')
        plt.show()
    Cp=8./lambd**(2)*Cp
    P_wind=0.5*m.pi*rho_inf*R**2*V_inf**3
    P_Cp=Cp*P_wind
    ########################################################################################################################
    """Printing/Plotting Part"""

    #print "Results:"
    #print "Wind power is "+str(P_wind/1000.)+"kW"
    #print "Power output is "+str(P_abs/1000)+"kW"
    #print "Power output is "+ str(P_abs_T/1000)+"kW"
    #print "P_Cp is "+str(P_Cp/1000.)+"kW"
    #print "Cp is "+ str(Cp)

    #plt.plot(rm,cm)
    #plt.plot(rm,gammad)
    #plt.plot(rm[:],ide[:])
    #plt.plot(rm,betad)
    #plt.show()
    return [P_wind/1000,P_Cp/1000,Cp]
