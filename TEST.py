"""
Created on Wed May 08 10:06:47 2018

@author: augus

"""
import numpy as np
import matplotlib.pyplot as plt
import math as m

#from Function_BEM_V3 import BEM
#from Function_BEM_V3_notwist_optimization import BEM
from Function_BEM_V4_notwist_optimization import BEM

from Read_Blade_Data import extract_hub_radius_and_total_radius


WT_char='NREL_5MW'
L_R=extract_hub_radius_and_total_radius(WT_char)
R_max = L_R[1]
#print R_max

#V_var=11.4                                      # Rated (reference NREL5MW)
V_var = 11.

number_blades=3.

#TSR_var=7.
#print TSR_var
#Rotor_Speed_rated=12.4                          # Rated (reference NREL5MW)
#Rotor_Speed_rated=12.126                        # sDWM Reference
Rotor_Speed_rated=10.609                         # put rpm from Flex5 BEM
TSR_var= Rotor_Speed_rated*m.pi/30*R_max/V_var   #

plot_bool=False
drag_bool=True
q_bool=True

N=100                # 100 is a good compromise to be accurate and not spend too much time
J=15                 # 15  is a good compromise to be accurate and not spend too much time

Pitch=0.             # Degree

L=BEM(V_var,WT_char,number_blades,TSR_var,plot_bool, drag_bool, q_bool,N,J,Pitch)
print L
print 'TSR rated is '+ str(TSR_var)
print 'P_aero_turbine is '+str(L[1])+' (kW).'
print 'Cp is '+str(L[2])
########################################################################################################################
#------------------------------POWER function of wind Speed for a fixed rotationnal Speed------------------------------#
"""
V_cut_in=10
V_cut_out=20
V_VAR=np.linspace(V_cut_in,V_cut_out,20)

Rotor_Speed_rated=12.126                        # sDWM Reference

plot_bool=False
drag_bool=True
q_bool=True
N=40
J=5

WindSpeed_Power = np.zeros((len(V_VAR),2))

for i in range(len(V_VAR)):
    V_var = V_VAR[i]
    TSR_var = Rotor_Speed_rated * m.pi / 30 * R_max / V_var  # Rated (reference NREL5MW)
    L = BEM(V_var,WT_char,number_blades,TSR_var,plot_bool, drag_bool, q_bool,N,J,Pitch)
    WindSpeed_Power[i,0] = V_var
    print L
    WindSpeed_Power[i,1] = L[1]
print 'WindSpeed_Power:'
print WindSpeed_Power
np.savetxt('C:/Users/augus/OneDrive/Documents/Stage/Codes/Comparison_BEM_sDWM/BEM_Power_function_of_WindSpeed.dat', WindSpeed_Power)

#"""
########################################################################################################################
#------------------------------FIND THE BEST TSR-----------------------------------------------------------------------#
"""
TSR_var=np.linspace(6.5,8.5,200)
P_wind=[]
P_out=[]
Cp=[]
Cp_max=0
tsr_best=0
for tsr in TSR_var:
    result=BEM(11,WT_char,3,tsr,False, True, True,200,J,Pitch)
    P_wind.append(result[0])
    P_out.append(result[1])
    Cp.append(result[2])
    if result[2]>Cp_max:
        Cp_max=result[2]
        tsr_best=tsr

plt.figure(1)
plt.title('Main Results function of TSR')
plt.subplot(121),plt.plot(TSR_var,P_wind,label='Wind Power'),plt.plot(TSR_var,P_out,label='WindTurbine Power'),plt.title('Power function of TSR'),plt.xlabel('Tip Speed Ratio'),plt.ylabel('Power (kW)')
plt.legend()
plt.subplot(122),plt.plot(TSR_var,Cp,label='Cp(TSR)'),plt.title('Power Coefficient function of TSR'),plt.xlabel('Tip Speed Ratio'),plt.ylabel('Cp')
plt.legend()
#plt.show()
print 'Result:'
print 'The Best Cp is: '+str(Cp_max)+' with TSR = '+str(tsr_best)
#"""
########################################################################################################################
#-----------------------------------PLOT CP FUNCTION OF U--------------------------------------------------------------#

#"""
V_inf = np.linspace(5.,20.,30)
Rotor_Speed_rated = 1.265    # rad/s                     # put rpm from Flex5 BEM


P_wind = []
P_out = []
Cp = []

for V in V_inf:
    tsr = R_max * Rotor_Speed_rated / V
    result = BEM(V, WT_char, 3, tsr, False, True, True, N,J,Pitch)
    P_wind.append(result[0])
    P_out.append(result[1])
    Cp.append(result[2])
plt.figure(2)
plt.plot(V_inf,Cp,label='Cp(V)')
plt.legend(), plt.title('CP function of Wind Speed')

#"""

########################################################################################################################
#------------------------------
"""
tsr=tsr_best
V_inf=np.linspace(1,25,50)
P_wind = []
P_out = []
Cp = []

for V in V_inf:
    result=BEM(V,WT_char,3,tsr,False,True,True,40,J,Pitch)
    P_wind.append(result[0])
    P_out.append(result[1])
    Cp.append(result[2])
plt.figure(3)
plt.subplot(121),plt.plot(V_inf,P_wind,label='Wind Power'),plt.plot(V_inf,P_out,label='P_out'),plt.xlabel('Wind Speed (m/s)'),plt.ylabel('Power (kW)'),plt.title('Power function of Wind Speed for TSR = '+str(tsr_best)),plt.legend()
plt.subplot(122),plt.plot(V_inf,Cp,label='Cp(WS)'),plt.xlabel('Wind Speed (m/s)'),plt.ylabel('Cp'),plt.title('Power Coeficient function of Wind Speed'),plt.legend()
#"""
########################################################################################################################
#----------------------------- NUMERICAL PARAMETERS ANALYSIS ----------------------------------------------------------#
#"""
N_list=range(5,501,10)
P_wind=[]
P_out=[]
Cp=[]

#tsr=tsr_best
tsr=7.39447236181 #NREL5MW (V3 Calculation)
#"""
"""
for n in N_list:
    result = BEM(11, WT_char, 3, tsr, False, True, True, n,J,Pitch)
    P_wind.append(result[0])
    P_out.append(result[1])
    Cp.append(result[2])
plt.figure(4),plt.plot(N_list,Cp,label='Cp(n)'),plt.xlabel('N parameter of blade discretization'),plt.ylabel('Cp')
plt.title('Numeric Parameters Analysis'),plt.legend()
#"""
"""
J_list=range(1,20)
P_wind=[]
P_out=[]
Cp=[]
for J in J_list:
    result = BEM(11, WT_char, 3, tsr, False, True, True, N, J, Pitch)
    P_wind.append(result[0])
    P_out.append(result[1])
    Cp.append(result[2])
plt.figure(5),plt.plot(J_list,Cp,label='Cp(J)'),plt.xlabel('J iterative parameters'),plt.ylabel('Cp')
plt.title('Numeric Parameters Analysis'),plt.legend()
#"""
########################################################################################################################
#------------------------------------------ PITCH ANALYSIS ------------------------------------------------------------#
"""
PITCH = np.linspace(-20,20,100)
P_wind=[]
P_out=[]
Cp=[]

#tsr=tsr_best
tsr=7.39447236181 #NREL5MW (V3 Calculation)

for Pitch in PITCH:
    result = BEM(11, WT_char, 3, tsr, False, True, True, N, J, Pitch)
    P_wind.append(result[0])
    P_out.append(result[1])
    Cp.append(result[2])

plt.figure(6),plt.plot(PITCH,Cp,label='Cp(Pitch)'),plt.xlabel('Pitch in Degree'),plt.ylabel('Cp')
plt.title('Pitch Analysis'),plt.legend()
#"""
#######################################################################################################################
plt.show()