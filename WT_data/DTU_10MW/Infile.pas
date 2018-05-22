sweep
**** Blade data ****
DTU10.vda                    Blade data filename,  '-'= data below
**** Rotor + Nacelle ***
3                            B  (number of blades)
1.0 1.0  1.0 1.0  1.0 1.0    K_flap_fak  K_edge_fak, pairs, blade 1..3 ||Ændres hvis man vil ændre stivhed af vinger
1.0 1.0 1.0	             M_fak, blade 1..3  ||Ændres hvis man vil ændre stivhed af vinger
0  0  0                      Pitch offset blade 1..B (deg)
0                            X_rod (offset blade data x) (m) || Hvis første punkt på vinge har radius=0, selvom det er en distance (x_rod) fra centrum
-2.5   5                     coning  tilt  (deg, deg)
7.1  2.7                     Znav  Zrn    (rotor overhang + pos. of shaft bend) 
105520  0             	     Mnav  Zgnav
2.0E5  2.0E5  325671   	     Ixnav  Iynav  Iznav
1.3E10 1.3E10 2.05E9       ? Kgax  Kgay  Ktors  (shaft stiffness, DOF 14,15,28)  ||Scale 2.67 i forhold til NREL5MW
1500.5  50.0           	     Igenerator  Ngear
446034  2.69  2.75           Mkab  Zgkab  XKK2
7.33E6  7.33E6  4.1E6        Ixkab  Iykab  Izkab
3.6E9  1.0E10                Ktx  KKy  (yaw- and tilt stiffness) 
30  20   0.0                 CdARz CdARxy ZlatR  aero drag of nav ||Cd*A i z,x,y retning og distance fra punkt R til krafts angrebspunkt, kun vigtigt i storm
30  105  2.0         	     CdAKz CdAKxy ZlatK  aero drag of nacelle ||Cd*A i z,x,y retning og distance fra punkt K til krafts angrebspunkt, kun vigtigt i storm
.03 .03 .03 .03 .1     	     Damp. DOF 11,12 + 14,15 + shaft tors. (log.decr.)
****  tower data ****
DTU10MW_tower.tda            Tower data filename,  '-'= data below
****  foundation data ****
Fund_v02_h0.fda              Foundation data filename,  '-'= data below
**** operational data  ****
1.29  9.81           Ro  g
1.005  10  1  1         Omega  Tetap  Generator-on/off ((1 = on, 0 = off)
180  0                Psi   Yaw  (rotorpos  yawpos, deg,deg)
8.0  0.14  0  0      Vnav  Vexp  Vdir  Vslope
0.0  0.8  0.5        Turb.intens (u), Rel.ti. (v) and (w) 0.12  0.8  0.5
..\Wind\vs88152.int 0    Turbulence-filename  T-offset  (u)
..\Wind\vs88153v.int 0   Turbulence-filename  T-offset  (v)
..\Wind\vs88154w.int 0   Turbulence-filename  T-offset  (w)
**** data for simulation ***
1 1 1 0             	Blade-dof: 1F 2F 1K 2K,   1 = active, 0=stiff
1 1 1 0 0 1          	DOF 11..15 + 28 (shafttors) 1=active, 0=stiff, DOF11=yaw
1 1 1 1              	DOF 7..10 (twr: L1, L2, T1, T2)  L = long., T=transv.   
0 0 0 0 0 0           	DOF 1..6 (foundat: Tx, Tz,Ry, Ty,Rz, Rx) (T=transl,R=rot)
0 0.02 5550       	Tstart  dt  Tmax
0  20   2   2        	Printop  Nprint  Filop  Nfil
1 0.11 0.05 40 10 4  	Stall dClda dCldaS AlfS ALrund TauFak
1    1.0             	Dynamisk wake (0/1)  DynWTCfak
1.0  1.0             	Towershadow-factor  Cd-factor (twr)
***  data for generator, brake, pitch etc. ****
-Gendat.inf             Generatordata filename, '-'= data below
-Brakedat.inf           Brakedata filename,     '-'= data below
-Pitchdat.inf           Pitchdata filename,     '-'= data below
-Yawdat.inf             Yaw data filename,      '-'= data below
-Contrdat,inf           Controldata filename,   '-'= data belov
-Winddat.inf         	Winddata filename,      '-'= data below, '0'=no data
-Initdat.inf            Initdata filename,      '-'= default values used
-Restart.rstF          Restart filename,       '-'= default values used   -Restart.rstF
*** Inputdata for Flex5 generator version 2.1 ***
DTU10MWVS 
10000  0.05  0.2        Pref (kW), tau_g (s)  tau_set (s)
5.0  0.7             	F0 (Hz), Ksi (-), 2.order bandpass RPM-filter
1800                    D_gen (Nm/rpm)  Generator damping  
120  300  600           P_loss (kW) at 0%, 50%, 100% Pref
60 480                	P_loss_mech Nref (mech. loss (kW) at RPM, gen. off)
*** Inputdata for Flex5 brake model version 2.0 ***
DTU10MW 
52.3                	Dynamic brake moment on generator shaft (kNm)
62.8                 	Static brake moment (kNm)
0.7                   	tau      (s)
0.0                    	T-delay  (s)
0.0  0.25  0.05        	V (deg), R/V, K0/Ktors  ( mainshaft play )
*** Inputdata for Flex5 pitch model version 1.2 ***
DTU 10MW pitch servo  : "Blade pitch servo and generator models are not included in this controller and should be modeled separately, if they are to be included in the simulations."
6.3  0.8         	OMres (rad/s)  Ksi-rel (< 1) OMres=resonance frequency (rad/s)
*** Inputdata for Flex5 yaw model version 1.0 ***
DTU 10MW
1.0  0.2               	Yaw-rate (deg/s)  Yaw-tau (s)
*** Inputdata for Flex5 control system version 3.2 ***
DTU 10MW, variable speed
401.6  6800         	N1, Pow1  (rpm, kW) point on cubic part of P(rpm) curve
300 480                	N_min, N_max (rpm, rpm)
10000                	Pmax (kW)
2000 6000               KI KP (Nm/s/rpm, Nm/rpm) const RPM control
0.4 0.7                 F0 (Hz), Ksi (-)  2.order RPM-filter    0.2 0.7
0.045  0.20  0.0        KI KP KD (deg/s/rpm, deg/rpm, deg/rpm*s) 0.03 0.10 0
6.0          	        KK (deg) (reduction of gain as function of pitch) 6.0
90  10           	Teta_max Pitchrate_max  (deg,deg/s)
10                   	TauV (sec)  (wind-averaging for minimum pitch)
9          No. of lines in table for Wind, minPitch, pow(init) (m/s,deg, kW)
4    2.58   200
5    1.80   600
6    0.79  1500
7    0.0   2500
8    0.0   3600
9    0.0   5100
10   0.0   6800
12   0.0   9800
13   0.0  10000
-0.2  2  10          TePstart TePstop TePsstop (deg/s)
0.05                 Tsamp  (sec)   control system
*** Inputdata for deterministic wind, version 1.1 ***
NM80 test 
94                   N  number of lines in interpolation table below
40      3        0           T  Vnav  Vdir    (s, m/s, deg)
          41           4           0
         160           4           0
         161           5           0
         280           5           0
         281           6           0
         400           6           0
         401           7           0
         520           7           0
         521           8           0
         640           8           0
         641           9           0
         760           9           0
         761          10           0
         880          10           0
         881          11           0
        1000          11           0
        1001          12           0
        1120          12           0
        1121          13           0
        1240          13           0
        1241          14           0
        1360          14           0
        1361          15           0
        1480          15           0
        1481          16           0
        1600          16           0
        1601          17           0
        1720          17           0
        1721          18           0
        1840          18           0
        1841          19           0
        1960          19           0
        1961          20           0
        2080          20           0
        2081          21           0
        2200          21           0
        2201          22           0
        2320          22           0
        2321          23           0
        2440          23           0
        2441          24           0
        2560          24           0
        2561          25           0
        2680          25           0
        2681          25           0
        2800          25           0
        2801          25           0
        2920          25           0
        2921          24           0
        3040          24           0
        3041          23           0
        3160          23           0
        3161          22           0
        3280          22           0
        3281          21           0
        3400          21           0
        3401          20           0
        3520          20           0
        3521          19           0
        3640          19           0
        3641          18           0
        3760          18           0
        3761          17           0
        3880          17           0
        3881          16           0
        4000          16           0
        4001          15           0
        4120          15           0
        4121          14           0
        4240          14           0
        4241          13           0
        4360          13           0
        4361          12           0
        4480          12           0
        4481          11           0
        4600          11           0
        4601          10           0
        4720          10           0
        4721           9           0
        4840           9           0
        4841           8           0
        4960           8           0
        4961           7           0
        5080           7           0
        5081           6           0
        5200           6           0
        5201           5           0
        5320           5           0
        5321           4           0
        5440           4           0
        5441           3           0
        5560           3           0
        5561           3           0

