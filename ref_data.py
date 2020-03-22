import math as m
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
print("Reference Data")
# stationary measurement series 1
#mass of the observers in kg
p1=95
p2=92
co_ord=74
o_1L=66
o_1R=61
o_2L=75
o_2R=78
o_3L=86
o_3R=68
fuel=4050*0.453592
aircraft=9165*0.453592
g=9.80665
L=-0.0065  #k/m
R=287.05
y=1.4
rho0=1.225  #kg/m3
s=30   #m2
b=15.911  #m
#finding the  pressure
P0=101325 #pa
T0=288.15 #K
A=b**2/s
W_s=60500   #from append b in Newton
chrd=2.0569
Lift=(p1+p2+co_ord+o_1L+o_1R+o_2L+o_2R+o_3L+o_3R+fuel+aircraft)*g
h_p=[5010,5020,5020,5030,5020,5110]                                                    #input heights in feet
P=[]
for i in h_p:
     hp = i * 0.3048
     z = (-g / (L * R))
     P=P+([P0*(1+(L*hp/T0))**z])

# finding mach number
#to make the equation easier its broken down into smaller equations
#v1=2/(y-1)
#v2=((y-1)*rho0)/(2*y*P0)
#v3=y/(y-1)
kinematic_factor=1.4207E-5
V_c=[249,221,192,163,130,118]    # input speeds
M=[]
Re=[]
for i in range(0,6):
    V_c[i]=V_c[i]*0.51444
    V_temp=V_c[i]**2   #speed square
    Re=Re+[(chrd*V_c[i])/kinematic_factor]
    c1=(y-1)/(2*y)
    c2=rho0/P0
    step_1=(1+c1*c2*V_temp)
    c3=y/(y-1)
    step_2=step_1**c3
    step_3=step_2-1
    c4=P0/P[i]
    step_4=c4*step_3
    step_5=step_4+1
    c5=(y-1)/y
    step_6=step_5**c5
    step_7=step_6-1
    c6=2/(y-1)
    step_8=c6*step_7
    M=M+[step_8**0.5]
    #v2=(v2*(V_c[i])**2+1)**v3-1
    #v4=((P0/P[i])*v2+1)**v3
    #v5=v1*(v4-1)
    #M=M+[v5**0.5]


#C_L=L/(0.5*rho*v**2*s)

#temperature correction
T_m=[12.5,10.5,8.8,7.2,6,5.2]
T=[]
for i in range(0,6):
    T_m[i]=T_m[i]+273.15
    T=T+[T_m[i]/(1+(y-1)/2*M[i]**2)]

#calculating speed of sound
a=[]
for i in range(0,6):
    a=a+[(y*R*T[i])**0.5]

#calculating true airspeed
vt=[]
for i in range(0,6):
    vt=vt+[M[i]*a[i]]
#calculating rho using perfect gas law
rho=[]
for i in range(0,6):
    rho=rho+[P[i]/(R*T[i])]

#calculating equivalent airspeed
v_eq=[]
for i in range(0,6):
    v_eq=v_eq+[vt[i]*((rho[i]/rho0)**0.5)]

# calculating Cl
f_u=[360,412,447,478,532,570]
c_l=[]
Lift_req=[]
for i in range(0,6):
    f_u[i]=f_u[i]*4.44822     # converting to nootuns
    Lift_1=Lift
    Lift_1=Lift-f_u[i]
    Lift_req=Lift_req+[Lift_1]
    c_l=c_l+[Lift_1/(0.5*rho[i]*v_eq[i]**2*s)]

 #calculating Cd
L_Thr=[3665.03,2995.38,2399.67,1863.38,1892.21,2208.82]
R_Thr=[3770.95,3057.27,2526.11,2015.87,2070.75,2405.27]
Tot_Thr=[]
c_d=[]
for i in range(0,6):
    Tot_Thr=Tot_Thr+[L_Thr[i]+R_Thr[i]]
    c_d=c_d+[Tot_Thr[i]/(0.5*rho[i]*v_eq[i]**2*s)]

c_lsqr=[]
for i in range(0,6):
    c_lsqr=c_lsqr+[c_l[i]**2]

e=(c_lsqr[-1]-c_lsqr[0])/(c_d[-1]-c_d[0])*(1/(m.pi*A))  #e
m=(1/(m.pi*A*e))
cd_0=c_d[0]-m*c_lsqr[0]     #cd_0


alph=[1.7,2.4,3.6,5.4,8.7,10.6]
plt.figure(2)
q=np.polyfit(c_d,c_l,2)
f=np.poly1d(q)
c_d_plot=np.array(c_d)
c_d_plot=np.sort(c_d_plot)    # for linspace
c_l_plot=np.array(c_l)
c_l_plot=np.sort(c_l_plot)
x_new = np.linspace(c_d_plot[0], c_d_plot[-1])
y_new = f(x_new)
plt.plot(c_d_plot,c_l_plot,'o',x_new,y_new)
plt.xlabel("C_D [-]")
plt.ylabel("C_L [-]")
plt.show()

plt.figure(3)
q=np.polyfit(alph,c_l,2)
f=np.poly1d(q)
c_d_plot=np.array(alph)
c_d_plot=np.sort(c_d_plot)    # for linspace
c_l_plot=np.array(c_l)
c_l_plot=np.sort(c_l_plot)
x_new = np.linspace(c_d_plot[0], c_d_plot[-1])
y_new = f(x_new)
plt.plot(c_d_plot,c_l_plot,'o',x_new,y_new)
plt.xlabel("alpha [°]")
plt.ylabel("C_L [-]")
plt.show()

c_l_alpha=(c_l[-1]-c_l[0])/(alph[-1]-alph[0])
print("Data for figure 2 and 3:")
print("Aircraft Configuration: Clean")
print("Reynolds number range:",min(Re),max(Re))
print("Mach Range:", min(M), max(M))
print()
#print(e,cd_0)


#  stationary measurement series 2
#  stationary measurement series 2
#  stationary measurement series 2
#  stationary measurement series 2
cm_tc=-0.0064
mfs=0.048  #kg/s
de=[0,-0.4,-0.9,-1.5,0.4,0.6,1,0,-0.5]
Tot_Thr2=[5494.68,5875.28,6289.03,6653.86,6411.24,6541.06,6517.09]
d_eng=0.686

h_p2=[6060,6350,6550,6880,6160,5810,5310,5730,5790]                                                    #input heights in feet
P2=[]
for i in h_p2:
     hp = i * 0.3048
     z = (-g / (L * R))
     P2=P2+([P0*(1+(L*hp/T0))**z])

V_c2=[161,150,140,130,173,179,192,161,161]    # input speeds
M2=[]
for i in range(0,9):
    V_c2[i]=V_c2[i]*0.51444
    V_temp=V_c2[i]**2   #speed square
    c1=(y-1)/(2*y)
    c2=rho0/P0
    step_1=(1+c1*c2*V_temp)
    c3=y/(y-1)
    step_2=step_1**c3
    step_3=step_2-1
    c4=P0/P2[i]
    step_4=c4*step_3
    step_5=step_4+1
    c5=(y-1)/y
    step_6=step_5**c5
    step_7=step_6-1
    c6=2/(y-1)
    step_8=c6*step_7
    M2=M2+[step_8**0.5]

T_m2=[5.5,4.5,3.5,2.5,5.0,6.2,8.2,5,5]
T2=[]
for i in range(0,9):
    T_m2[i]=T_m2[i]+273.15
    T2=T2+[T_m2[i]/(1+(y-1)/2*M2[i]**2)]

a2=[]
for i in range(0,9):
    a2=a2+[(y*R*T2[i])**0.5]

#calculating true airspeed
vt2=[]
for i in range(0,9):
    vt2=vt2+[M2[i]*a2[i]]
#calculating rho using perfect gas law
rho2=[]
for i in range(0,9):
    rho2=rho2+[P2[i]/(R*T2[i])]

#calculating equivalent airspeed
v_eq2=[]
for i in range(0,9):
    v_eq2=v_eq2+[vt2[i]*((rho2[i]/rho0)**0.5)]
# calculating Cl2
f_u2=[664,694,730,755,798,825,846,881,910]   #last 2 values for find Cn
c_l2=[]
Lift_req2=[]
for i in range(0,9):
    f_u2[i]=f_u2[i]*4.44822     # converting to nootuns
    Lift_1=Lift
    Lift_1=Lift-f_u2[i]
    Lift_req2=Lift_req2+[Lift_1]
    c_l2=c_l2+[Lift_1/(0.5*rho2[i]*v_eq2[i]**2*s)]
#print(Lift_req,Lift_req2)

v_tilda=[]
for i in range(0,7):
    v_tilda=v_tilda+[v_eq2[i]*(W_s/Lift_req2[i])**0.5]    #input for the graph

#print(v_tilda)
#print(v_eq2,a2,rho2,vt2)
Tot_Thr_s2=[2678.12,2801.72,2914.08,3040.54,2593.0,2508.5,2355.64]

TC2=[]
Tcs2=[]
for i in range(0,7):
    TC2=TC2+[2*Tot_Thr2[i]/(rho2[i]*(vt2[i]**2)*d_eng**2)]
    Tcs2 = Tcs2 + [2 * Tot_Thr_s2[i] / (rho2[i] * (vt2[i] ** 2) * d_eng ** 2)]

delta_xcg2=(134-288)*0.0254   #metres
cn2=c_l2[-1]
cm_delta2=-1/(de[-1])*cn2*(delta_xcg2/chrd)
#print(cm_delta2)

de_str2=[]
for i in range(0,7):
    de_str2=de_str2+[de[i]-1/cm_delta2*cm_tc*(Tcs2[i]-TC2[i])]

#Data from meas:
#Ve delta e:
v_tilda1_meas= [ 69.14524483 , 73.82340752 , 79.03515151  ,84.684758   , 90.76076007, 96.09362594, 101.00179236]
de_str2_meas= [-1.9909169 , -1.49242438, -0.99365751 ,-0.594742  , -0.19465784 , 0.20498867,0.50501256]
x_new1_meas= [ 69.14524483 , 69.79537846,  70.44551208 , 71.0956457  , 71.74577932,
  72.39591295 , 73.04604657  ,73.69618019 , 74.34631382 , 74.99644744,
  75.64658106 , 76.29671468  ,76.94684831  ,77.59698193 , 78.24711555,
  78.89724918 , 79.5473828   ,80.19751642,  80.84765005 , 81.49778367,
  82.14791729 , 82.79805091  ,83.44818454 , 84.09831816 , 84.74845178,
  85.39858541 , 86.04871903,  86.69885265  ,87.34898628 , 87.9991199,
  88.64925352 , 89.29938714 , 89.94952077,  90.59965439 , 91.24978801,
  91.89992164 , 92.55005526  ,93.20018888 , 93.85032251 , 94.50045613,
  95.15058975 , 95.80072337,  96.450857    ,97.10099062 , 97.75112424,
  98.40125787 , 99.05139149,  99.70152511 ,100.35165873 ,101.00179236]    
y_new1_meas= [-1.96074272, -1.89628794, -1.83243869 ,-1.76919497, -1.70655676, -1.64452409,
 -1.58309693, -1.52227531 ,-1.4620592  ,-1.40244862, -1.34344357, -1.28504404,
 -1.22725003, -1.17006155 ,-1.11347859 ,-1.05750116, -1.00212926, -0.94736287,
 -0.89320201, -0.83964668 ,-0.78669687 ,-0.73435258, -0.68261382, -0.63148059,
 -0.58095288, -0.53103069 ,-0.48171403 ,-0.43300289, -0.38489727, -0.33739718,
 -0.29050262, -0.24421358 ,-0.19853006 ,-0.15345207, -0.10897961, -0.06511266,
 -0.02185125,  0.02080465 , 0.06285502 , 0.10429986,  0.14513918,  0.18537298,
  0.22500125,  0.26402399 , 0.30244122 , 0.34025291,  0.37745909,  0.41405973,
  0.45005486,  0.48544446]    
#Ve Fe:
v_tilda_meas= [ 69.14524483 , 73.82340752,  79.03515151,  84.684758,    90.76076007,  96.09362594, 101.00179236]
Fe_air_eff_meas= [-39.10086568, -24.26634863, -12.64696884,   2.10147738 , 32.81117267,  64.66456428  ,90.39522629]
x_new_meas= [ 69.14524483 , 69.79537846 , 70.44551208,  71.0956457  , 71.74577932,
  72.39591295 , 73.04604657 , 73.69618019 , 74.34631382 , 74.99644744,
  75.64658106 , 76.29671468 , 76.94684831 , 77.59698193 , 78.24711555,
  78.89724918 , 79.5473828  , 80.19751642 , 80.84765005 , 81.49778367,
  82.14791729 , 82.79805091 , 83.44818454 , 84.09831816 , 84.74845178,
  85.39858541 , 86.04871903 , 86.69885265 , 87.34898628 , 87.9991199,
  88.64925352 , 89.29938714 , 89.94952077 , 90.59965439 , 91.24978801,
  91.89992164 , 92.55005526 , 93.20018888 , 93.85032251 , 94.50045613,
  95.15058975 , 95.80072337 , 96.450857   , 97.10099062 , 97.75112424,
  98.40125787 , 99.05139149 , 99.70152511 ,100.35165873 ,101.00179236]
y_new_meas= [-36.94666264 ,-35.89015747 ,-34.76833876 ,-33.58120653, -32.32876077,
 -31.01100148 ,-29.62792866, -28.17954231 ,-26.66584243, -25.08682902,
 -23.44250208 ,-21.73286161, -19.95790762 ,-18.11764009 ,-16.21205904,
 -14.24116446 ,-12.20495635, -10.10343471 , -7.93659954 , -5.70445084,
  -3.40698861 , -1.04421285,   1.38387644 ,  3.87727925 ,  6.4359956,
   9.06002547 , 11.74936887,  14.5040258  , 17.32399626 , 20.20928025,
  23.15987777 , 26.17578882,  29.2570134  , 32.4035515  , 35.61540314,
  38.8925683  , 42.235047  ,  45.64283922 , 49.11594497 , 52.65436425,
  56.25809706 , 59.9271434 ,  63.66150327 , 67.46117667 , 71.32616359,
  75.25646405 , 79.25207803,  83.31300555 , 87.43924659 , 91.63080116]

#print(de_str2)
#print(v_tilda)
plt.figure(0)
q1=np.polyfit(v_tilda,de_str2,2)
f1=np.poly1d(q1)
v_tilda1=np.array(v_tilda)
v_tilda1=np.sort(v_tilda1)    # for linspace sorted in ascending order, for this graph v_tilda is called v_tilda1
de_str2=np.array(de_str2)
de_str2=np.sort(de_str2)
#print(de_str2)
x_new1 = np.linspace(v_tilda1[0], v_tilda1[-1])
y_new1 = f1(x_new1)
plt.plot(v_tilda1,de_str2,'o', x_new1, y_new1,color="green",label="Reference Data")
plt.plot(v_tilda1_meas,de_str2_meas,'o',x_new1_meas, y_new1_meas,color="orange",label="Measurement Data")
plt.legend()
#plt.xlim([v_tilda[0]-1, v_tilda[-1] + 1 ])
plt.xlabel("ṽ_e [m/s]")
plt.ylabel("δ_e [°]")

#plt.scatter(v_tilda1,de_str2)
plt.gca().invert_yaxis()
plt.show()


Fe_ref2=[0,-23,-29,-46,26,40,83] #Fe ref2 in N
alpha2=[5.3,6.3,7.3,8.5,4.5,4.1,3.4] #Alpha ref2 in deg
de2=[0,-0.4,-0.9,-1.5,0.4,0.6,1]    #de ref2 in deg
slope_de_dalpha=(min(de2)-max(de2))/(max(alpha2)-min(alpha2))
#print(slope_de_dalpha)
cm_alpha2=-slope_de_dalpha*cm_delta2
#print(cm_alpha2)
Lift_req2=Lift_req2[0:-2]     # because the last 2 values where for finding Cn
Fe_air_eff=[]
for i in range (0,7):
    Fe_air_eff= Fe_air_eff+[Fe_ref2[i]*W_s/Lift_req2[i]]
    
    


plt.figure(1)
plt.gca().invert_yaxis()
q=np.polyfit(v_tilda,Fe_air_eff,2)
f=np.poly1d(q)
v_tilda=np.array(v_tilda)
v_tilda=np.sort(v_tilda)    # for linspace
Fe_air_eff=np.array(Fe_air_eff)
Fe_air_eff=np.sort(Fe_air_eff)
x_new = np.linspace(v_tilda[0], v_tilda[-1])
y_new = f(x_new)
plt.plot(v_tilda,Fe_air_eff,'o', x_new, y_new,color='green',label="Reference Data")
plt.plot(v_tilda_meas,Fe_air_eff_meas,'o',x_new_meas, y_new_meas,color="orange",label="Measurement Data")
plt.legend()
plt.xlim([v_tilda[0]-1, v_tilda[-1] + 1 ])
plt.xlabel("ṽ_e [m/s]")
plt.ylabel("F_e [N]")
plt.show()
#plt.scatter(alpha2,de2)
#plt.show()
#plt.scatter(v_tilda,de_str2)
#plt.gca().invert_yaxis()
#plt.show()

#Print required values:
print("Values required:")
print("Oswald factor:", e)
print("CD_0:",cd_0)
print("C_L_alpha=", c_l_alpha)
print("C_m_alpha=", cm_alpha2)
print("C_m_delta=", cm_delta2)
