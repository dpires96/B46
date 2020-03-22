import math as m
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
print("Measurement Data")
# stationary measurement series 1
#mass of the observers in kg
p1=102
p2=80
co_ord=183
o_1L=74
o_1R=80
o_2L=81
o_2R=91
o_3L=71
o_3R=99
fuel=2561*0.453592
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
h_p=[9000,9000,9000,9010,9030,9080]                                                    #input heights in feet
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
V_c=[250,225,187,161,133,114]    # input speeds
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
T_m=[8.5,6.1,3.2,2.8,1.5,0.5]
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
f_u=[433,472,530,553,582,613]
c_l=[]
Lift_req=[]
for i in range(0,6):
    f_u[i]=f_u[i]*4.44822     # converting to nootuns
    Lift_1=Lift
    Lift_1=Lift-f_u[i]
    Lift_req=Lift_req+[Lift_1]
    c_l=c_l+[Lift_1/(0.5*rho[i]*v_eq[i]**2*s)]

 #calculating Cd
L_Thr=[3440.3,2920.46,2106.24,1926.26,1510.33,1983.92]
R_Thr=[3780.36,3241.28,2458.8,2189.41,1923.65,2312.05]
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


alph=[1.4,2.0,3.4,4.9,7.7,10.6]

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

c_l_alpha=(max(c_l)-min(c_l))/(max(alph)-min(alph))
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
de=[-0.6,-1.0,-1.5,-2.0,-0.2,0.2,0.5,-0.5,-1.1]
Tot_Thr2=[5864.64,6135.92,6380.2,6662.9,6262.55,6311.75,6613.49]
d_eng=0.686
chrd=2.0569
h_p2=[8090,8430,8640,9000,8310,7990,7480,8470,8560]                                                    #input heights in feet
P2=[]
for i in h_p2:
     hp = i * 0.3048
     z = (-g / (L * R))
     P2=P2+([P0*(1+(L*hp/T0))**z])

V_c2=[161,150,140,131,172,182,191,162,162]    # input speeds
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

T_m2=[3.2,2.2,1.5,0.8,3.5,4.2,2.2,2.8,3.0]
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
f_u2=[680,719,733,754,774,794,835,887,918]   #last 2 values for find Cn
c_l2=[]
Lift_req2=[]
for i in range(0,9):
    f_u2[i]=f_u2[i]*4.44822     # converting to newtons
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
Tot_Thr_s2=[2848.48,2976.46,3091.32,3209.14,2768.16,2659.54,2572.3]

TC2=[]
Tcs2=[]
for i in range(0,7):
    TC2=TC2+[2*Tot_Thr2[i]/(rho2[i]*(vt2[i]**2)*d_eng**2)]
    Tcs2 = Tcs2 + [2 * Tot_Thr_s2[i] / (rho2[i] * (vt2[i] ** 2) * d_eng ** 2)]

delta_xcg2=(131-288)*0.0254   #metres
cn2=c_l2[-1]
cm_delta2=-1/(de[-1]-de[7])*cn2*(delta_xcg2/chrd)
#print(cm_delta2)

de_str2=[]
for i in range(0,7):
    de_str2=de_str2+[de[i]-1/cm_delta2*cm_tc*(Tcs2[i]-TC2[i])]

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
plt.plot(v_tilda1,de_str2,'o', x_new1, y_new1)
#plt.xlim([v_tilda[0]-1, v_tilda[-1] + 1 ])
plt.xlabel("ṽ_e [m/s]")
plt.ylabel("δ_e [°]")

#plt.scatter(v_tilda1,de_str2)
plt.gca().invert_yaxis()
plt.show()

#print("Stuff graphs ve delta e:")
#print("v_tilda1",v_tilda1)
#print("de_str2", de_str2)
#print("x_new1", x_new1)
#print("y_new1", y_new1)

Fe_ref2=[2,-12,-23,-37,31,61,85] #Fe ref2 in N
alpha2=[5.0,5.8,6.7,7.7,4.2,3.6,3.2] #Alpha ref2 in deg
de2=[-0.6,-1.0,-1.5,-2.0,-0.2,0.2,0.5]    #de ref2 in deg
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
plt.plot(v_tilda,Fe_air_eff,'o', x_new, y_new)
plt.xlim([v_tilda[0]-1, v_tilda[-1] + 1 ])
plt.xlabel("ṽ_e [m/s]")
plt.ylabel("F_e [N]")
plt.show()

#print("Stuff graphs ve Fe:")
#print("v_tilda",v_tilda1)
#print("Fe_air_eff", Fe_air_eff)
#print("x_new", x_new)
#print("y_new", y_new)
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