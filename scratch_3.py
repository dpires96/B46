import math as m
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# stationary measurement series 1
#mass of the observers in kg
p1=102
p2=80
o_1L=74
o_1R=80
o_2L=81
o_2R=91
o_3L=71
o_3R=99
co_ord=183
M_pass=[p1,p2,o_1L,o_1R,o_2L,o_2R,o_3L,o_3R,co_ord]
fuel=2561*0.453592
aircraft=9165*0.453592
aircraft_pounds=aircraft/0.453592
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
Mass_Pass_tot_pound=(p1+p2+co_ord+o_1L+o_1R+o_2L+o_2R+o_3L+o_3R)/0.453592
Lift=(p1+p2+co_ord+o_1L+o_1R+o_2L+o_2R+o_3L+o_3R+fuel+aircraft)*g
h_p=[9000,9000,9000,9010,9030,9080]  #input heights in feet
V_c=[250,225,187,161,133,114]    # input speeds
T_m=[8.5,6.1,3.2,2.8,1.5,0.5]
f_u=[433,472,530,553,582,613]
L_Thr=[3440.3,2920.46,2106.24,1926.26,1510.33,1983.92]
R_Thr=[3780.36,3241.28,2458.8,2189.41,1923.65,2312.05]
alph=[1.4,2.0,3.4,4.9,7.7,10.6]
de=[-0.6,-1.0,-1.5,-2.0,-0.2,0.2,0.5,-0.5,-1.1]
Tot_Thr2=[5864.64,6135.92,6380.2,6662.9,6262.55,6311.75,6613.49]
h_p2=[8090,8430,8640,9000,8310,7990,7480,8470,8560]
V_c2=[161,150,140,131,172,182,191,162,162]    # input speeds
T_m2=[3.2,2.2,1.5,0.8,3.5,4.2,2.2,2.8,3.0]
f_u2=[680,719,733,754,774,794,835,887,918]   #last 2 values for find Cn
Tot_Thr_s2=[2848.48,2976.46,3091.32,3209.14,2768.16,2659.54,2572.3] #there is an input that needs to be changed in line 242
Fe_ref2=[2,-12,-23,-37,31,61,85] #Fe ref2 in N
alpha2=[5.0,5.8,6.7,7.7,4.2,3.6,3.2] #Alpha ref2 in deg
de2=[-0.6,-1.0,-1.5,-2.0,-0.2,0.2,0.5]    #de ref2 in deg
LEMAC=261.56
MAC=80.98
Time=[0,1157,1297,1426,1564,1787,1920,2239,2351,2484,2576,2741,2840,2920,3062,3166]
# input at line 250 needs to be updated if the delta xcg changes. This delta cg is computed at the end of the code
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
M=[]
for i in range(0,6):
    V_c[i]=V_c[i]*0.51444
    V_temp=V_c[i]**2   #speed square
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
c_l=[]
Lift_req=[]
for i in range(0,6):
    f_u[i]=f_u[i]*4.44822     # converting to nootuns
    Lift_1=Lift
    Lift_1=Lift-f_u[i]
    Lift_req=Lift_req+[Lift_1]
    c_l=c_l+[Lift_1/(0.5*rho0*v_eq[i]**2*s)]

 #calculating Cd
Tot_Thr=[]
c_d=[]
for i in range(0,6):
    Tot_Thr=Tot_Thr+[L_Thr[i]+R_Thr[i]]
    c_d=c_d+[Tot_Thr[i]/(0.5*rho0*v_eq[i]**2*s)]

c_lsqr=[]
for i in range(0,6):
    c_lsqr=c_lsqr+[c_l[i]**2]

e=(c_lsqr[-1]-c_lsqr[0])/(c_d[-1]-c_d[0])*(1/(m.pi*A))  #e
m=(1/(m.pi*A*e))
cd_0=c_d[0]-m*c_lsqr[0]     #cd_0
print("e=",e)
print("cdo",cd_0)
alph_rad=[]
for i in range(0,6):
    alph_rad=alph_rad+[alph[i]*3.14/180]
#print(alph_rad)

cl_alpha=(c_l[-1]-c_l[0])/(alph_rad[-1]-alph_rad[0])
#g=np.polyfit(alph_rad,c_l,1)
#g1=np.poly1d(g)
print("cl_alpha",cl_alpha)

print("cd-cl stuff")
q=np.polyfit(c_d,c_l,2)
f=np.poly1d(q)
c_d_plot=np.array(c_d)
c_d_plot=np.sort(c_d_plot)    # for linspace
c_l_plot=np.array(c_l)
c_l_plot=np.sort(c_l_plot)
x_new = np.linspace(c_d_plot[0], c_d_plot[-1])
y_new = f(x_new)
print("c_d_plot", c_d_plot)
print("c_l_plot",c_l_plot)
print("x_new",x_new)
print("y_new",y_new)
print("end of cd-cl stuff")

q=np.polyfit(alph_rad,c_l,1)
c_l_alpha=q[0]
f=np.poly1d(q)
c_d_plot=np.array(alph_rad)
c_d_plot=np.sort(c_d_plot)    # for linspace
c_l_plot=np.array(c_l)
c_l_plot=np.sort(c_l_plot)
x_new = np.linspace(c_d_plot[0], c_d_plot[-1])
y_new = f(x_new)
print("stuff for cl-alpha")
print("alpha_plot", c_d_plot)
print("c_l_plot",c_l_plot)
print("x_new",x_new)
print("y_new",y_new)
print("end of alpha-cl stuff")
#v_tilda=np.array(v_tilda)
#v_tilda=np.sort(v_tilda)    # for linspace
#Fe_air_eff=np.array(Fe_air_eff)
#Fe_air_eff=np.sort(Fe_air_eff)
#x_new = np.linspace(v_tilda[0], v_tilda[-1])
#y_new = f(x_new)
#plt.plot(v_tilda,Fe_air_eff,'o', x_new, y_new)
#plt.xlim([v_tilda[0]-1, v_tilda[-1] + 1 ])
#print(cl_alpha)
#plt.plot(c_lsqr,c_d)
#plt.scatter(alph_rad,c_l)
#plt.plot(alph,c_d)
#plt.show()

#  stationary measurement series 2
#  stationary measurement series 2
#  stationary measurement series 2
#  stationary measurement series 2
cm_tc=-0.0064
mfs=0.048  #kg/s
d_eng=0.686
chrd=2.0569
P2=[]
for i in h_p2:
     hp = i * 0.3048
     z = (-g / (L * R))
     P2=P2+([P0*(1+(L*hp/T0))**z])

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
c_l2=[]
Lift_req2=[]
for i in range(0,9):
    f_u2[i]=f_u2[i]*4.44822     # converting to nootuns
    Lift_1=Lift
    Lift_1=Lift-f_u2[i]
    Lift_req2=Lift_req2+[Lift_1]
    c_l2=c_l2+[Lift_1/(0.5*rho0*v_eq2[i]**2*s)]
#print(Lift_req,Lift_req2)

v_tilda=[]
for i in range(0,7):
    v_tilda=v_tilda+[v_eq2[i]*(W_s/Lift_req2[i])**0.5]    #input for the graph

#print(v_tilda)
#print(v_eq2,a2,rho2,vt2)

TC2=[]
Tcs2=[]
for i in range(0,7):
    TC2=TC2+[2*Tot_Thr2[i]/(rho2[i]*(vt2[i]**2)*d_eng**2)]
    Tcs2 = Tcs2 + [2 * Tot_Thr_s2[i] / (rho2[i] * (vt2[i] ** 2) * d_eng ** 2)]
deg_rad=3.14/180
de1=(de[-1]-de[-2])*deg_rad
delta_xcg2=(-2.696)*0.0254   #metres 2,696
cn2=c_l2[-1]
cm_delta2=-1/(de1)*cn2*(delta_xcg2/chrd)
print("cm_delta",cm_delta2)

de_str2=[]
for i in range(0,7):
    de_str2=de_str2+[de[i]-1/cm_delta2*cm_tc*(Tcs2[i]-TC2[i])]

#print(de_str2)
#print(v_tilda)
#q1=np.polyfit(v_tilda,de_str2,2)
#f1=np.poly1d(q1)
#v_tilda1=np.array(v_tilda)
#v_tilda1=np.sort(v_tilda1)    # for linspace sorted in ascending order, for this graph v_tilda is called v_tilda1
#de_str2=np.array(de_str2)
#de_str2=np.sort(de_str2)
##print(de_str2)
#x_new1 = np.linspace(v_tilda1[0], v_tilda1[-1])
#y_new1 = f1(x_new1)
#print("delta_e-ve stuff")
#print("v_tilda", v_tilda1)
#print("delta_e", de_str2)
#print("x_new", x_new1)
#print("y_new", y_new1)
#print("end of delta_e-ve-sutff")
#plt.plot(v_tilda1,de_str2,'o', x_new1, y_new1)
##plt.xlim([v_tilda[0]-1, v_tilda[-1] + 1 ])
#plt.xlabel("ṽe")
#plt.ylabel("δ*e")

#plt.scatter(v_tilda1,de_str2)
plt.gca().invert_yaxis()
plt.show()



slope_de_dalpha=(min(de2)-max(de2))/(max(alpha2)-min(alpha2))
#print(slope_de_dalpha)
cm_alpha2=-slope_de_dalpha*cm_delta2
Lift_req2=Lift_req2[0:-2]     # because the last 2 values where for finding Cn
Fe_air_eff=[]
for i in range (0,7):
    Fe_air_eff= Fe_air_eff+[Fe_ref2[i]*W_s/Lift_req2[i]]
print("cm_alpha",cm_alpha2)
#plt.scatter(alpha2,de2)
#plt.show()
#plt.scatter(v_tilda,Fe_air_eff)
#plt.gca().invert_yaxis()
#q=np.polyfit(v_tilda,Fe_air_eff,2)
#f=np.poly1d(q)
#v_tilda=np.array(v_tilda)
#v_tilda=np.sort(v_tilda)    # for linspace
#Fe_air_eff=np.array(Fe_air_eff)
#Fe_air_eff=np.sort(Fe_air_eff)
#x_new = np.linspace(v_tilda[0], v_tilda[-1])
#y_new = f(x_new)
#print("fe-ve stuff")
#print("v_tilda", v_tilda)
#print("fe", Fe_air_eff)
#print("x_new", x_new)
#print("y_new", y_new)
#print("end of fe-ve-sutff")
#plt.plot(v_tilda,Fe_air_eff,'o', x_new, y_new)
#plt.xlim([v_tilda[0]-1, v_tilda[-1] + 1 ])
#plt.xlabel("ṽe")
#plt.ylabel("F*e")
#plt.show()

##### finding mass balance #####
inch_metre=0.0254
pound_to_kg=0.453592
kg_to_pound=2.20462
#xcg of aircraft BEM
#using data from Mass and Balance report on BS
#2672953.5
#67874.26553958168
xcg_BEM=291.65 #inches
#xcg_BEM_tot=2672953.5
#xcg_BEM_tot_metre=xcg_BEM_tot*inch_metre*pound_to_kg
#xcg of pilots,co-ordinater and students in inches
xcg_p1=131
xcg_p2=131
xcg_o1L=214
xcg_o1R=214
xcg_o2L=251
xcg_o2R=251
xcg_o3L=288
xcg_o3R=288
xcg_cord=170
xcg_pass=[xcg_p1,xcg_p2,xcg_o1L,xcg_o1R,xcg_o2L,xcg_o2R,xcg_o3L,xcg_o3R,xcg_cord]
xcg_pass_tot=0
M1=0
#xcg for both cases 1)normal 2)when o3 moves to 131 (cockpit)
Mass_pass_pound=[]
for i in range(0,9):
    Mass_pass_pound =Mass_pass_pound+[M_pass[i]*kg_to_pound]
    xcg_pass_tot=(M1*xcg_pass_tot+Mass_pass_pound[i]*xcg_pass[i])/(M1+Mass_pass_pound[i])  #inches for case 1
    M1=M1+Mass_pass_pound[i]
    #print(xcg_pass_tot)
mzfw=aircraft_pounds+Mass_Pass_tot_pound
M_123=M1-(o_3L*kg_to_pound)
xcg_pass_tot2=((xcg_pass_tot*M1)-(288*o_3R*kg_to_pound))/(M_123)  #change this line
xcg_pass_tot2=((xcg_pass_tot2*M_123)+(131*o_3R*kg_to_pound))/(M1)  #case 2 xcg
xcg_zfw=((M1*xcg_pass_tot)+(xcg_BEM*aircraft_pounds))/(aircraft_pounds+M1)
xcg_zfw2=((M1*xcg_pass_tot2)+(xcg_BEM*aircraft_pounds))/(aircraft_pounds+M1)
print(xcg_zfw2,xcg_zfw)
#print(xcg_zfw)
#print(xcg_pass_tot)
#xcg of aircraft fuel
f_u_totl=f_u+f_u2
fuel_left=[fuel]
for i in range(len(f_u_totl)):
    f_u_totl[i]=f_u_totl[i]/g    #converting back from newtons to kg
    fuel_left=fuel_left+[fuel-f_u_totl[i]] #in kg

fuel_left_pounds=[]             #for mass balance form
for i in range(len(fuel_left)):
    fuel_left_pounds=fuel_left_pounds+[fuel_left[i]*kg_to_pound]

moment_pounds=[298.16,591.18,879.08,1165.42,1448.40,1732.53,2014.80,2298.84,2581.92,2866.30,3150.18,3434.52,3718.52,4003.23,4287.76,4572.24,4856.56,5141.16,5425.64,5709.90,5994.04,6278.47,6562.82,6846.96,7131,7415.33,7699.60,7984.34,8269.06,8554.05,8839.04,9124.80,9410.62,9696.97,9983.40,10270.08,10556.80,10843.87,11131,11418.2,11705.5,11993.31,12281.18,12569.04,12856.86,13144.73,13432.48,13720.56,14008.46,14320.34]
mass_pounds=[]
for i in range(1,51):
    mass_pounds=mass_pounds+[i*100]
q_1=np.polyfit(mass_pounds,moment_pounds,1)  #for finding equation of mass_pounds to moment_pounds
#equation for exact moment_pounds
moment_pounds_final=[]
for i in range(len(fuel_left_pounds)):
    moment_pounds_final=moment_pounds_final+[(q_1[0]*fuel_left_pounds[i]+q_1[1])*100]
xcg_fuel_tot_inch=[]
for i in range(len(moment_pounds_final)):
        xcg_fuel_tot_inch=xcg_fuel_tot_inch+[(moment_pounds_final[i]/fuel_left_pounds[i])*100]
#finding the total xcg at each point of the measurement
xcg_tot_inch=[]
xcg_tot_inch2=[]
mass_aircraft_pounds=[]
for i in range(0,16):
    mass_aircraft_pounds=mass_aircraft_pounds+[mzfw+fuel_left_pounds[i]]
    xcg_tot_inch=xcg_tot_inch+[((xcg_zfw*mzfw)+(moment_pounds_final[i]))/mass_aircraft_pounds[i]]
    xcg_tot_inch2 = xcg_tot_inch2 + [((xcg_zfw2 * mzfw) + (moment_pounds_final[i])) / mass_aircraft_pounds[i]]
LEMAC=261.56
MAC=80.98
Time=[0,1157,1297,1426,1564,1787,1920,2239,2351,2484,2576,2741,2840,2920,3062,3166]
#####taking into account the change in xcg due to position change of the student####
xcg_pass_new=[]
xcg_tot_final_inch=xcg_tot_inch+[xcg_tot_inch2[-2]]+[xcg_tot_inch2[-1]]
xcg_tot_final_inch=np.array(xcg_tot_final_inch)-LEMAC
xcg_tot_inch=np.array(xcg_tot_inch)-LEMAC
xcg_tot_inch2=np.array(xcg_tot_inch2)-LEMAC
xcg_tot_met=xcg_tot_inch*inch_metre
xcg_tot_met2=xcg_tot_inch2*inch_metre
xcg_tot_final_met=xcg_tot_final_inch*inch_metre
Time_final=Time+[Time[-2]]+[Time[-1]]
#per_mac=xcg_tot_inch/MAC
plt.plot(Time,xcg_tot_met)
plt.show()
print(xcg_tot_final_met)
print(Time_final)
plt.scatter(Time_final,xcg_tot_final_met)
plt.show()
#print(mzfw)
#print(moment_pounds_final[0])
#print(mass_aircraft_pounds[0])
#print(per_mac)
#print(moment_pounds_final)
#print(fuel_left_pounds)
delta_xcg2=(xcg_tot_final_inch[-3]-xcg_tot_final_inch[-1])
print("Delta_xcg",delta_xcg2)