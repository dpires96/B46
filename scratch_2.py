import math as m
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# stationary measurement series 1
#mass of the observers in kg
p1=95
p2=92
o_1L=66
o_1R=61
o_2L=75
o_2R=78
o_3L=86
o_3R=68
co_ord=74
M_pass=[p1,p2,o_1L,o_1R,o_2L,o_2R,o_3L,o_3R,co_ord]
fuel=4050*0.453592
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
V_c=[249,221,192,163,130,118]    # input speeds
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
    c_l=c_l+[Lift_1/(0.5*rho0*v_eq[i]**2*s)]  #for v_eq the density used is rho0

 #calculating Cd
L_Thr=[3665.03,2995.38,2399.67,1863.38,1892.21,2208.82]
R_Thr=[3770.95,3057.27,2526.11,2015.87,2070.75,2405.27]
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
print("e",e)
print("cd0",cd_0)

alph=[1.7,2.4,3.6,5.4,8.7,10.6]
alph_rad=[]
for i in range(0,6):
    alph_rad=alph_rad+[alph[i]*3.14/180]
#print(alph_rad)

cl_alpha=(c_l[-1]-c_l[0])/(alph_rad[-1]-alph_rad[0])
#g=np.polyfit(alph_rad,c_l,1)
#g1=np.poly1d(g)
print("cl_alpha",cl_alpha)


#cd-cl
cd_plot_meas=[0.02407944, 0.02530621, 0.02705309, 0.03283982, 0.04007996, 0.06817726]
cd_cl_meas=[0.19567677, 0.24027442, 0.3451729,  0.46393924, 0.6771203,  0.91850414]
cd_x_new_meas=[0.02407944, 0.0249794,  0.02587935, 0.02677931, 0.02767926, 0.02857922,
 0.02947917, 0.03037913, 0.03127909, 0.03217904, 0.033079  , 0.03397895,
 0.03487891, 0.03577886, 0.03667882, 0.03757877, 0.03847873, 0.03937869,
 0.04027864, 0.0411786 , 0.04207855, 0.04297851, 0.04387846, 0.04477842,
 0.04567837, 0.04657833, 0.04747829, 0.04837824, 0.0492782 , 0.05017815,
 0.05107811, 0.05197806, 0.05287802, 0.05377797, 0.05467793, 0.05557789,
 0.05647784, 0.0573778 , 0.05827775, 0.05917771, 0.06007766, 0.06097762,
 0.06187757, 0.06277753, 0.06367749, 0.06457744, 0.0654774 , 0.06637735,
 0.06727731, 0.06817726]
cd_y_new_meas=[0.2057175,  0.23796059, 0.26946655, 0.30023538, 0.33026708, 0.35956165,
 0.3881191 , 0.41593941, 0.4430226 , 0.46936867, 0.4949776 , 0.51984941,
 0.54398408, 0.56738163, 0.59004206, 0.61196535, 0.63315151, 0.65360055,
 0.67331246, 0.69228724, 0.7105249 , 0.72802542, 0.74478882, 0.76081509,
 0.77610423, 0.79065625, 0.80447113, 0.81754889, 0.82988952, 0.84149302,
 0.8523594 , 0.86248864, 0.87188076, 0.88053575, 0.88845361, 0.89563435,
 0.90207795, 0.90778443, 0.91275378, 0.916986  , 0.92048109, 0.92323906,
 0.9252599 , 0.92654361, 0.92709019, 0.92689964, 0.92597197, 0.92430717,
 0.92190524, 0.91876618]

plt.figure(2)
q=np.polyfit(c_d,c_l,2)
f=np.poly1d(q)
c_d_plot=np.array(c_d)
c_d_plot=np.sort(c_d_plot)    # for linspace
c_l_plot=np.array(c_l)
c_l_plot=np.sort(c_l_plot)
x_new = np.linspace(c_d_plot[0], c_d_plot[-1])
y_new = f(x_new)
plt.plot(c_d_plot,c_l_plot,'o',x_new,y_new,color="green",label="Reference Data")
plt.plot(cd_plot_meas,cd_cl_meas,'o',cd_x_new_meas, cd_y_new_meas,color="orange",label="Measurement Data")
plt.legend()
plt.xlabel("C_D [-]")
plt.ylabel("C_L [-]")
plt.show()

#cl-alpha
alpha_plot_meas=[0.02442222, 0.03488889, 0.05931111, 0.08547778, 0.13432222, 0.18491111]
alpha_c_l_plot_meas= [0.19567677, 0.24027442, 0.3451729,  0.46393924, 0.6771203,  0.91850414]
alpha_x_new_meas= [ 0.02442222, 0.02769751, 0.03097279, 0.03424807, 0.03752336, 0.04079864,
 0.04407392, 0.04734921, 0.05062449, 0.05389977, 0.05717506, 0.06045034,
 0.06372562, 0.06700091, 0.07027619, 0.07355147, 0.07682676, 0.08010204,
 0.08337732, 0.08665261, 0.08992789, 0.09320317, 0.09647846, 0.09975374,
 0.10302902, 0.10630431, 0.10957959, 0.11285488, 0.11613016, 0.11940544,
 0.12268073, 0.12595601, 0.12923129, 0.13250658, 0.13578186, 0.13905714,
 0.14233243, 0.14560771, 0.14888299, 0.15215828, 0.15543356, 0.15870884,
 0.16198413, 0.16525941, 0.16853469, 0.17180998, 0.17508526, 0.17836054,
 0.18163583, 0.18491111]
alpha_y_new_meas= [0.19151532, 0.20621929, 0.22092326, 0.23562723, 0.2503312,  0.26503517,
 0.27973914, 0.29444311, 0.30914708, 0.32385105, 0.33855502, 0.35325899,
 0.36796296, 0.38266693, 0.3973709 , 0.41207487, 0.42677884, 0.44148281,
 0.45618678, 0.47089075, 0.48559472, 0.50029869, 0.51500266, 0.52970663,
 0.5444106 , 0.55911457, 0.57381854, 0.58852251, 0.60322648, 0.61793045,
 0.63263442, 0.64733839, 0.66204236, 0.67674633, 0.6914503 , 0.70615427,
 0.72085824, 0.73556221, 0.75026618, 0.76497015, 0.77967412, 0.79437809,
 0.80908206, 0.82378603, 0.83849   , 0.85319397, 0.86789794, 0.88260191,
 0.89730588, 0.91200985]

plt.figure(3)
q=np.polyfit(alph_rad,c_l,1)
c_l_alpha=q[0]
f=np.poly1d(q)
c_d_plot=np.array(alph_rad)
c_d_plot=np.sort(c_d_plot)    # for linspace
c_l_plot=np.array(c_l)
c_l_plot=np.sort(c_l_plot)
x_new = np.linspace(c_d_plot[0], c_d_plot[-1])
y_new = f(x_new)
plt.plot(c_d_plot,c_l_plot,'o',x_new,y_new,color="green",label="Reference Data")
plt.plot(alpha_plot_meas,alpha_c_l_plot_meas,'o',alpha_x_new_meas, alpha_y_new_meas,color="orange",label="Measurement Data")
plt.legend()
plt.xlabel("alpha [rad]")
plt.ylabel("C_L [-]")
plt.show()

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
de=[0,-0.4,-0.9,-1.5,0.4,0.6,1,0,-0.5]
Tot_Thr2=[5494.68,5875.28,6289.03,6653.86,6411.24,6541.06,6517.09]
d_eng=0.686
chrd=2.0569
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
    f_u2[i]=f_u2[i]*4.44822 # converting to nootuns
    Lift_1=Lift
    Lift_1=Lift-f_u2[i]
    Lift_req2=Lift_req2+[Lift_1]
    c_l2=c_l2+[Lift_1/(0.5*rho0*v_eq2[i]**2*s)]

#print(Lift,Lift_req2,f_u2)

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
deg_rad=3.14/180
de1=de[-1]*deg_rad
delta_xcg2=(-1.7)*0.0254   #metres
cn2=c_l2[-1]
print(rho2[-1],V_c2[-1])
cm_delta2=-1/(de1)*cn2*(delta_xcg2/chrd)
print("cm_delta=",cm_delta2)
#print(cn2)
de_str2=[]
for i in range(0,7):
    de_str2=de_str2+[de[i]-1/cm_delta2*cm_tc*(Tcs2[i]-TC2[i])]

#Data from meas:
#Ve delta e:
v_tilda1_meas= [ 69.14524483,  73.82340752,  79.03515151,  84.684758,    90.76076007,
  96.09362594, 101.00179236]
de_str2_meas= [-1.98803304, -1.49001913, -0.99164379, -0.5930726,  -0.19296171,  0.20657256,
  0.50660404]
x_new1_meas= [ 69.14524483,  69.79537846,  70.44551208,  71.0956457,   71.74577932,
  72.39591295,  73.04604657,  73.69618019,  74.34631382,  74.99644744,
  75.64658106,  76.29671468,  76.94684831,  77.59698193,  78.24711555,
  78.89724918,  79.5473828 ,  80.19751642,  80.84765005,  81.49778367,
  82.14791729,  82.79805091,  83.44818454,  84.09831816,  84.74845178,
  85.39858541,  86.04871903,  86.69885265,  87.34898628,  87.9991199,
  88.64925352,  89.29938714,  89.94952077,  90.59965439,  91.24978801,
  91.89992164,  92.55005526,  93.20018888,  93.85032251,  94.50045613,
  95.15058975,  95.80072337,  96.450857  ,  97.10099062,  97.75112424,
  98.40125787,  99.05139149,  99.70152511, 100.35165873, 101.00179236]    
y_new1_meas= [-1.95789023, -1.89350097, -1.82971554, -1.76653393, -1.70395616, -1.64198221,
 -1.58061208, -1.51984579, -1.45968332, -1.40012468, -1.34116987, -1.28281888,
 -1.22507172, -1.16792839, -1.11138889, -1.05545321, -1.00012136, -0.94539333,
 -0.89126914, -0.83774877, -0.78483223, -0.73251952, -0.68081063, -0.62970557,
 -0.57920434, -0.52930693, -0.48001335, -0.4313236 , -0.38323768, -0.33575558,
 -0.28887731, -0.24260287, -0.19693226, -0.15186547, -0.10740251, -0.06354337,
 -0.02028807,  0.02236341,  0.06441106,  0.10585489,  0.14669488,  0.18693105,
  0.2265634 ,  0.26559191,  0.3040166 ,  0.34183746,  0.37905449,  0.4156677,
  0.45167708,  0.48708263]    
#Ve Fe:
v_tilda_meas= [ 69.14524483,  73.82340752,  79.03515151,  84.684758,    90.76076007,
  96.09362594, 101.00179236]
Fe_air_eff_meas= [-39.10086568, -24.26634863, -12.64696884,   2.10147738,  32.81117267,
  64.66456428,  90.39522629]
x_new_meas= [ 69.14524483,  69.79537846,  70.44551208,  71.0956457,   71.74577932,
  72.39591295,  73.04604657,  73.69618019,  74.34631382,  74.99644744,
  75.64658106,  76.29671468,  76.94684831,  77.59698193,  78.24711555,
  78.89724918,  79.5473828 ,  80.19751642,  80.84765005,  81.49778367,
  82.14791729,  82.79805091,  83.44818454,  84.09831816,  84.74845178,
  85.39858541,  86.04871903,  86.69885265,  87.34898628,  87.9991199,
  88.64925352,  89.29938714,  89.94952077,  90.59965439,  91.24978801,
  91.89992164,  92.55005526,  93.20018888,  93.85032251,  94.50045613,
  95.15058975,  95.80072337,  96.450857  ,  97.10099062,  97.75112424,
  98.40125787,  99.05139149,  99.70152511, 100.35165873, 101.00179236]
y_new_meas= [-36.94666264, -35.89015747, -34.76833876, -33.58120653, -32.32876077,
 -31.01100148, -29.62792866, -28.17954231, -26.66584243, -25.08682902,
 -23.44250208, -21.73286161, -19.95790762, -18.11764009, -16.21205904,
 -14.24116446, -12.20495635, -10.10343471,  -7.93659954,  -5.70445084,
  -3.40698861,  -1.04421285,   1.38387644,   3.87727925,   6.4359956,
   9.06002547,  11.74936887,  14.5040258 ,  17.32399626,  20.20928025,
  23.15987777,  26.17578882,  29.2570134 ,  32.4035515 ,  35.61540314,
  38.8925683 ,  42.235047  ,  45.64283922,  49.11594497,  52.65436425,
  56.25809706,  59.9271434 ,  63.66150327,  67.46117667,  71.32616359,
  75.25646405,  79.25207803,  83.31300555,  87.43924659,  91.63080116]

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
print(slope_de_dalpha)
cm_alpha2=-slope_de_dalpha*cm_delta2
Lift_req2=Lift_req2[0:-2]     # because the last 2 values where for finding Cn
Fe_air_eff=[]
for i in range (0,7):
    Fe_air_eff= Fe_air_eff+[Fe_ref2[i]*W_s/Lift_req2[i]]
print("cm_alpha2",cm_alpha2)
#plt.scatter(alpha2,de2)
#plt.show()
#plt.scatter(v_tilda,Fe_air_eff)


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
#xcg_pass_tot=
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
#print(xcg_tot_final_met)
#print(Time_final)
plt.scatter(Time_final,xcg_tot_final_met)
plt.show()
#print(mzfw)
#print(moment_pounds_final[0])
#print(mass_aircraft_pounds[0])
#print(per_mac)
#print(moment_pounds_final)
#print(fuel_left_pounds)
print(xcg_tot_final_inch[-3]-xcg_tot_final_inch[-1])
