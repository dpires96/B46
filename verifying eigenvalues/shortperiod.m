% Citation 550 - Linear simulation


% Symmetric: phugoid 0:56 (index around 33460)  Short period 0:59 (index around 35735) 
% Asymmetric: Aper. Roll 1:01 Dutch Roll 1:03  Dutch Roll YD 1:04  Spiral 1:05

fd =  load('FTISxprt-20200311_flight2');
asym = 0; %0 if symmetric, 1 if assymetric

startindex = 35734; %(Second-9)*10
refindex = 35734; %3 seconds before

%
%     

% Stationary flight condition
hp0    = fd.flightdata.Dadc1_alt.data(refindex)*0.3048;      	  % pressure altitude in the stationary flight condition [m]
V0     = fd.flightdata.Dadc1_tas.data(refindex)*0.514444;            % true airspeed in the stationary flight condition [m/sec]
alpha0 = fd.flightdata.vane_AOA.data(refindex)*pi/180;       	  % angle of attack in the stationary flight condition [rad]
th0    = fd.flightdata.Ahrs1_Pitch.data(refindex)*pi/180;       	  % pitch angle in the stationary flight condition [rad]


% Aircraft mass
m      = 4157.174 + 2022.45125 - (fd.flightdata.lh_engine_FU.data(startindex) + fd.flightdata.rh_engine_FU.data(startindex))/2.205;         	  % mass [kg]

% aerodynamic properties

e      = 0.908;            % Oswald factor [ ]
CD0    = 0.0289;            % Zero lift drag coefficient [ ]
CLa    = 5.958761;            % Slope of CL-alpha curve [ ]

% Longitudinal stability
Cma    = -1.037;            % longitudinal stabilty [ ]
Cmde   = -1.866;            % elevator effectiveness [ ]

% Aircraft geometry

S      = 30.00;	          % wing area [m^2]
Sh     = 0.2*S;           % stabiliser area [m^2]
Sh_S   = Sh/S;	          % [ ]
lh     = 0.71*5.968;      % tail length [m]
c      = 2.0569;	  % mean aerodynamic cord [m]
lh_c   = lh/c;	          % [ ]
b      = 15.911;	  % wing span [m]
bh     = 5.791;	          % stabilser span [m]
A      = b^2/S;           % wing aspect ratio [ ]
Ah     = bh^2/Sh;         % stabilser aspect ratio [ ]
Vh_V   = 1;		  % [ ]
ih     = -2*pi/180;       % stabiliser angle of incidence [rad]


%cg

xcg = 0.25*c;


% Constant values concerning atmosphere and gravity

rho0   = 1.2250;          % air density at sea level [kg/m^3] 
lambda = -0.0065;         % temperature gradient in ISA [K/m]
Temp0  = 288.15;          % temperature at sea level in ISA [K]
R      = 287.05;          % specific gas constant [m^2/sec^2K]
g      = 9.81;            % [m/sec^2] (gravity constant)

rho    = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
W      = m*g;				                        % [N]       (aircraft weight)

% Constant values concerning aircraft inertia

muc    = m/(rho*S*c);
mub    = m/(rho*S*b);
KX2    = 0.019;
KZ2    = 0.042;
KXZ    = 0.002;
KY2    = 1.25*1.114;

% Aerodynamic constants

Cmac   = 0;                     % Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa;   		        % Wing normal force slope [ ]
CNha   = 2*pi*Ah/(Ah+2);        % Stabiliser normal force slope [ ]
depsda = 4/(A+2);               % Downwash gradient [ ]

% Lift and drag coefficient

CL = 2*W/(rho*V0^2*S);               % Lift coefficient [ ]
CD = CD0 + (CLa*alpha0)^2/(pi*A*e);  % Drag coefficient [ ]

% Stabiblity derivatives

CX0    = W*sin(th0)/(0.5*rho*V0^2*S);
CXu    = -0.095;
CXa    = 0.47966;
CXadot = +0.08330;
CXq    = -0.28170;
CXde   = -0.03728;

CZ0    = -W*cos(th0)/(0.5*rho*V0^2*S);
CZu    = -0.37616;
CZa    = -5.74340;
CZadot = -0.00350;
CZq    = -5.66290;
CZde   = -0.69612;

Cmu    = +0.06990;
Cmadot = +0.17800;
Cmq    = -8.79415;

CYb    = -0.7500;
CYbdot =  0     ;
CYp    = -0.0304;
CYr    = +0.8495;
CYda   = -0.0400;
CYdr   = +0.2300;

Clb    = -0.10260;
Clp    = -0.71085;
Clr    = +0.23760;
Clda   = -0.23088;
Cldr   = +0.03440;

Cnb    =  +0.1348;
Cnbdot =   0     ;
Cnp    =  -0.0602;
Cnr    =  -0.2061;
Cnda   =  -0.0120;
Cndr   =  -0.0939;

% Nominal speed
V = fd.flightdata.Dadc1_tas.data(startindex)*0.514444;
u = (V-V0)/V0;

if asym == 0

    Cs1 = [(-2*muc*c/V) 0 0 0; 0 (CZadot-2*muc)*c/V 0 0; 0 0 -c/V 0;
        0 Cmadot*c/V 0 (-2*muc*(KY2)*c/V)];
    Cs2 = [CXu CXa CZ0 CXq; CZu CZa -CX0 CZq+2*muc; 0 0 0 1; Cmu Cma 0 Cmq];
    Cs3 = [CXde; CZde; 0; Cmde];
    As = -Cs1\Cs2;
    Bs = -Cs1\Cs3;
    Cs = eye(4);
    Ds = zeros([4 1]);

    syss = ss(As,Bs,Cs,Ds);
end


if asym == 1
    Ca1 = [(CYbdot-2*mub)*b/V 0 0 0; 0 -0.5*b/V 0 0;
        0 0 -4*mub*KX2*b/V 4*mub*KXZ*b/V; Cnbdot*b/V 0 4*mub*KXZ*b/V -4*mub*KZ2*b/V ];
    Ca2 = [CYb CL CYp CYr-4*mub; 0 0 1 0; Clb 0 Clp Clr; Cnb 0 Cnp Cnr];
    Ca3 = [CYda CYdr; 0 0; Clda Cldr; Cnda Cndr];
    Aa = -Ca1\Ca2;
    Ba = -Ca1\Ca3;
    Ca = eye(4);
    Da = zeros([4 2]);


    sysa = ss(Aa,Ba,Ca,Da);
end



% From body axis to stab axis
alphab = fd.flightdata.vane_AOA.data(refindex)*pi/180;
thb = fd.flightdata.Ahrs1_Pitch.data(refindex)*pi/180;
qb = fd.flightdata.Ahrs1_bPitchRate.data(startindex)*pi/180;
alpha = alphab - alpha0;
th = thb - th0;

Tfinal = 21; %Seconds
des = fd.flightdata.delta_e.data;
de = des(startindex:startindex+Tfinal*10);
x0 = [u ; alpha; th; 0];
t = 0:0.1:Tfinal;
t = transpose(t);
y = lsim(syss, de*pi/180, t, x0);

truealphas = fd.flightdata.vane_AOA.data(startindex:startindex+Tfinal*10)*pi/180-alpha0;
trueths = fd.flightdata.Ahrs1_Pitch.data(startindex:startindex+Tfinal*10)*pi/180-th0;
trueqs = fd.flightdata.Ahrs1_bPitchRate.data(startindex:startindex+Tfinal*10)*pi/180-qb;

% plot(t,de);
plot(t,y(:,4));%For degrees multiply by 180/pi
hold on
plot(t,trueqs);%For degrees multiply by 180/pi
hold off