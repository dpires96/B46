% Citation 550 - Linear simulation


% Symmetric: phugoid 0:56 (index around 33460)  Short period 0:59 
% Asymmetric: Aper. Roll 1:01 Dutch Roll 1:03  Dutch Roll YD 1:04  Spiral 1:05

fd =  load('matlab');
asym = 0; %0 if symmetric, 1 if assymetric

startindex = 33460; %(Second-9)*10
refindex = 33460; %3 seconds before

%
%     

% Stationary flight condition
hp0    = fd.flightdata.Dadc1_alt.data(refindex)*0.3048;      	  % pressure altitude in the stationary flight condition [m]
V0     = 59.9;            % true airspeed in the stationary flight condition [m/sec]
alpha0 = fd.flightdata.vane_AOA.data(refindex)*pi/180;       	  % angle of attack in the stationary flight condition [rad]
th0    = fd.flightdata.Ahrs1_Pitch.data(refindex)*pi/180;       	  % pitch angle in the stationary flight condition [rad]


% Aircraft mass
m      = 4547.8;

% aerodynamic properties

e      = 0.908;            % Oswald factor [ ]
CD0    = 0.0289;            % Zero lift drag coefficient [ ]
CLa    = 5.958761;            % Slope of CL-alpha curve [ ]

% Longitudinal stability
Cma    = -0.43;            % longitudinal stabilty [ ]
Cmde   = -1.553;            % elevator effectiveness [ ]

% Aircraft geometry

S      = 24.2;	          % wing area [m^2]
Sh     = 0.2*S;           % stabiliser area [m^2]
Sh_S   = Sh/S;	          % [ ]
lh     = 5.5;      % tail length [m]
c      = 2.022;	  % mean aerodynamic cord [m]
lh_c   = lh/c;	          % [ ]
b      = 13.36;	  % wing span [m]
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

muc    = 102.7;
mub    = 15.5;
KX2    = 0.012;
KZ2    = 0.037;
KXZ    = 0.002;
KY2    = 0.98;

% Aerodynamic constants

Cmac   = 0;                     % Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa;   		        % Wing normal force slope [ ]
CNha   = 2*pi*Ah/(Ah+2);        % Stabiliser normal force slope [ ]
depsda = 4/(A+2);               % Downwash gradient [ ]

% Lift and drag coefficient

CL = 1.136;               % Lift coefficient [ ]
CD = CD0 + (CLa*alpha0)^2/(pi*A*e);  % Drag coefficient [ ]

% Stabiblity derivatives

CX0    = 0;
CXu    = -0.2199;
CXa    = 0.4653;
CXadot = 0;
CXq    = 0;
CXde   = 0;

CZ0    = -1.136;
CZu    = -2.272;
CZa    = -5.16;
CZadot = -1.43;
CZq    = -3.86;
CZde   = -0.6238;

Cmu    = 0;
Cmadot = -3.7;
Cmq    = -7.04;

CYb    = -0.9896;
CYbdot =  0     ;
CYp    = -0.087;
CYr    = +0.430;
CYda   = -0.0;
CYdr   = +0.3037;

Clb    = -0.0772;
Clp    = -0.3444;
Clr    = +0.280;
Clda   = -0.2349;
Cldr   = +0.0286;

Cnb    =  +0.1638;
Cnbdot =   0     ;
Cnp    =  -0.0108;
Cnr    =  -0.193;
Cnda   =  -0.0286;
Cndr   =  -0.1261;

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
    e=eig(As);
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
alpha = alphab - alpha0;
th = thb - th0;

Tfinal = 200; %Seconds
des = fd.flightdata.delta_e.data;
de = des(startindex:startindex+Tfinal*10);
q = fd.flightdata.Ahrs1_bPitchRate.data(startindex);
x0 = [u ; alpha; th; 0];
t = 0:0.1:Tfinal;
t = transpose(t);
y = lsim(syss, de*pi/180, t, x0); 

truealphas = fd.flightdata.vane_AOA.data(startindex:startindex+Tfinal*10)*pi/180-alpha0;
trueths = fd.flightdata.Ahrs1_Pitch.data(startindex:startindex+Tfinal*10)*pi/180-th0;
trueqs = fd.flightdata.Ahrs1_bPitchRate.data(startindex:startindex+Tfinal*10)-q;

%plot(t,de);
plot(t,y(:,3)*180/pi);%For degrees multiply by 180/pi
hold on
plot(t,trueths*180/pi); %For degrees multiply by 180/pi
hold off


