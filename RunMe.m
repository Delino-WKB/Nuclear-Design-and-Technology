% Be sure to have in the same directory also all the needed files( Data.mat
% LM.mat, MOX88.mat, MOX95.mat, MOX98.mat, Results.mat) that are essential
% for this code to copile

%% Initialization     

tic 
clear 
close all
clc
load Data.mat
load MOX88.mat
load MOX95.mat
load MOX98.mat
load LM.mat
N     = 1000;     % Vector dimension

 %% (I)                  
  %% a)               
% Task:             
% Determine the coolant mass flow rate, average velocity, and pressure drop
% along the generic subchannel specified in Figure 1. In the calculation of the
% pressure drop, neglect the inlet/outlet losses and those due to the grid
% spacers. Employ the following correlation to calculate the Darcy friction
% factor:

% a) 1 → coolant mass flow rate 
m_dot       = P_fa/(C_fl*(Tfl_out-Tfl_in));         % Total coolant mass flow rate [kg/s]

% a) 2 → Average velocity
A_fl        = L_fa^2 - 289*pi*(D/2)^2;              % Flowing area [m^2]
v_fl        = m_dot/(rho_fl*A_fl);                  % Average velocity of the coolant [m/s]
A_ch        = p^2-pi*(D/2)^2;                       % Area of a subchannel [m^2]
m_dot_ch    = m_dot*A_ch/A_fl;                      % Mass flow rate in a subchannel [kg/s]

% a) 3 → Pressure drop
Dh          = D*(4/pi*(p/D)^2-1);                   % Hydraulic diameter [m]
Re          = rho_fl*v_fl*Dh/mu_fl;                 % Reynolds number [-]
f_D         = 0.184*Re^(-0.2);                      % Darcy friction factor [-]
dp_f        = f_D*rho_fl/2*v_fl^2/Dh*H_tot;         % Frictional pressure drop along the channel [Pa] 
dp_g        = 9.806*rho_fl*H_tot;                   % Gravitational pressure gradient [Pa]
dp          = dp_f+dp_g;                            % Total pressure drop along the core [Pa]

  %% b)               
% Task:             
% Estimate the average heat transfer coefficient between the coolant and the 
% fuel pin.

% b) 1 → heat transfer coefficient
Pr          = C_fl*mu_fl/k_fl;                      % Prandtl number [-]
Pe          = Pr*Re;                                % Peclet number [-]
x           = p/D;                                  % [#] choosing Nusselt correlation number

if Pe>=80 && Pe<=4000 && x>=1.1 && x<=1.5 % First correlation
    Nu      = 0.58*(Dhy/D_pin)^.55*Pe^.45; % Nusselt number [-]
else 
    Nu      = 4.5+0.014*Pe^.8;  % Nu from second correlation
    h       = Nu*k_fl/Dh; % Heat transfer coefficient [W/mK]
end
   
 %% (II)
  %% a)               
% Task:             
% Evaluate the average linear heat rate (q_av) and the extrapolated length of
% the fuel pin.
syms H_e    % It the code blocks here, you need to install the add-on: Symbolic Math Toolbox
eq          = Tfl_in + 0.5*(2*H_e*q1_0/pi/m_dot_ch/C_fl*sin(pi*H_a/2/H_e)*(1+sin(pi*H_a/2/H_a)/sin(pi*H_a/2/H_e))) == Tfl_out ;
H_e         = vpasolve(eq, H_e);
H_e         = double(H_e);                          % Extrapolated length [m]

q1_av       = q1_0*2/H_a*H_e/pi*sin(pi*H_a/2/H_e);  % Average linear heat rate [W/m]

  %% b)               
% Task:             
% Draw the radial temperature profile at the midplane of the fuel pin active
% length (z = 0). In particular, evaluate the cladding outer and inner
% temperatures, and the fuel outer and maximum temperatures as well.

z           = 0;                                    % Midplane axial coordinate 
deltaT      = 2*H_e*q1_0/pi/m_dot_ch/C_fl*sin(pi*H_a/2/H_e); % Temperature variation along the channel
Tfl_0       = Tfl_in + 0.5*deltaT*(1+sin(pi*z./H_a)/sin(pi*H_a/2/H_e)); % Coolant T [°C]

r_fl        = linspace(R_co, 1.1*R_co,N);
Tfl         = Tfl_0 + q1_0/h/2/pi./r_fl;            % Coolant radial T profile[°C]

Tco_0       = Tfl(1);                               % Outer cladding T [°C]

% Cladding
r_c         = linspace(R_ci, R_co,N);
Tc          = Tco_0 + q1_0*log(R_co./r_c)./2./pi./k_c; % Cladding radial T profile [°C]
Tci_0       = Tc(1);                                   % Inner cladding T [°C]

% Gap
deltaG      = R_ci-R_fo;                            % Gap thickness [m] % cosidering also roughness, +10e-6 [m]
theta       = Tci_0+273.15;                         % Initialization of the average gap T [K]
e           = 1;
while e>1e-5
    k_g     = 15.8e-4*theta^0.79;                   % Gap thermal conductivity [W/mK]
    Tfo_0   = Tci_0+q1_0*deltaG/(2*pi*R_fo*k_g);    % Fuel outer T [°C]
    theta   = [theta, (Tfo_0+Tci_0)/2+273.15];      % Update of iteration variable
    e       = abs(theta(end)-theta(end-1));         % Error
    theta   = theta(end);
end

r_g         = linspace(R_fo, R_ci,N);
Tg          = Tci_0+q1_0*(R_ci-r_g)./(2*pi*r_g*k_g);% Gap radial T profile [°C]
Tfo_0       = Tg(1);                                % Outer fuel T [°C]

% Fuel
I_k88       = MOX88;                                % Points selected from the conductivity integral plot
P_Ik        = polyfit(I_k88(:,1), I_k88(:,2),6);    % Interpolation of the graph I vs. T
P_T         = polyfit(I_k88(:,2), I_k88(:,1),6);    % Auxiliary interpolation: T vs. I

r_f         = linspace(0,R_fo,N);
I_r         = -q1_0*(r_f.^2-R_fo^2)./4/pi/R_fo^2/100 + polyval(P_Ik,Tfo_0); % Radial dependence of the Conductivity integral

Tf          = polyval(P_T, I_r);                    % Extrapolation of the fuel T profile from I(r)
Tmax        = Tf(1);                                % Maximum fuel T [°C]

% Radial profile
r_profileAvg  = [r_f, r_g, r_c, r_fl];
T_profileAvg  = [Tf, Tg, Tc, Tfl];                    % Temperature profile plotted in the results

  %% c)               
% Task:             
% Draw the axial profile of the cladding temperature at the outer and inner
% radius.

% Coolant
z           = linspace(-H_a/2, H_a/2,N);            % Active length axial domain [m]
q1          = q1_0*cos(pi*z/H_e);                   % Axial linear heat rate [W/m]
Tfl_ax_avg  = Tfl_in + 0.5*deltaT*(1+sin(pi*z./H_a)/sin(pi*H_a/2/H_e)); % Coolant T axial profile [°C]

% Cladding
Tco_ax_avg  = Tfl_ax_avg + q1/h/2/pi/R_co;              % Outer cladding T axial profile [°C]
Tci_ax_avg  = Tco_ax_avg + q1*log(R_co/R_ci)/2/pi/k_c;  % Inner cladding T axial profile [°C]


  %% d)               
% Task:             
% Determine the axial coordinate z with respect to the fuel pin midplane where
% the maximum cladding outer temperature is reached, and the corresponding
% cladding outer temperature.

syms z_sym
q1          = q1_0*cos(pi*z_sym/H_e);                   % Axial linear heat rate [W/m]
Tfl_ax_max  = Tfl_in + 0.5*deltaT*(1+sin(pi*z_sym./H_a)/sin(pi*H_a/2/H_e)); % Coolant T axial profile [°C]

Tco_ax      = Tfl_ax_max + q1/h/2/pi/R_co;              % Outer cladding T axial profile [°C]

dTco        = diff(Tco_ax, z_sym)==0;                   % Derivative set to 0
zmax        = vpasolve(dTco, z_sym);
Tco_max     = subs(Tco_ax, zmax);

zmax        = double(zmax);                         % Axial position of the maximal outer cladding T [m]
Tco_max     = double(Tco_max);                      % Maximal outer cladding T [°C]

  %% e)               
% Task              
% Estimate the minimum value of the linear heat rate, at the midplane of the
% active length of the fuel pin, leading to fuel pellet cracking due to thermal
% stresses. Indicate the radial positions where the fracture stress is reached
% and where the brittle-ductile transition occurs.

% Analysis in average core conditions
kfuel       = 2.5;                                    % Fuel thermal conductivity [W/m/°C]
q3_0        = q1_0/pi/R_fo^2;                         % Volumetric heat rate 
T           = @(r) Tfo_0+q3_0.*(R_fo^2-r.^2)/4/kfuel; % Temperature profile in the pellet [°C]
Tavg        = 1/R_fo*integral(T,0,R_fo);              % Mean radial temperature in the fuel [°C]
Theta_r     = @(r) (T(r)-Tavg).*r;                    % Funcion of T used in the expression
Theta       = @(r) T(r)-Tavg;                         % Another funcion of T used in the expression

sigmaR_old  = zeros(N,1); sigmaTheta_old = zeros(N,1); sigmaZ_old = zeros(N,1); % Initialization of the stresses
i = 1;
for r = r_f    % Vectors for the radial stress profiles
    sigmaR_old(i)       = alpha_f.*E_f./(1-nu_f).*(1/R_fo^2.*integral(Theta_r,0,R_fo)-1./r.^2.*integral(Theta_r,0,r));
    sigmaTheta_old(i)   = alpha_f.*E_f./(1-nu_f).*(1/R_fo^2.*integral(Theta_r,0,R_fo)+1./r.^2.*integral(Theta_r,0,r)-Theta(r));
    sigmaZ_old(i)       = alpha_f.*E_f./(1-nu_f).*(2/R_fo^2.*integral(Theta_r,0,R_fo)-Theta(r));
    i=1+i;
end

% Function handels of the stresses of interest
sigma_r_old    = @(r) (9.81.*(13.6+T(r)./4500).*(T(r)<=1733)) + (9.81.*(66-0.03.*T(r)).*(T(r)>1733)); % breaking point 
sigmaTheta_fun = @(r) alpha_f.*E_f./(1-nu_f).*(1/R_fo^2.*integral(Theta_r,0,R_fo)+1./r.^2.*integral(Theta_r,0,r)-Theta(r));
sigmaR_fun     = @(r) alpha_f.*E_f./(1-nu_f).*(1/R_fo^2.*integral(Theta_r,0,R_fo)-1./r.^2.*integral(Theta_r,0,r));

% Radial position (for higher R, stress has occured)
R_fract     = fzero(@(r) sigma_r_old(r)-abs(sigmaTheta_fun(r)-sigmaR_fun(r)),0.5*R_fo);
             
q1_min = q1_0; e=1; % Initialization
while e>1e-5
    % Updated temperature profile
    q3_min   = q1_min/pi/R_fo^2; 
    deltaT_r = 2*H_e*q1_min/pi/m_dot_ch/C_fl*sin(pi*H_a/2/H_e);

    Tfl_r    = Tfl_in + 0.5*deltaT_r*(1+sin(pi*0/H_a)/sin(pi*H_a/2/H_e)); 
    Tco_r    = Tfl_r + q1_min/h/2/pi./R_co;
    Tci_r    = Tco_r + q1_min*log(R_co./R_ci)./2./pi./k_c;

    e = 1; theta = Tci_r+273.15; % Initialization of the average gap T [K]
    while e>1e-5
       k_g_r = 15.8e-4*theta^0.79; % Gap thermal conductivity [W/mK]
       Tfo_r = Tci_r+q1_min*deltaG/(2*pi*R_fo*k_g_r); % Fuel outer T [°C]
       theta = [theta, (Tfo_r+Tci_r)/2+273.15]; % Update of iteration variable
       e     = abs(theta(end)-theta(end-1)); % Error
       theta = theta(end);
    end

    Tfo_r    = Tci_r+q1_min*deltaG/(2*pi*R_fo*k_g_r);
    Tmax_r   = Tfo_r+q3_min*R_fo^2/4/kfuel;
    
    % Updated radial fuel temperature profile
    T_r      = @(r) Tmax_r-q3_min.*r.^2/4/kfuel;
    Tavg     = 1/R_fo*integral(T_r,0,R_fo);
    Theta_r  = @(r) (T_r(r)-Tavg).*r;
    Theta    = @(r) T_r(r)-Tavg;

    sigma_r  = @(r) (9.81.*(13.6+T_r(r)./4500).*(T_r(r)<=1733)) + (9.81.*(66-0.03.*T_r(r)).*(T_r(r)>1733));

    q3_min_new = sigma_r(R_fo)*(1-nu_f)*8*kfuel/alpha_f/E_f/R_fo^2; % From heat dependence of sigma_theta
    q1_min_new = q3_min_new*pi*R_fo^2;

    e = abs(q1_min_new-q1_min);
    q1_min = q1_min_new;
end

sigmaR  = zeros(N,1); sigmaTheta = zeros(N,1); sigmaZ = zeros(N,1); % Initialization of the stresses
i=1; 
for r = r_f % Updated radial profiles 
    sigmaR(i)     = alpha_f.*E_f./(1-nu_f).*(1/R_fo^2.*integral(Theta_r,0,R_fo)-1./r.^2.*integral(Theta_r,0,r));
    sigmaTheta(i) = alpha_f.*E_f./(1-nu_f).*(1/R_fo^2.*integral(Theta_r,0,R_fo)+1./r.^2.*integral(Theta_r,0,r)-Theta(r));
    sigmaZ(i)     = alpha_f.*E_f./(1-nu_f).*(2/R_fo^2.*integral(Theta_r,0,R_fo)-Theta(r));
    i=1+i;
end

% Brittle - Ductile transition:
NDTT   = (Tf_melt+273.15)/2;
R_NDTT = fzero(@(r) T(r)+273.15-NDTT,R_fo);


 %% (III)
  %% a)               
% Task: Hot channel outlet coolant T, axial coordinate where the maximum
% cladding T is reached, compatibility with the corrosion of the employed
% steel.

% Outer coolant T
q1_max     = 1.5*q1_av; % Maximal linear heat flux in hot channel conditions [W/m]
deltaT_HC  = 2*H_e*q1_max/pi/m_dot_ch/C_fl*sin(pi*H_a/2/H_e);

Tfl_out_HC = Tfl_in + 0.5*deltaT_HC*(1+sin(pi*H_a/(2*H_a))/sin(pi*H_a/2/H_e)); % Coolant outlet T in the hot channel [°C]

% Maximum cladding outer T
syms z
q1_HC      = q1_max*cos(pi*z/H_e); % Axial linear heat rate in the hot channel [W/m]
Tfl_ax_HC  = Tfl_in + 0.5*deltaT_HC*(1+sin(pi*z./H_a)/sin(pi*H_a/2/H_e)); % Coolant T axial profile in the hot channel [°C]
Tco_ax_HC  = Tfl_ax_HC + q1_HC/h/2/pi/R_co; % Outer cladding T axial profile in the hot channel [°C]

dTco_HC    = diff(Tco_ax_HC, z)==0; % Derivative set to 0 to find maximum
zmax       = vpasolve(dTco_HC, z);
Tco_max_HC = subs(Tco_ax_HC, zmax); 

zmax_HC    = double(zmax); % Axial position of the maximal outer cladding T in the hot channel [m]
Tco_max_HC = double(Tco_max_HC); % Maximal outer cladding T in the hot channel [°C]

% Analysis of the oxidated layer
k_ox = 1; % Thermal conductivity of the oxide layer [W/mK]
q1_ox_HC = q1_max*cos(pi*zmax_HC/H_e); % Linear heat flux at zmax [W/m]
Tfl_ax_zmax = Tfl_in + 0.5*deltaT_HC*(1+sin(pi*zmax_HC/H_a)/sin(pi*H_a/2/H_e)); % Coolant T axial profile in the hot channel [°C]

% Oxide thickness driven by cladding-oxide interface T:
% e=1; % Initialization of the error
% Tox_max_HC=Tco_max_HC; %  First guess of the T at the cladding-oxide interface [°C]
% while e>1e-5
%     deltaC = (0.26*(Tox_max_HC-380)+1)*10^(-6); % Thickness of the oxide layer [m]
%     R_ox = R_co-deltaC;  % Radial position of the cladding-oxide interface
%     Tox_max_HC = [Tox_max_HC, Tco_max_HC+q1_ox_HC*log(R_co/R_ox)/(2*pi*k_ox)]; % T at the cladding-oxide interface [°C]
%     e = abs(Tox_max_HC(end)-Tox_max_HC(end-1)); % Error
%     Tox_max_HC=Tox_max_HC(end);
% end
% Tox_max_HC % Oxide-cladding interface T @ zmax in HC conditions [°C]
% deltaC % Thickness of the oxide layer @ zmax [m]

% Oxide thickness formed both inside and outside of the cladding
% PB = 2.07; % Pilling-Bedworth ratio
% e=1;
% Tox_out=Tco_max_HC; % initial guess
% while e>1e-5
%     deltaC = (0.26*(Tox_out-380)+1)*10^(-6);
%     deltaC_int = deltaC/PB;
%     deltaC_out = deltaC-deltaC_int;
%     Tox_out_new = Tfl_ax_zmax + q1_ox_HC/2/pi/h/(R_co+deltaC_out);
%     e = abs(Tox_out_new-Tox_out);
%     Tox_out=Tox_out_new;
% end
% Tox_in = Tox_out + q1_ox_HC*log((R_co+deltaC_out)/(R_co-deltaC_int))/2/pi/k_ox;
% Tci_max_HC = Tox_in+q1_ox_HC*log((R_co-deltaC_int)/R_ci)/(2*pi*k_c); % Inner cladding T @ zmax [°C]

% Oxide thickness driven by outer oxide T:
deltaC_max = (0.26*(Tco_max_HC-380)+1)*10^(-6); % Thickness of the oxide layer [m]
R_ox       = R_co-deltaC_max;  % Radial position of the cladding-oxide interface
Tox_max_HC = Tco_max_HC+q1_ox_HC*log(R_co/R_ox)/(2*pi*k_ox); % T at the cladding-oxide interface [°C]
Tci_max_HC = Tox_max_HC+q1_ox_HC*log(R_ox/R_ci)/(2*pi*k_c); % Inner cladding T @ zmax [°C]


  %% b)               
% Task: Evaluation of an average core throttle

% Case 1:  Hot channel coolant outer T equal to average channel coolant outer T
deltaT_Th1 = 2*(Tfl_out-Tfl_in)/(1+sin(pi/2)/sin(pi*H_a/2/H_e)); 
m_dot_ch1  = 2*H_e*q1_max/pi/deltaT_Th1/C_fl*sin(pi*H_a/2/H_e); % Mass flow rate in a single channel [kg/s]

m_dot1     = m_dot_ch1*A_fl/A_ch; % Mass flow rate to satisfy case 1 [kg/s]

% Case 2:  Hot channel maximum outer cladding T equal to average channel outer cladding  T
e=1; % Initialization of the error
m_dot2 = m_dot1; % Initialization of the mass flow rate
while e>1e-5
    v_fl_Th = m_dot2/(rho_fl*A_fl); % Average velocity of the coolant with throttle [m/s]
    Re_Th   = rho_fl*v_fl_Th*Dh/mu_fl; % Reynolds number [-]
    Pe_Th   = Pr*Re_Th; % Peclet number [-]
    Nu_Th   = 4.5+0.014*Pe_Th^0.8; % Nusselt number [-]
    h_Th    = Nu_Th*k_fl/Dh; % Heat transfer coefficient in throttle conditions [W/m^2K]
    
    q1_zmax_HC   = q1_max*cos(pi*zmax_HC/H_e); % Linear heat transfer rate @ zmax [W/m]
    Tfl_zmax_avg = Tco_max-q1_zmax_HC/(2*pi*h_Th*R_co); % Coolant T at zmax [°C]
    deltaT_Th2   = 2*(Tfl_zmax_avg-Tfl_in)/(1+sin(zmax_HC*pi/H_a)/sin(pi*H_a/2/H_e));
    m_dot_ch2    = 2*H_e*q1_max/pi/deltaT_Th2/C_fl*sin(pi*H_a/2/H_e); % Mass flow rate in a single channel [kg/s]
    
    m_dot2_new   = m_dot_ch2*A_fl/A_ch;
    
    e      = abs(m_dot2_new-m_dot2);
    m_dot2 = m_dot2_new;  % Mass flow rate to satisfy case 2 [kg/s]
end
 
  %% c)               
% Task: Radial T profile at fuel midplane in hot channel conditions

% Fluid - Cladding:
z        = 0; % Midplane axial coordinate 
Tfl_0_HC = Tfl_in + 0.5*deltaT_HC*(1+sin(pi*z./H_a)/sin(pi*H_a/2/H_e)); % Coolant T [°C]
Tfl_HC   = Tfl_0_HC + q1_max/h/2/pi./r_fl; % Coolant radial T profile [°C]
Tco_0_HC = Tfl_HC(1); % Outer cladding T [°C]

% Oxide thickness driven by outer oxide T:
k_ox     = 1; % Thermal conductivity of the oxide layer [W/mK]
deltaC_0 = (0.26*(Tco_0_HC-380)+1)*10^(-6); % Thickness of the oxide layer [m]
R_ox_0   = R_co-deltaC_0; % Radial position of the cladding-oxide interface
r_ox_HC  = linspace(R_ox_0,R_co,N); 
Tox_HC   = Tco_0_HC + q1_max*log(R_co./r_ox_HC)./2./pi./k_ox; % T profile in the oxide layer [°C]
Tox_0_HC = Tox_HC(1);

% Cladding:
r_c_HC   = linspace(R_ci, R_ox_0,N);
Tc_HC    = Tox_0_HC + q1_max*log(R_ox_0./r_c_HC)./2./pi./k_c; % Cladding radial T profile [°C]
Tci_0_HC = Tc_HC(1); % Inner T of the oxide layer [°C]

% Gap:
theta = Tci_0_HC+273.15; % Initialization of the average gap T [K]
e = 1; % Initialization of the error
while e>1e-5
    k_g_HC   = 15.8e-4*theta^0.79; % Gap thermal conductivity [W/mK]
    Tfo_0_HC = Tci_0_HC+q1_max*deltaG/(2*pi*R_fo*k_g_HC); % Fuel outer T [°C]
    theta_new= (Tfo_0_HC+Tci_0_HC)/2+273.15; % Update of iteration variable
    e        = abs(theta-theta_new); % Error
    theta    = theta_new;
end
Tg_HC        = Tci_0_HC+q1_max*(R_ci-r_g)./(2*pi*r_g*k_g_HC); % Gap radial T profile [°C]
Tfo_0_HC     = Tg_HC(1); % Outer fuel T [°C]
 
% Fuel
I_r_HC  = -q1_max*(r_f.^2-R_fo^2)./4/pi/R_fo^2/100 + polyval(P_Ik,Tfo_0_HC); % Radial dependence of the Conductivity integral [W/cm]
Tf_HC   = polyval(P_T, I_r_HC); % Extrapolation of the fuel T profile from I(r) [°C]
Tmax_HC = Tf_HC(1); % Maximum fuel T [°C]

r_HC    = [r_f, r_g ,r_c_HC, r_ox_HC, r_fl];
T_HC    = [Tf_HC, Tg_HC, Tc_HC,Tox_HC, Tfl_HC]; 


  %% d)               
% Task: Axial T profile
z     = linspace(-H_a/2, H_a/2,N); % Active length axial domain [m]
q1_HC = q1_max*cos(pi*z/H_e);% Axial linear heat rate in the hot channel [W/m]

% Coolant:
Tfl_ax_HC = Tfl_in + 0.5*deltaT_HC*(1+sin(pi*z./H_a)/sin(pi*H_a/2/H_e)); % Coolant T axial profile [°C]

% Oxide:
Tco_ax_HC = Tfl_ax_HC + q1_HC/h/2/pi/R_co; % Outer cladding T axial profile [°C]
deltaC_HC = (0.26.*(Tco_ax_HC-380)+1)*10^(-6); % Oxide thickness driven by outer oxide T
R_ox_HC   = R_co-deltaC_HC; % Radial position of the cladding-oxide interface [m]
Tox_ax_HC = Tco_ax_HC+q1_HC.*log(R_co./R_ox_HC)./(2*pi*k_ox); % Cladding-oxide interface T axial profile [m]

% Cladding:
Tci_ax_HC = Tox_ax_HC+q1_HC.*log(R_ox_HC/R_ci)/(2*pi*k_c); % Axial profile of the inner oxide layer [°C]


  %% e)               
% Task: Linear heat rate leading to fuel melting, axial position where it
% first occurs and sensitivity analysis
q1_melt = zeros(1,3);
for disc = [1,2,3]
    e_q         = 1; % Initialization of the error
    q1_max_melt = q1_max; % First guess maximal linear heat transfer rate [W/m]
    while e_q>1e-5
        % Update of the axial profile
        z           = linspace(-H_a/2, H_a/2,N); % Active length axial domain [m]
        q1_melt_z   = q1_max_melt*cos(pi*z/H_e);% Axial linear heat rate in the hot channel [W/m]
        
        deltaT_melt = 2*H_e*q1_max_melt/pi/m_dot_ch/C_fl*sin(pi*H_a/2/H_e); % Delta T updated with new q1_max 
        Tfl_ax_melt = Tfl_in + 0.5*deltaT_melt*(1+sin(pi*z./H_a)/sin(pi*H_a/2/H_e)); % Coolant T axial profile [°C]
     
        Tco_ax_melt = Tfl_ax_melt + q1_melt_z/h/2/pi/R_co; % Outer cladding T axial profile [°C]
        
        for i=1:length(z)    
            % Oxide thickness driven by outer oxide T:
            deltaC_melt(i) = (0.26*(Tco_ax_melt(i)-380)+1)*10^(-6); % Axial profile of the oxide layer thickness [m]
            R_ox_melt      = R_co-deltaC_melt(i); % Radial position of the cladding-oxide interface [m]
            Tox_ax_melt(i) = Tco_ax_melt(i)+q1_melt_z(i)*log(R_co/R_ox_melt)/(2*pi*k_ox);

            Tci_ax_melt(i) = Tox_ax_melt(i)+q1_melt_z(i).*log(R_ox_melt/R_ci)/(2*pi*k_c); % Inner cladding T axial profile [°C]
            
            theta = Tci_ax_melt(i)+273.15; % Initialization of the average gap T [K]
            e     = 1; % Initialization of the error
            while e>1e-5
                k_g_melt(i)    = 15.8e-4*theta^0.79; % Gap thermal conductivity [W/mK]
                Tfo_ax_melt(i) = Tci_ax_melt(i)+q1_melt_z(i)*deltaG/(2*pi*R_fo*k_g_melt(i)); % Fuel outer T [°C]
                theta_new      = (Tfo_ax_melt(i)+Tci_ax_melt(i))/2+273.15; % Update of iteration variable
                e              = abs(theta_new-theta); % Error
                theta          = theta_new;
            end
        end
        [Tfo_max_melt, imelt]  = max(Tfo_ax_melt); % Maximal fuel outer T and relative index
        zmelt                  = z(imelt); % Axial position corresponding to the index [m]
        
        
        if disc == 1
            I_melt               = 78; % Integral of the thermal conductivity towards the melting of the fuel [W/cm]
            q1_melt(disc)        = 4*pi*(I_melt-polyval(P_Ik, Tfo_max_melt))*100; % Linear heat rate in the melting position [W/m]
            q1_max_melt_Integral = q1_melt(disc)/cos(zmelt*pi/H_e); % Corresponing value of q1 @ midplane [W/m]
            
            e_q                  = abs(q1_max_melt_Integral-q1_max_melt); % Updated value of q1_max
            q1_max_melt          = q1_max_melt_Integral;
            zmelt_Integral       = zmelt;
            
        elseif disc == 2
            x                 = 0; % Deviation from stoichiometry [-]
            wtPu              = 0.2; % Plutonium content, [Pu] (conversion from wt% to at% in function 'Magni') 
            P                 = 1-rho_f/100; % Porosity [-]
            bu                = 0; % Burn-up [GWd/tHM]
            [~,Tmelt_Magni]   = Magni(1,x,wtPu,P,bu); 
            k_Magni           = @(T) Magni(T,x,wtPu,P,bu);

            q1_melt(disc)     = 4*pi*integral(k_Magni, Tfo_max_melt+273.15,Tmelt_Magni);
            Tmelt_Magni       = Tmelt_Magni-273.15; % Conversion from [K] to [°C]
            q1_max_melt_Magni = q1_melt(disc)/cos(zmelt*pi/H_e); % Corresponing value of q1 @ midplane [W/m]
            
            e_q               = abs(q1_max_melt_Magni-q1_max_melt); % Updated value of q1_max
            q1_max_melt       = q1_max_melt_Magni;
            zmelt_Magni       = zmelt;

        elseif disc == 3
            x                = 0; % Deviation from stoichiometry [-]
            wtPu             = 0.2; % Plutonium content, [Pu]: [0, 45] at.%. 
            P                = 1-rho_f/100; % Porosity [-]
            bu               = 0; % Burn-up [GWd/tHM]
            [~,Tmelt_Even]   = Even(1,x,wtPu,P,bu);
            k_Even           = @(T) Even(T,x,wtPu,P,bu);

            q1_melt(disc)    = 4*pi*integral(k_Even, Tfo_max_melt+273.15,Tmelt_Even);
            Tmelt_Even       = Tmelt_Even-273.15; % Conversion from [K] to [°C]
            q1_max_melt_Even = q1_melt(disc)/cos(zmelt*pi/H_e); % Corresponing value of q1 @ midplane [W/m]
            
            e_q              = abs(q1_max_melt_Even-q1_max_melt); % Updated value of q1_max
            q1_max_melt      = q1_max_melt_Even;
            zmelt_Even       = zmelt;

        end
    end
end

% Summary of the results:
q1_max_melt = [q1_max_melt_Integral,q1_max_melt_Magni,q1_max_melt_Even];
zmelt       = [zmelt_Integral, zmelt_Magni, zmelt_Even];
Tmelt_comp  = [polyval(P_T, I_melt), Tmelt_Magni, Tmelt_Even];

% DI SEGUITO (DA RIGA 485 A RIGA 812) DIPENDENZA DELLE CORRELAZIONI DA: BU, Pu E X

% Sensitivity analysis on the correlations varying x, bu:
% Magni: x=[0,0.04], bu = [0, 110], [Pu]: [0, 45]; Even: x=[0,0.06], bu = [0, 60], [Pu] = [0.15, 03]
% x                = 0;
% bu_val           = linspace(0,60,10);
% wtPu_val         = linspace(0.15,0.3,10);
% q1_max_melt_sens = q1_max; % First guess maximal linear heat transfer rate [W/m]
% 
% % Initialization
% Tmelt_Magni_sens = zeros(10,10); Tmelt_Even_sens  = zeros(10,10);
% q1_max_Magni     = zeros(10,10); q1_max_Even      = zeros(10,10);
% 
% for j = 1:length(bu_val)
%     for k = 1:length(wtPu_val)
%         e_q         = 1; % Initialization of the error       
%         while e_q>1e-5
%             q1_melt_z   = q1_max_melt_sens.*cos(pi*z/H_e);% Axial linear heat rate in the hot channel [W/m]
%             deltaT_melt = 2*H_e*q1_max_melt_sens/pi/m_dot_ch/C_fl*sin(pi*H_a/2/H_e); % Delta T updated with new q1_max 
%             Tfl_ax_melt = Tfl_in + 0.5*deltaT_melt*(1+sin(pi*z./H_a)/sin(pi*H_a/2/H_e)); % Coolant T axial profile [°C]
%      
%             Tco_ax_melt = Tfl_ax_melt + q1_melt_z/h/2/pi/R_co; % Outer cladding T axial profile [°C]
%         
%             for i=1:length(z)    
%                 % Oxide thickness driven by outer oxide T:
%                 deltaC_melt(i) = (0.26*(Tco_ax_melt(i)-380)+1)*10^(-6); % Axial profile of the oxide layer thickness [m]
%                 R_ox_melt      = R_co-deltaC_melt(i); % Radial position of the cladding-oxide interface [m]
%                 Tox_ax_melt(i) = Tco_ax_melt(i)+q1_melt_z(i)*log(R_co/R_ox_melt)/(2*pi*k_ox);
%     
%                 Tci_ax_melt(i) = Tox_ax_melt(i)+q1_melt_z(i).*log(R_ox_melt/R_ci)/(2*pi*k_c); % Inner cladding T axial profile [°C]
%                 
%                 theta = Tci_ax_melt(i)+273.15; % Initialization of the average gap T [K]
%                 e     = 1; % Initialization of the error
%                 while e>1e-5
%                     k_g_melt(i)    = 15.8e-4*theta^0.79; % Gap thermal conductivity [W/mK]
%                     Tfo_ax_melt(i) = Tci_ax_melt(i)+q1_melt_z(i)*deltaG/(2*pi*R_fo*k_g_melt(i)); % Fuel outer T [°C]
%                     theta          = [theta, (Tfo_ax_melt(i)+Tci_ax_melt(i))/2+273.15]; % Update of iteration variable
%                     e              = abs(theta(end)-theta(end-1)); % Error
%                     theta          = theta(end);
%                 end
%             end
%             [Tfo_max_melt, imelt]     = max(Tfo_ax_melt); % Maximal fuel outer T and relative index
%             zmelt                     = z(imelt); % Axial position corresponding to the index [m]
%       
%             [~,Tmelt_Magni_sens(j,k)] = Magni(1,x,wtPu_val(k),P,bu_val(j));
%             k_Magni                   = @(T) Magni(T,x,wtPu_val(k),P,bu_val(j));
% 
%             q1_melt                   = 4*pi*integral(k_Magni, Tfo_max_melt+273.15,Tmelt_Magni_sens(j,k));
%             q1_max_Magni(j,k)         = q1_melt/cos(zmelt*pi/H_e);
%             e_q                       = abs(q1_max_Magni(j,k)-q1_max_melt_sens); % Updated value of q1_max
%             q1_max_melt_sens          = q1_max_Magni(j,k);  
%         end
% 
%         % Even
%         e_q         = 1; % Initialization of the error
%         while e_q>1e-5
%             q1_melt_z   = q1_max_melt_sens*cos(pi*z/H_e);% Axial linear heat rate in the hot channel [W/m]
%             deltaT_melt = 2*H_e*q1_max_melt_sens/pi/m_dot_ch/C_fl*sin(pi*H_a/2/H_e); % Delta T updated with new q1_max 
%             Tfl_ax_melt = Tfl_in + 0.5*deltaT_melt*(1+sin(pi*z./H_a)/sin(pi*H_a/2/H_e)); % Coolant T axial profile [°C]
%      
%             Tco_ax_melt = Tfl_ax_melt + q1_melt_z/h/2/pi/R_co; % Outer cladding T axial profile [°C]
%         
%             for i=1:length(z)    
%                 % Oxide thickness driven by outer oxide T:
%                 deltaC_melt(i) = (0.26*(Tco_ax_melt(i)-380)+1)*10^(-6); % Axial profile of the oxide layer thickness [m]
%                 R_ox_melt      = R_co-deltaC_melt(i); % Radial position of the cladding-oxide interface [m]
%                 Tox_ax_melt(i) = Tco_ax_melt(i)+q1_melt_z(i)*log(R_co/R_ox_melt)/(2*pi*k_ox);
%     
%                 Tci_ax_melt(i) = Tox_ax_melt(i)+q1_melt_z(i).*log(R_ox_melt/R_ci)/(2*pi*k_c); % Inner cladding T axial profile [°C]
%                 
%                 theta = Tci_ax_melt(i)+273.15; % Initialization of the average gap T [K]
%                 e     = 1; % Initialization of the error
%                 while e>1e-5
%                     k_g_melt(i)    = 15.8e-4*theta^0.79; % Gap thermal conductivity [W/mK]
%                     Tfo_ax_melt(i) = Tci_ax_melt(i)+q1_melt_z(i)*deltaG/(2*pi*R_fo*k_g_melt(i)); % Fuel outer T [°C]
%                     theta          = [theta, (Tfo_ax_melt(i)+Tci_ax_melt(i))/2+273.15]; % Update of iteration variable
%                     e              = abs(theta(end)-theta(end-1)); % Error
%                     theta          = theta(end);
%                 end
%             end
%             [Tfo_max_melt, imelt]    = max(Tfo_ax_melt); % Maximal fuel outer T and relative index
%             zmelt                    = z(imelt); % Axial position corresponding to the index [m]
%       
%             [~,Tmelt_Even_sens(j,k)] = Even(1,x,wtPu_val(k),P,bu_val(j));
%             k_Magni                  = @(T) Even(T,x,wtPu_val(k),P,bu_val(j));
% 
%             q1_melt                  = 4*pi*integral(k_Even, Tfo_max_melt+273.15,Tmelt_Even_sens(j,k));
%             q1_max_Even(j,k)         = q1_melt/cos(zmelt*pi/H_e); % Corresponing value of q1 @ midplane [W/m]
%             
%             e_q                      = abs(q1_max_Even(j,k)-q1_max_melt_sens); % Updated value of q1_max
%             q1_max_melt_sens         = q1_max_Even(j,k);      
%         end
%     end
% end
% 
% 
% figure(11)
% % Magni
% plot(bu_val, (q1_max_Magni(:,1)-q1_max_melt_Magni)/q1_max_melt_Magni*100,'Color',[0 0.4470 0.7410],'LineWidth',2) % Pu = 0.15
% hold on 
% grid on 
% grid minor
% plot(bu_val, (q1_max_Magni(:,4)-q1_max_melt_Magni)/q1_max_melt_Magni*100,'Color',[0.9290 0.6940 0.1250],'LineWidth',2) % Pu = 0.2
% plot(bu_val, (q1_max_Magni(:,end)-q1_max_melt_Magni)/q1_max_melt_Magni*100,'Color',[0.8500 0.3250 0.0980],'LineWidth',2) % Pu = 0.3
% % Even
% plot(bu_val, (q1_max_Even(:,1)-q1_max_melt_Even)/q1_max_melt_Even*100,'--','Color',[0 0.4470 0.7410],'LineWidth',2) % Pu = 0.15
% plot(bu_val, (q1_max_Even(:,4)-q1_max_melt_Even)/q1_max_melt_Even*100,'--','Color',[0.9290 0.6940 0.1250],'LineWidth',2) % Pu = 0.2
% plot(bu_val, (q1_max_Even(:,end)-q1_max_melt_Even)/q1_max_melt_Even*100,'--','Color',[0.8500 0.3250 0.0980],'LineWidth',2) % Pu = 0.3
% ylabel('$\%$ Variation of $q^{\prime}$','Interpreter','latex','FontSize',15)
% xlabel('Burnup','Interpreter','latex','FontSize',15)
% legend('${[Pu]=0.15}$','${[Pu]=0.2}$','${[Pu]=0.3}$', 'Location','northeast','Interpreter','latex','FontSize',12)
% title('x=0 fixed','Interpreter','latex','FontSize',15)
% 
% 
% % Magni: x=[0,0.04], bu = [0, 110], [Pu]: [0, 45]; Even: x=[0,0.06], bu = [0, 60], [Pu] = [0.15, 03]
% x_val            = linspace(0,0.04,10);
% bu_val           = linspace(0,60,10);
% wtPu             = 0.2;
% q1_max_melt_sens = q1_max;% First guess maximal linear heat transfer rate [W/m]
% Tmelt_Magni_sens = zeros(10,10); Tmelt_Even_sens  = zeros(10,10);
% q1_max_Magni     = zeros(10,10); q1_max_Even      = zeros(10,10);
% 
% for j = 1:length(bu_val)
%     for k = 1:length(x_val)
%         e_q         = 1; % Initialization of the error       
%         while e_q>1e-5
%             q1_melt_z   = q1_max_melt_sens.*cos(pi*z/H_e);% Axial linear heat rate in the hot channel [W/m]
%             deltaT_melt = 2*H_e*q1_max_melt_sens/pi/m_dot_ch/C_fl*sin(pi*H_a/2/H_e); % Delta T updated with new q1_max 
%             Tfl_ax_melt = Tfl_in + 0.5*deltaT_melt*(1+sin(pi*z./H_a)/sin(pi*H_a/2/H_e)); % Coolant T axial profile [°C]
%      
%             Tco_ax_melt = Tfl_ax_melt + q1_melt_z/h/2/pi/R_co; % Outer cladding T axial profile [°C]
%         
%             for i=1:length(z)    
%                 % Oxide thickness driven by outer oxide T:
%                 deltaC_melt(i) = (0.26*(Tco_ax_melt(i)-380)+1)*10^(-6); % Axial profile of the oxide layer thickness [m]
%                 R_ox_melt      = R_co-deltaC_melt(i); % Radial position of the cladding-oxide interface [m]
%                 Tox_ax_melt(i) = Tco_ax_melt(i)+q1_melt_z(i)*log(R_co/R_ox_melt)/(2*pi*k_ox);
%     
%                 Tci_ax_melt(i) = Tox_ax_melt(i)+q1_melt_z(i).*log(R_ox_melt/R_ci)/(2*pi*k_c); % Inner cladding T axial profile [°C]
%                 
%                 theta = Tci_ax_melt(i)+273.15; % Initialization of the average gap T [K]
%                 e     = 1; % Initialization of the error
%                 while e>1e-5
%                     k_g_melt(i)    = 15.8e-4*theta^0.79; % Gap thermal conductivity [W/mK]
%                     Tfo_ax_melt(i) = Tci_ax_melt(i)+q1_melt_z(i)*deltaG/(2*pi*R_fo*k_g_melt(i)); % Fuel outer T [°C]
%                     theta          = [theta, (Tfo_ax_melt(i)+Tci_ax_melt(i))/2+273.15]; % Update of iteration variable
%                     e              = abs(theta(end)-theta(end-1)); % Error
%                     theta          = theta(end);
%                 end
%             end
%             [Tfo_max_melt, imelt]     = max(Tfo_ax_melt); % Maximal fuel outer T and relative index
%             zmelt                     = z(imelt); % Axial position corresponding to the index [m]
%       
%             [~,Tmelt_Magni_sens(j,k)] = Magni(1,x_val(k),wtPu,P,bu_val(j));
%             k_Magni                   = @(T) Magni(T,x_val(k),wtPu,P,bu_val(j));
% 
%             q1_melt                   = 4*pi*integral(k_Magni, Tfo_max_melt+273.15,Tmelt_Magni_sens(j,k));
%             q1_max_Magni(j,k)         = q1_melt/cos(zmelt*pi/H_e);
%             e_q                       = abs(q1_max_Magni(j,k)-q1_max_melt_sens); % Updated value of q1_max
%             q1_max_melt_sens          = q1_max_Magni(j,k);  
%         end
% 
%         % Even
%         e_q         = 1; % Initialization of the error
%         while e_q>1e-5
%             q1_melt_z   = q1_max_melt_sens*cos(pi*z/H_e);% Axial linear heat rate in the hot channel [W/m]
%             deltaT_melt = 2*H_e*q1_max_melt_sens/pi/m_dot_ch/C_fl*sin(pi*H_a/2/H_e); % Delta T updated with new q1_max 
%             Tfl_ax_melt = Tfl_in + 0.5*deltaT_melt*(1+sin(pi*z./H_a)/sin(pi*H_a/2/H_e)); % Coolant T axial profile [°C]
%      
%             Tco_ax_melt = Tfl_ax_melt + q1_melt_z/h/2/pi/R_co; % Outer cladding T axial profile [°C]
%         
%             for i=1:length(z)    
%                 % Oxide thickness driven by outer oxide T:
%                 deltaC_melt(i) = (0.26*(Tco_ax_melt(i)-380)+1)*10^(-6); % Axial profile of the oxide layer thickness [m]
%                 R_ox_melt      = R_co-deltaC_melt(i); % Radial position of the cladding-oxide interface [m]
%                 Tox_ax_melt(i) = Tco_ax_melt(i)+q1_melt_z(i)*log(R_co/R_ox_melt)/(2*pi*k_ox);
%     
%                 Tci_ax_melt(i) = Tox_ax_melt(i)+q1_melt_z(i).*log(R_ox_melt/R_ci)/(2*pi*k_c); % Inner cladding T axial profile [°C]
%                 
%                 theta = Tci_ax_melt(i)+273.15; % Initialization of the average gap T [K]
%                 e     = 1; % Initialization of the error
%                 while e>1e-5
%                     k_g_melt(i)    = 15.8e-4*theta^0.79; % Gap thermal conductivity [W/mK]
%                     Tfo_ax_melt(i) = Tci_ax_melt(i)+q1_melt_z(i)*deltaG/(2*pi*R_fo*k_g_melt(i)); % Fuel outer T [°C]
%                     theta          = [theta, (Tfo_ax_melt(i)+Tci_ax_melt(i))/2+273.15]; % Update of iteration variable
%                     e              = abs(theta(end)-theta(end-1)); % Error
%                     theta          = theta(end);
%                 end
%             end
%             [Tfo_max_melt, imelt]    = max(Tfo_ax_melt); % Maximal fuel outer T and relative index
%             zmelt                    = z(imelt); % Axial position corresponding to the index [m]
%       
%             [~,Tmelt_Even_sens(j,k)] = Even(1,x_val(k),wtPu,P,bu_val(j));
%             k_Magni                  = @(T) Even(T,x_val(k),wtPu,P,bu_val(j));
% 
%             q1_melt                  = 4*pi*integral(k_Even, Tfo_max_melt+273.15,Tmelt_Even_sens(j,k));
%             q1_max_Even(j,k)         = q1_melt/cos(zmelt*pi/H_e); % Corresponing value of q1 @ midplane [W/m]
%             
%             e_q                      = abs(q1_max_Even(j,k)-q1_max_melt_sens); % Updated value of q1_max
%             q1_max_melt_sens         = q1_max_Even(j,k);      
%         end
%     end
% end
% 
% 
% figure(12)
% % Magni
% plot(bu_val, (q1_max_Magni(:,1)-q1_max_melt_Magni)/q1_max_melt_Magni*100,'Color',[0 0.4470 0.7410],'LineWidth',2) % x = 0
% hold on 
% grid on 
% grid minor
% plot(bu_val, (q1_max_Magni(:,5)-q1_max_melt_Magni)/q1_max_melt_Magni*100,'Color',[0.9290 0.6940 0.1250],'LineWidth',2) % x = 0.2
% plot(bu_val, (q1_max_Magni(:,end)-q1_max_melt_Magni)/q1_max_melt_Magni*100,'Color',[0.8500 0.3250 0.0980],'LineWidth',2) % x = 0.4
% % Even
% plot(bu_val, (q1_max_Even(:,1)-q1_max_melt_Even)/q1_max_melt_Even*100,'--','Color',[0 0.4470 0.7410],'LineWidth',2) % Pu = 0.15
% plot(bu_val, (q1_max_Even(:,5)-q1_max_melt_Even)/q1_max_melt_Even*100,'--','Color',[0.9290 0.6940 0.1250],'LineWidth',2) % Pu = 0.2
% plot(bu_val, (q1_max_Even(:,end)-q1_max_melt_Even)/q1_max_melt_Even*100,'--','Color',[0.8500 0.3250 0.0980],'LineWidth',2) % Pu = 0.3
% ylabel('$\%$ Variation of $q^{\prime}$','Interpreter','latex','FontSize',15)
% xlabel('Burnup','Interpreter','latex','FontSize',15)
% legend('${x=0}$','${x=0.02}$','${x=0.04}$', 'Location','northeast','Interpreter','latex','FontSize',12)
% title('wtPu = 0.2 fixed','Interpreter','latex','FontSize',15)
% 
% % Magni: x=[0,0.04], bu = [0, 110], [Pu]: [0, 45]; Even: x=[0,0.06], bu = [0, 60], [Pu] = [0.15, 03]
% x_val            = linspace(0,0.04,10);
% bu               = 0;
% wtPu_val         = linspace(0.15,0.3,10);
% q1_max_melt_sens = q1_max;% First guess maximal linear heat transfer rate [W/m]
% Tmelt_Magni_sens = zeros(10,10); Tmelt_Even_sens  = zeros(10,10);
% q1_max_Magni     = zeros(10,10); q1_max_Even      = zeros(10,10);
% 
% for j = 1:length(wtPu_val)
%     for k = 1:length(x_val)
%         e_q         = 1; % Initialization of the error       
%         while e_q>1e-5
%             q1_melt_z   = q1_max_melt_sens.*cos(pi*z/H_e);% Axial linear heat rate in the hot channel [W/m]
%             deltaT_melt = 2*H_e*q1_max_melt_sens/pi/m_dot_ch/C_fl*sin(pi*H_a/2/H_e); % Delta T updated with new q1_max 
%             Tfl_ax_melt = Tfl_in + 0.5*deltaT_melt*(1+sin(pi*z./H_a)/sin(pi*H_a/2/H_e)); % Coolant T axial profile [°C]
%      
%             Tco_ax_melt = Tfl_ax_melt + q1_melt_z/h/2/pi/R_co; % Outer cladding T axial profile [°C]
%         
%             for i=1:length(z)    
%                 % Oxide thickness driven by outer oxide T:
%                 deltaC_melt(i) = (0.26*(Tco_ax_melt(i)-380)+1)*10^(-6); % Axial profile of the oxide layer thickness [m]
%                 R_ox_melt      = R_co-deltaC_melt(i); % Radial position of the cladding-oxide interface [m]
%                 Tox_ax_melt(i) = Tco_ax_melt(i)+q1_melt_z(i)*log(R_co/R_ox_melt)/(2*pi*k_ox);
%     
%                 Tci_ax_melt(i) = Tox_ax_melt(i)+q1_melt_z(i).*log(R_ox_melt/R_ci)/(2*pi*k_c); % Inner cladding T axial profile [°C]
%                 
%                 theta = Tci_ax_melt(i)+273.15; % Initialization of the average gap T [K]
%                 e     = 1; % Initialization of the error
%                 while e>1e-5
%                     k_g_melt(i)    = 15.8e-4*theta^0.79; % Gap thermal conductivity [W/mK]
%                     Tfo_ax_melt(i) = Tci_ax_melt(i)+q1_melt_z(i)*deltaG/(2*pi*R_fo*k_g_melt(i)); % Fuel outer T [°C]
%                     theta          = [theta, (Tfo_ax_melt(i)+Tci_ax_melt(i))/2+273.15]; % Update of iteration variable
%                     e              = abs(theta(end)-theta(end-1)); % Error
%                     theta          = theta(end);
%                 end
%             end
%             [Tfo_max_melt, imelt]     = max(Tfo_ax_melt); % Maximal fuel outer T and relative index
%             zmelt                     = z(imelt); % Axial position corresponding to the index [m]
%       
%             [~,Tmelt_Magni_sens(j,k)] = Magni(1,x_val(k),wtPu_val(j),P,bu);
%             k_Magni                   = @(T) Magni(T,x_val(k),wtPu_val(j),P,bu);
% 
%             q1_melt                   = 4*pi*integral(k_Magni, Tfo_max_melt+273.15,Tmelt_Magni_sens(j,k));
%             q1_max_Magni(j,k)         = q1_melt/cos(zmelt*pi/H_e);
%             e_q                       = abs(q1_max_Magni(j,k)-q1_max_melt_sens); % Updated value of q1_max
%             q1_max_melt_sens          = q1_max_Magni(j,k);  
%         end
% 
%         % Even
%         e_q         = 1; % Initialization of the error
%         while e_q>1e-5
%             q1_melt_z   = q1_max_melt_sens*cos(pi*z/H_e);% Axial linear heat rate in the hot channel [W/m]
%             deltaT_melt = 2*H_e*q1_max_melt_sens/pi/m_dot_ch/C_fl*sin(pi*H_a/2/H_e); % Delta T updated with new q1_max 
%             Tfl_ax_melt = Tfl_in + 0.5*deltaT_melt*(1+sin(pi*z./H_a)/sin(pi*H_a/2/H_e)); % Coolant T axial profile [°C]
%      
%             Tco_ax_melt = Tfl_ax_melt + q1_melt_z/h/2/pi/R_co; % Outer cladding T axial profile [°C]
%         
%             for i=1:length(z)    
%                 % Oxide thickness driven by outer oxide T:
%                 deltaC_melt(i) = (0.26*(Tco_ax_melt(i)-380)+1)*10^(-6); % Axial profile of the oxide layer thickness [m]
%                 R_ox_melt      = R_co-deltaC_melt(i); % Radial position of the cladding-oxide interface [m]
%                 Tox_ax_melt(i) = Tco_ax_melt(i)+q1_melt_z(i)*log(R_co/R_ox_melt)/(2*pi*k_ox);
%     
%                 Tci_ax_melt(i) = Tox_ax_melt(i)+q1_melt_z(i).*log(R_ox_melt/R_ci)/(2*pi*k_c); % Inner cladding T axial profile [°C]
%                 
%                 theta = Tci_ax_melt(i)+273.15; % Initialization of the average gap T [K]
%                 e     = 1; % Initialization of the error
%                 while e>1e-5
%                     k_g_melt(i)    = 15.8e-4*theta^0.79; % Gap thermal conductivity [W/mK]
%                     Tfo_ax_melt(i) = Tci_ax_melt(i)+q1_melt_z(i)*deltaG/(2*pi*R_fo*k_g_melt(i)); % Fuel outer T [°C]
%                     theta          = [theta, (Tfo_ax_melt(i)+Tci_ax_melt(i))/2+273.15]; % Update of iteration variable
%                     e              = abs(theta(end)-theta(end-1)); % Error
%                     theta          = theta(end);
%                 end
%             end
%             [Tfo_max_melt, imelt]    = max(Tfo_ax_melt); % Maximal fuel outer T and relative index
%             zmelt                    = z(imelt); % Axial position corresponding to the index [m]
%       
%             [~,Tmelt_Even_sens(j,k)] = Even(1,x_val(k),wtPu_val(j),P,bu);
%             k_Magni                  = @(T) Even(T,x_val(k),wtPu_val(j),P,bu);
% 
%             q1_melt                  = 4*pi*integral(k_Even, Tfo_max_melt+273.15,Tmelt_Even_sens(j,k));
%             q1_max_Even(j,k)         = q1_melt/cos(zmelt*pi/H_e); % Corresponing value of q1 @ midplane [W/m]
%             
%             e_q                      = abs(q1_max_Even(j,k)-q1_max_melt_sens); % Updated value of q1_max
%             q1_max_melt_sens         = q1_max_Even(j,k);      
%         end
%     end
% end
% 
% 
% figure(13)
% % Magni
% plot(wtPu_val, (q1_max_Magni(:,1)-q1_max_melt_Magni)/q1_max_melt_Magni*100,'Color',[0 0.4470 0.7410],'LineWidth',2) % x = 0
% hold on 
% grid on 
% grid minor
% plot(wtPu_val, (q1_max_Magni(:,5)-q1_max_melt_Magni)/q1_max_melt_Magni*100,'Color',[0.9290 0.6940 0.1250],'LineWidth',2) % x = 0.2
% plot(wtPu_val, (q1_max_Magni(:,end)-q1_max_melt_Magni)/q1_max_melt_Magni*100,'Color',[0.8500 0.3250 0.0980],'LineWidth',2) % x = 0.4
% % Even
% plot(wtPu_val, (q1_max_Even(:,1)-q1_max_melt_Even)/q1_max_melt_Even*100,'--','Color',[0 0.4470 0.7410],'LineWidth',2) % Pu = 0.15
% plot(wtPu_val, (q1_max_Even(:,5)-q1_max_melt_Even)/q1_max_melt_Even*100,'--','Color',[0.9290 0.6940 0.1250],'LineWidth',2) % Pu = 0.2
% plot(wtPu_val, (q1_max_Even(:,end)-q1_max_melt_Even)/q1_max_melt_Even*100,'--','Color',[0.8500 0.3250 0.0980],'LineWidth',2) % Pu = 0.3
% ylabel('$\%$ Variation of $q^{\prime}$','Interpreter','latex','FontSize',15)
% xlabel('Plutionium content','Interpreter','latex','FontSize',15)
% legend('${x=0}$','${x=0.02}$','${x=0.04}$', 'Location','northeast','Interpreter','latex','FontSize',12)
% title('bu=0 fixed','Interpreter','latex','FontSize',15)


  %% f)               
% Task: Considering q=750 W/cm and Tmelt=2730, determine if melting occurs
% and estimate ratio of molten fuel area

q1_0_melt    = 75000; % Maximal linear heat flux in hot channel conditions [W/m]
deltaT_qmelt = 2*H_e*q1_0_melt/pi/m_dot_ch/C_fl*sin(pi*H_a/2/H_e);

% Fluid - Cladding:
Tfl_0_qmelt  = Tfl_in + 0.5*deltaT_qmelt*(1+sin(pi*0/H_a)/sin(pi*H_a/2/H_e)); % Coolant T [°C]
Tfl_qmelt    = Tfl_0_qmelt + q1_0_melt/h/2/pi./r_fl; % Coolant radial T profile [°C]
Tco_0_qmelt  = Tfl_qmelt(1); % Outer cladding T [°C]

% Oxide layer:
deltaC_qmelt = (0.26*(Tco_0_qmelt-380)+1)*10^(-6); % Oxide layer thickness @ midplane [m]
R_ox_qmelt   = R_co-deltaC_qmelt; % Radial position of the cladding-oxide interface [m]
r_ox_qmelt   = linspace(R_ox_qmelt,R_co); 
Tox_qmelt    = Tco_0_qmelt + q1_0_melt*log(R_co./r_ox_qmelt)./2./pi./k_ox; % T profile in the oxide layer [°C]
Tox_0_qmelt  = Tox_qmelt(1);

% Cladding:
r_c_qmelt    = linspace(R_ci, R_ox_qmelt,N);
Tc_qmelt     = Tox_0_qmelt + q1_0_melt*log(R_ox_qmelt./r_c_qmelt)./2./pi./k_c; % Cladding radial T profile [°C]
Tci_0_qmelt  = Tc_qmelt(1); % Inner T of the oxide layer [°C]

% Gap:
theta = Tci_0_qmelt+273.15; % Initialization of the average gap T [K]
e     = 1; % Initialization of the error
while e>1e-5
    k_g_qmelt   = 15.8e-4*theta^0.79; % Gap thermal conductivity [W/mK]
    Tfo_0_qmelt = Tci_0_qmelt+q1_0_melt*deltaG/(2*pi*R_fo*k_g_qmelt); % Fuel outer T [°C]
    theta_new   = (Tfo_0_qmelt+Tci_0_qmelt)/2+273.15; % Update of iteration variable
    e           = abs(theta-theta_new); % Error
    theta       = theta_new;
end
Tg_qmelt    = Tci_0_qmelt+q1_0_melt*(R_ci-r_g)./(2*pi*r_g*k_g_qmelt); % Gap radial T profile [°C]
Tfo_0_qmelt = Tg_qmelt(1); % Outer fuel T [°C]

% Fraction of molten area:
R_melt     = R_fo*sqrt(1-4*pi/q1_0_melt*100*(I_melt-polyval(P_Ik, Tfo_0_qmelt))); % Radial position where melting starts to occur [m]
MoltenArea = R_melt^2/R_fo^2*100; % Ratio of the molten area with respect to total fuel area [%]


  %% g)               
% Task: find linear heat rate that causes yielding of the cladding.
% Estimate axial and radial position where yielding first occurs
e=1;
q1_0_y  = 40785; % First guess
sy_corr = @(T) 536.1 -4.878e-1.*T +1.6e-3.*T.^2 -3e-6.*T.^3 +8e-10.*T.^4; % yield strength [MPa]
z       = linspace(-H_a/2, H_a/2,N);
while e>1e-5 % Update of the temperature profile
    q1_y     = q1_0_y*cos(pi*z./H_e);
    deltaT_y = 2*H_e*q1_0_y/pi/m_dot_ch/C_fl*sin(pi*H_a/2/H_e);

    Tfl_y    = Tfl_in + 0.5*deltaT_y.*(1+sin(pi*z./H_a)/sin(pi*H_a/2/H_e)); % Coolant T [°C]
    Tco_y    = Tfl_y + q1_y./h/2/pi/R_co;

    deltaC_y = (0.26.*(Tco_y-380)+1).*10^(-6); % Axial variation of the oxide thickness [m]
    R_ox_y   = R_co-deltaC_y; % Radial position of the cladding-oxide interface [m]
    
    Tox_y    = Tco_y + q1_y.*log(R_co./R_ox_y)./(2*pi*k_ox); % Cladding-oxide interface T axial profile [m]
    Tci_y    = Tox_y + q1_y.*log(R_ox_y./R_ci)/(2*pi*k_c); % Axial profile of the inner oxide layer [°C]

    sy           = sy_corr(Tci_y); % Yield strength at the inner cladding radius [MPa]
    [sy_min,i_y] = min(sy); % Value of first yielding [MPa]
    
    % Radial temperature profile
    T_rad = @(r) Tox_y(i_y) + q1_y(i_y).*log(R_ox_y(i_y)./r)./2/pi/k_c;
    T_avg = 1/(R_ox_y(i_y)-R_ci)*integral(T_rad, R_ci, R_ox_y(i_y));
    Theta = @(r) T_rad(r)-T_avg;
    Theta_r = @(r) (T_rad(r)-T_avg).*r;
    
    % Circumferential stress
    sigmaTheta_y   = alpha_c*E_c/(1-nu_c)*(2/(R_ox_y(i_y)^2-R_ci^2)*integral(Theta_r,R_ci,R_ox_y(i_y))-Theta(R_ci));
    sigmaTheta_abs = abs(sigmaTheta_y);

    err     = (sy_min-sigmaTheta_abs);
    err_rel = (sy_min-sigmaTheta_abs)/sigmaTheta_abs/10;
    
    % Update of the heat according to the value of the relative error
    q1_0_y_new = q1_0_y + err_rel*q1_0_y;
    e          = abs(err);
    q1_0_y     = q1_0_y_new ;    
end

  %% h)               
% Task: Estimate fuel pin internal pressure that causes yielding of the cladding
sy_HC                   = sy_corr(Tci_ax_HC); % Yield strength in hot channel conditions at Rci [MPa]
[sy_mech_min, i_mech_y] = min(sy_HC); % Value of first yielding [MPa]

p_y_mech                = (R_ox_HC(i_mech_y)^2-R_ci^2)/R_ox_HC(i_mech_y)^2*sy_mech_min/2; % Pressure [MPa]


  %% i)               
% Task: Estimate He content after 1 year of operation in T91, AISI-304.

% Additional Data
phi_Th       = 1e13; % Thermal neutron flux [n/cm2s]
phi_F        = 1e15; % Fast neutron flux [n/cm2s]
sigmaTh_10B  = 3.837e3*10^(-24); % Effective cross section (th. flux) for B-10 [cm2]
sigmaTh_58Ni = 4.4*10^(-24); % Effective cross section (th. flux) for Ni-58 [cm2]
sigmaTh_59Ni = 13*10^(-24); % Effective cross section (th. flux) for Ni-59 [cm2]
sigmaF_10B   = 623e-3*10^(-24); % Effective cross section (f. flux) for B-10 [cm2]
sigmaF_Fe    = 0.23e-3*10^(-24); % Effective cross section (f. flux) for Fe [cm2]
sigmaF_Cr    = 0.2e-3*10^(-24); % Effective cross section (f. flux) for Cr [cm2]
sigmaF_Ni    = 4.2e-3*10^(-24); % Effective cross section (f. flux) for Ni [cm2]
rho_AISI     = 8.32; % Density of AISI-304 [g/cm3]
comp_AISI    = [70,55.487; 19,51.996; 9,58.690; 0.8,54.938; 0.5,28.086; ...
                0.2,95.94; 0.06,12.011; 0.0005,10.811];  % Composition of AISI: [wt%, Molar mass [g/mol]]
imp_AISI     = 0.4395; % wt of impurities in AISI
rho_T91      = 7.76; % Density of T91 [g/cm3]
comp_T91     = [89.155,55.487; 8.5,51.996; 0.13,58.690; 0.45,54.938; 0.35,28.086; ...
                0.925,95.94; 0.1,12.011; 0.22,50.94];  % Composition of T91: [wt%, Molar mass [g/mol]]
imp_T91      = 0.17; % wt of impurities in T91
at_58Ni      = 68.3; % Isotopic abundance of Ni-58 [at%]
at_10B       = 19.8; % Isotopic abundance of B-10 [at%]
Nav          = 6.0221409e+23; % Avogadro's number
t_year       = 365*24*60*60; % Seconds in an one-year long irradiation period [s]

% Solution for AISI-304
Mavg_AISI    = (100-imp_AISI)/sum(comp_AISI(:,1)./comp_AISI(:,2)); % Average molar mass [g/mol]
Vmol_AISI    = Mavg_AISI/rho_AISI; % Volume of a mole [cm3/mol]
n_AISI       = Nav*rho_AISI*(comp_AISI(:,1)./100./comp_AISI(:,2)); % Number of nuclei per unit volume
n_AISI(3)    = n_AISI(3)*0.683; % Atomic fraction
n_AISI(8)    = n_AISI(8)*0.198; % Atomic fraction
comp_AISI    = [comp_AISI; imp_AISI,Mavg_AISI]; % Composition considering that impurities have M=M_avg
% Helium produced by only fast reactions
n_alpha_AISI = n_AISI(1:2).*(1-exp(-[sigmaF_Fe; sigmaF_Cr].*phi_F*t_year)); % only fast flux
% Helium produced by fast and thermal reactions (with 10B)
n_alpha_AISI = [n_alpha_AISI; n_AISI(end)*(1-exp(-(sigmaF_10B*phi_F+sigmaTh_10B*phi_Th)*t_year))];
% Assumption: Small consumption of Ni atoms, thermal and fast components
% can be decoupled:
x            = (sigmaF_Ni*phi_F+sigmaTh_58Ni*phi_Th)*t_year; % Since x<<1, decoupling is valid / Ni-58 can be assumed to be constant
% Helium produced by Ni-58 (fast reactions):
n_alpha_AISI = [n_alpha_AISI; n_AISI(3)*sigmaF_Ni*phi_F*t_year];
% Helium produced by Ni-59 (thermal reactions):
n_alpha_AISI = [n_alpha_AISI; n_AISI(3)*sigmaTh_58Ni*phi_Th*(t_year+(exp(-sigmaTh_59Ni*phi_Th*t_year)-1)/(sigmaTh_59Ni*phi_Th))]; 
% Total amount of Helium in AISI-304
nHe_AISI_tot = sum(n_alpha_AISI); % total concentration of He
Ni_appm_AISI  = n_AISI(3)*Vmol_AISI*1/Nav*10^6;
He_appm_AISI = nHe_AISI_tot*Vmol_AISI*1/Nav*10^6; % He content in atoms per ppm after an irradiation time of 1 year

% Solution for T-91
Mavg_T91     = (100-imp_T91)/sum(comp_T91(:,1)./comp_T91(:,2)); % Average molar mass [g/mol]
Vmol_T91     = Mavg_T91/rho_T91; % Volume of a mole [cm3/mol]
n_T91        = Nav*rho_T91*(comp_T91(:,1)./100./comp_T91(:,2)); % Number of nuclei per unit volume
n_T91(3)     = n_T91(3)*0.683; % Atomic fraction
comp_T91     = [comp_T91; imp_T91,Mavg_T91]; % Composition considering that impurities have M=M_avg
% Helium produced by only fast reactions
n_alpha_T91  = n_T91(1:2).*(1-exp(-[sigmaF_Fe; sigmaF_Cr].*phi_F*t_year)); 
% Helium produced by Ni-58 (fast reactions):
n_alpha_T91  = [n_alpha_T91; n_T91(3)*sigmaF_Ni*phi_F*t_year];
% Helium produced by Ni-59 (thermal reactions):
n_alpha_T91  = [n_alpha_T91; n_T91(3)*sigmaTh_58Ni*phi_Th*(t_year+(exp(-sigmaTh_59Ni*phi_Th*t_year)-1)/(sigmaTh_59Ni*phi_Th))]; 
% Total amount of Helium in T-91
nHe_T91_tot  = sum(n_alpha_T91); % total concentration of He
Ni_appm_T91  = n_T91(3)*Vmol_T91*1/Nav*10^6;
He_appm_T91  = nHe_T91_tot*Vmol_T91*1/Nav*10^6; % He content in atoms per ppm after an irradiation time of 1 year

 %% (IV)
  %% a)               
% Task: Evaluate three-zones fuel restructuring and the maximal T reached
% by the fuel varying rho0, q1_max

Teqax = 1600; Tclm  = 1800; % Imposed boundary temperatures [°C]
rho_eqax = 95; rho_clm = 98; rho_0 = rho_f; % Densities of the three zones [%TD]
A0 = 3.08e-2; B0 = 2.516e-4; C = 4.715e9; D = 16361; % Parameters for the correlation
k_z1 = @(T,P) (1./(A0+B0.*(T+273.15)) + C./(T+273.15).^2.*exp(-D./(T+273.15))).*(1-P).^2.5;

% T profile without restructuring
P0           = 1-rho_0/100; % Initial porosity
q1_res       = q1_max;
[Tf_noRes,~] = T_NoRes(q1_res, P0); % Unrestructured radial temperature profile

% First zone (unrestuctured fuel):
%R_eqax   = R_fo*sqrt(1-4*pi/q1_res*integral(@(T) k_z1(T,P0),Tf_noRes(end), Teqax));
[~,i_eqax]= min(abs(Tf_noRes-Teqax));
R_eqax    = r_f(i_eqax);
r_z1      = linspace(R_eqax,R_fo,N);
k_z1_avg  = 1/(Tf_noRes(1)-Tf_noRes(end))*integral(@(T) k_z1(T,P0),Tf_noRes(end), Tf_noRes(1));% Average thermal conductivity
T_z1      = Tf_noRes(end) + q1_res/4/pi/k_z1_avg.*(1-(r_z1./R_fo).^2);

% Second zone (equiaxed grains)
I_k95    = MOX95; % Points selected from the conductivity integral plot
P_Ik95   = polyfit(I_k95(:,2), I_k95(:,1),6); % Interpolation of the graph I vs. T
P_T95    = polyfit(I_k95(:,1), I_k95(:,2),6); % Auxiliary interpolation: T vs. I
F        = @(R_clm) q1_res/4/pi*(rho_eqax/rho_0)*(R_eqax/R_fo)^2*(1-(R_clm/R_eqax)^2-(rho_eqax-rho_0)/rho_eqax*log((R_eqax/R_clm)^2))/100-polyval(P_Ik95,Tclm)+polyval(P_Ik95,Teqax);
R_clm    = fsolve(F, R_eqax,optimset('Display','off'));  % if the code blocks here, you need the add-on: Optimization toolbox
% [~,i_clm]= min(abs(Tf_noRes-Tclm)); R_clm = r_f(i_clm);
r_z2     = linspace(R_clm,R_eqax,N);
I_r95    = q1_res/4/pi*(rho_eqax/rho_0)*(R_eqax/R_fo)^2.*(1-(r_z2./R_eqax).^2-(rho_eqax-rho_0)/rho_eqax.*log((R_eqax./r_z2).^2))./100+polyval(P_Ik95,Teqax);
T_z2     = polyval(P_T95, I_r95);

% Third zone (columnar grains):
I_k98    = MOX98; % Points selected from the conductivity integral plot
P_Ik98   = polyfit(I_k98(:,2), I_k98(:,1),6); % Interpolation of the graph I vs. T
P_T98    = polyfit(I_k98(:,1), I_k98(:,2),6); % Auxiliary interpolation: T vs. I
R_v      = sqrt((rho_clm-rho_eqax)/rho_clm*R_clm^2 + (rho_eqax-rho_0)/rho_clm*R_eqax^2);
r_z3     = linspace(R_v,R_clm,N);
I_r98    = q1_res/4/pi*(rho_clm/rho_0).*(R_clm^2-r_z3.^2)./R_fo^2.*(1-log((R_clm./r_z3).^2)/((R_clm./r_z3).^2-1))./100+polyval(P_Ik98, Tclm);
T_z3     = polyval(P_T98, I_r98);
Tmax_res = T_z3(1);

% Variable linear power, porosity
N_var        = 11; % Number of points to evaluate
rho0_vect    = linspace(80,94,N_var+4);
q1_vect      = linspace(0.9*q1_0, q1_max_melt_Integral,N_var);
Tmax_res_val = zeros(N_var); Tmax_HC_vect = zeros(N_var);
i = 1;
for rho0_val = rho0_vect
    j = 1;
    for q1_res_val = q1_vect
        % First zone (unrestuctured fuel):
        P0_val          = 1-rho0_val/100;
        [Tf_NoRes_val,~]= T_NoRes(q1_res_val,P0_val);
        Tfo_0_NoRes     = Tf_NoRes_val(end);
        R_eqax_val      = R_fo*sqrt(1-4*pi/q1_res_val*integral(@(T)k_z1(T,P0_val),Tfo_0_NoRes, Teqax));
%       [~,i_eqax]      = min(abs(Tf_NoRes_val-Teqax)); R_eqax_val     = r_f(i_eqax);  
        
        % Second zone (equiaxed grains
        F          = @(R_clm) q1_res_val/4/pi*(rho_eqax/rho0_val)*(R_eqax_val/R_fo)^2*(1-(R_clm/R_eqax_val)^2-(rho_eqax-rho0_val)/rho_eqax*log((R_eqax_val/R_clm)^2))/100-polyval(P_Ik95,Tclm)+polyval(P_Ik95,Teqax);
        R_clm_val  = fsolve(F, R_eqax_val,optimset('Display','off'));
%       [~,i_clm]  = min(abs(Tf_NoRes_val-Tclm)); R_clm_val     = r_f(i_clm);  
        I_Rclm95   = q1_res_val/4/pi*(rho_eqax/rho0_val)*(R_eqax_val/R_fo)^2.*(1-(R_clm_val/R_eqax_val).^2-(rho_eqax-rho0_val)/rho_eqax.*log((R_eqax_val./R_clm_val).^2))./100+polyval(P_Ik95,Teqax);
        Tclm_val   = Tclm; %polyval(P_T95, I_Rclm95);
        
        % Third zone (columnar grains):
        R_v_val  = sqrt((rho_clm-rho_eqax)/rho_clm*R_clm_val^2 + (rho_eqax-rho0_val)/rho_clm*R_eqax_val^2);
        I_Rv98   = q1_res_val/4/pi*(rho_clm/rho0_val)*(R_clm_val^2-R_v_val^2)/R_fo^2*(1-log((R_clm_val/R_v_val)^2)/((R_clm_val/R_v_val)^2-1))/100+polyval(P_Ik98, Tclm_val);
         
        % Comparison of fuel maximal temperatures
        Tmax_res_val(i,j)     = polyval(P_T98, I_Rv98);
        [~,Tmax_HC_vect(i,j)] = T_NoRes(q1_res_val,P0_val);
        j = j+1;
    end
    i = i+1;
end

  %% b)               
% Task: Estimation of the reduction of the gap due to differential thermal
% expansion of fuel and cladding

z          = linspace(-H_a/2, H_a/2,N); % Active length axial domain [m]
deltaG_exp = (R_ci-R_fo)*ones(N,1);
R_fo_exp = R_fo*ones(N,1); R_ci_exp = R_ci*ones(N,1); R_co_exp = R_co*ones(N,1); R_ox_exp = R_ox_HC;
T0 = 25;
Tco_ax_exp = zeros(N,1); Tox_ax_exp = zeros(N,1); Tci_ax_exp = zeros(N,1); Tfo_ax_exp = zeros(N,1); 
deltaC_exp = zeros(N,1); k_g_exp = zeros(N,1);

Tfl_ax_exp = Tfl_ax_HC; % Inner cladding T axial profile [°C]
for i=1:length(z)  
    e=1;
    while e>1e-5
        Tco_ax_exp(i) = Tfl_ax_exp(i) + q1_HC(i)/h/2/pi/R_co_exp(i); % Outer cladding T axial profile [°C]
        
        deltaC_exp(i) = (0.26.*(Tco_ax_exp(i)-380)+1)*10^(-6); % Oxide thickness driven by outer oxide T
        R_ox_exp(i)   = R_co_exp(i)-deltaC_exp(i); % Radial position of the cladding-oxide interface [m]
        Tox_ax_exp(i) = Tco_ax_exp(i)+q1_HC(i)*log(R_co_exp(i)./R_ox_exp(i))./(2*pi*k_ox); % Cladding-oxide interface T axial profile [m]
        
        Tci_ax_exp(i) = Tox_ax_exp(i)+q1_HC(i).*log(R_ox_exp(i)/R_ci_exp(i))/(2*pi*k_c); % Inner cladding T axial profile [°C]

        theta = Tci_ax_exp(i)+273.15; % Initialization of the average gap T [K]
        e_g   = 1;                    % Initialization of the error
        while e_g>1e-5
            k_g_exp(i)    = 15.8e-4*theta^0.79; % Gap thermal conductivity [W/mK]
            Tfo_ax_exp(i) = Tci_ax_exp(i)+q1_HC(i)*deltaG_exp(i)/(2*pi*R_fo_exp(i)*k_g_exp(i)); % Fuel outer T [°C]
            theta_new     = (Tfo_ax_exp(i)+Tci_ax_exp(i))/2+273.15; % Update of iteration variable
            e_g           = abs(theta-theta_new); % Error
            theta         = theta_new;
        end
        
        % T radial profiles
        r_f_exp  = linspace(0, R_fo_exp(i),N); % Updated fuel region
        I_r_exp  = -q1_max*(r_f_exp.^2-R_fo_exp(i)^2)./4/pi/R_fo_exp(i)^2/100 + polyval(P_Ik,Tfo_ax_exp(i)); % Radial dependence of the Conductivity integral [W/cm]
        Tf_exp   = polyval(P_T, I_r_exp); % Extrapolation of the fuel T profile from I(r) [°C]
        
        r_c_exp  = linspace(R_ci_exp(i), R_ox_exp(i),N); % Updated cladding region
        Tc_exp   = Tox_ax_exp(i) + q1_HC.*log(R_ox_exp(i)./r_c_exp)./2./pi./k_c; % Radial temperature profile in the expanded cladding

        r_ox_exp = linspace(R_ox_exp(i), R_co_exp(i),N); % Updated oxide region
        Tox_exp  = Tco_ax_exp(i)+q1_HC(i)*log(R_co_exp(i)./R_ox_exp(i))./(2*pi*k_ox); % Radial temperature profile in the expanded oxide layer

        % Mean temperature values:
        Tf_exp_avg    = mean(Tf_exp);
        Tc_exp_avg    = mean(Tc_exp);
        Tox_exp_avg   = mean(Tox_exp);
        
        % Update of the radii due to thermal expansion
        R_fo_exp_new  = R_fo + alpha_f*(Tf_exp_avg-T0)*R_fo;
        R_ci_exp_new  = R_ci + alpha_c*(Tc_exp_avg-T0)*R_ci;
        R_ox_exp_new  = R_ox_HC(i) + alpha_c*(Tox_exp_avg-T0)*R_ox_HC(i);
        R_co_exp_new  = R_ox_exp_new +deltaC_exp(i);
        
        % New gap width
        deltaG_new    = R_ci_exp_new-R_fo_exp_new;
        e             = abs((deltaG_new-deltaG_exp(i)));
        
        % Expanded geometry characteristics:
        deltaG_exp(i) = deltaG_new;
        R_fo_exp(i)   = R_fo_exp_new;
        R_ci_exp(i)   = R_ci_exp_new;
        R_ox_exp(i)   = R_ox_exp_new;
        R_co_exp(i)   = R_co_exp_new;

    end

    % Radial fuel temperature profile at the midplane
    if i == length(z)/2
        A0 = 3.08e-2; B0 = 2.516e-4; C = 4.715e9; D = 16361;
        k_fuel = @(T,P) (1./(A0+B0.*(T+273.15)) + C./(T+273.15).^2.*exp(-D./(T+273.15))).*(1-P).^2.5;
        q3 = q1_max/pi/R_fo_exp(i)^2;
        e=1;
        Tmax_exp = Tmax_HC;
        Tfo_0_exp = Tf_exp(end);
        while e>1e-5
            k_exp = 1/(Tmax_exp-Tfo_0_exp)*integral(@(T) k_fuel(T,P),Tfo_0_exp,Tmax_exp);
            Tmax_new = Tfo_0_HC + q3*R_fo^2/4/k_exp;
            e = abs(Tmax_new-Tmax_exp);
            Tmax_exp = Tmax_new;
        end
        rf_0_exp = linspace(0,R_fo_exp(i),N);
        Tf_0_exp = Tmax_exp - q3.*rf_0_exp.^2/4/k_exp;
    end
end

 %% (V)
  %% a)
  %% b)                 
% Task: Evaluation of thermal creep and irradiation
t_EL       = 24*365*4; % Expected life
p_int      = 12;       % External pressure
PLM        = LM(:,1)*1000;              % Larson Miller parameter                  
sigmaPLM   = 10.^LM(:,2);               % Applied stress [MPa] 
P_PLM      = polyfit(PLM,sigmaPLM,6);   % Interpolation of the Kakasaki curve
P_sigmaPLM = polyfit(sigmaPLM,PLM,6);   % Inverse interpolation

% Mariotte approach - Stresses in case of external pressure
% Case 1: no oxide
Tco_ax_noOx             = Tco_ax_HC;
Tci_ax_noOx             = Tco_ax_noOx+q1_HC.*log(R_co/R_ci)/(2*pi*k_c); % Temperature profile without the oxide layer
[Tc_creep_noOx, i_noOx] = max((Tci_ax_noOx+Tco_ax_noOx)/2+273.15); % Maximal cladding T
thick_noOx              = R_co-R_ci; % Cladding thickness without oxide
% Stress profiles
T_noOx = @(r) Tci_ax_noOx(i_noOx) + q1_HC(i_noOx).*log(R_co./r)./2/pi/k_c;
T_avg_noOx = 1/(R_co-R_ci)*integral(T_noOx, R_ci, R_co);
Theta_noOx = @(r) T_noOx(r)-T_avg_noOx;
Theta_r_noOx = @(r) (T_noOx(r)-T_avg_noOx).*r;
i = 1;
sigmaTheta_th=zeros(N,1);sigmaR_th=zeros(N,1);sigmaZ_th=zeros(N,1);sigmaTheta_mec=zeros(N,1);sigmaR_mec=zeros(N,1);sigmaZ_mec=zeros(N,1);
for r_noOx = linspace(R_ci,R_co,N)
    sigmaTheta_th(i) = alpha_c*E_c/(1-nu_c)*((1+(R_ci./r_noOx)^2)./(R_co^2-R_ci^2)*integral(Theta_r_noOx,R_ci,R_co)+1/r_noOx^2*integral(Theta_r_noOx,R_ci,r_noOx)-Theta_noOx(r_noOx));
    sigmaR_th(i)     = alpha_c*E_c/(1-nu_c)*((1-(R_ci./r_noOx)^2)./(R_co^2-R_ci^2)*integral(Theta_r_noOx,R_ci,R_co)-1/r_noOx^2*integral(Theta_r_noOx,R_ci,r_noOx));
    sigmaZ_th(i)     = alpha_c*E_c/(1-nu_c)*(2/(R_co^2-R_ci^2)*integral(Theta_r_noOx,R_ci,R_co)-Theta_noOx(r_noOx));
    sigmaTheta_mec(i)= R_ci^2/(R_co^2-R_ci^2)*p_int*(R_co^2/r_noOx^2+1);
    sigmaR_mec(i)    = -R_ci^2/(R_co^2-R_ci^2)*p_int*(R_co^2/r_noOx^2-1);
    sigmaZ_mec(i)    = 0.5*p_int*R_ci/thick_noOx; % Mariotte
    i = i+1;
end

sigmaR_noOx     = mean(sigmaTheta_th+sigmaTheta_mec);
sigmaTheta_noOx = mean(sigmaR_th+sigmaR_mec);
sigmaZ_noOx     = mean(sigmaZ_th+sigmaZ_mec);

sigma_vM_noOx   = 1/sqrt(2)*sqrt((sigmaR_noOx-sigmaTheta_noOx)^2+(sigmaR_noOx-sigmaZ_noOx)^2+(sigmaTheta_noOx-sigmaZ_noOx)^2); % Von Mises equivalent stress
PLM_noOx        = polyval(P_sigmaPLM, sigma_vM_noOx); % Larson-Miller parameter, obtained from graph

% Time margin:  
t_rupt_noOx     = 10^(PLM_noOx/Tc_creep_noOx-29.1146); % from Tc_creep_avg*(29.1146+log10(trupt_ox))==PLM
margin_t_noOx   = t_rupt_noOx/t_EL;

% Temperature margin:
Trupt_noOx      = PLM_noOx/(29.1146+log10(t_EL));
margin_T_noOx   = Trupt_noOx-Tc_creep_noOx;

% Stress margin:
PLM_margin_noOx      = Tc_creep_noOx*(29.1146+log10(t_EL)); % Result: 29.7388
sigma_vM_margin_noOx = polyval(P_PLM, PLM_margin_noOx);
margin_stress_noOx = sigma_vM_margin_noOx/sigma_vM_noOx; 


% Case 2: with Oxide layer 
[Tc_creep_Ox, i_Ox] = max((Tci_ax_HC+Tox_ax_HC)/2+273.15); % Maximal cladding T
thick_Ox      =  R_ox_HC(i_Ox)-R_ci;
% Stress profiles
T_Ox       = @(r) Tci_ax_HC(i_noOx) + q1_HC(i_Ox).*log(R_ox_HC(i_Ox)./r)./2/pi/k_c;
T_avg_Ox   = 1/(R_ox_HC(i_Ox)-R_ci)*integral(T_Ox, R_ci, R_ox_HC(i_Ox));
Theta_Ox   = @(r) T_Ox(r)-T_avg_Ox;
Theta_r_Ox = @(r) (T_Ox(r)-T_avg_Ox).*r;
i = 1;
sigmaTheta_th=zeros(N,1);sigmaR_th=zeros(N,1);sigmaZ_th=zeros(N,1);sigmaTheta_mec=zeros(N,1);sigmaR_mec=zeros(N,1);sigmaZ_mec=zeros(N,1);
for r_Ox = linspace(R_ci,R_ox_HC(i_Ox),N)
    sigmaTheta_th(i) = alpha_c*E_c/(1-nu_c)*((1+(R_ci./r_Ox)^2)./(R_ox_HC(i_Ox)^2-R_ci^2)*integral(Theta_r_Ox,R_ci,R_ox_HC(i_Ox))+1/r_Ox^2*integral(Theta_r_Ox,R_ci,r_Ox)-Theta_Ox(r_Ox));
    sigmaR_th(i)     = alpha_c*E_c/(1-nu_c)*((1-(R_ci./r_Ox)^2)./(R_ox_HC(i_Ox)^2-R_ci^2)*integral(Theta_r_Ox,R_ci,R_ox_HC(i_Ox))-1/r_Ox^2*integral(Theta_r_Ox,R_ci,r_Ox));
    sigmaZ_th(i)     = alpha_c*E_c/(1-nu_c)*(2/(R_ox_HC(i_Ox)^2-R_ci^2)*integral(Theta_r_Ox,R_ci,R_ox_HC(i_Ox))-Theta_Ox(r_Ox));
    sigmaTheta_mec(i)= R_ci^2/(R_ox_HC(i_Ox)^2-R_ci^2)*p_int*(R_ox_HC(i_Ox)^2/r_Ox^2+1);
    sigmaR_mec(i)    = -R_ci^2/(R_ox_HC(i_Ox)^2-R_ci^2)*p_int*(R_ox_HC(i_Ox)^2/r_Ox^2-1);
    sigmaZ_mec(i)    = 0.5*p_int*R_ci/thick_Ox; % Mariotte
    i = i+1;
end

sigmaR_Ox     = mean(sigmaTheta_th+sigmaTheta_mec);
sigmaTheta_Ox = mean(sigmaR_th+sigmaR_mec);
sigmaZ_Ox     = mean(sigmaZ_th+sigmaZ_mec);
 
sigma_vM_Ox   = 1/sqrt(2)*sqrt((sigmaR_Ox-sigmaTheta_Ox)^2+(sigmaR_Ox-sigmaZ_Ox)^2+(sigmaTheta_Ox-sigmaZ_Ox)^2); % Von Mises equivalent stress
PLM_Ox        = polyval(P_sigmaPLM, sigma_vM_Ox); % Larson-Miller parameter, obtained from graph

% Time margin: 
t_rupt_Ox     = 10^(PLM_Ox/Tc_creep_Ox-29.1146); % from Tc_creep_avg*(29.1146+log10(trupt_ox))==PLM
margin_t_Ox   = t_rupt_Ox/t_EL;

% Temperature margin:
Trupt_Ox      = PLM_Ox/(29.1146+log10(t_EL));
margin_T_Ox   = Trupt_Ox-Tc_creep_Ox;

% Stress margin:
PLM_margin_Ox      = Tc_creep_Ox*(29.1146+log10(t_EL)); % Result: 29.7388
sigma_vM_margin_Ox = polyval(P_PLM, PLM_margin_Ox);
margin_stress_Ox = sigma_vM_margin_Ox/sigma_vM_Ox; 

% save('Results.mat')
%% Results              
clc
%load Results.mat
   %% (I)
    %% a)
    %% b)
   %% (II)
    %% a)
    %% b)               
       
figure('Position', [1920/2, 950, 1200,500])
plot(1000*r_profileAvg,T_profileAvg,'LineWidth',2)
hold on
grid on 
grid minor
ylabel('Temperature [$^{\circ}$C]','Interpreter','latex','FontSize',15)
xlabel('Radius [mm]','Interpreter','latex','FontSize',15)
plot(1000*R_co*ones(100,1),linspace(Tfl_0,Tmax,100),'k--')
plot(1000*R_ci*ones(100,1),linspace(Tfl_0,Tmax,100),'k--')
plot(1000*R_fo*ones(100,1),linspace(Tfl_0,Tmax,100),'k--')
text(4.4,2000,'Coolant','Interpreter','latex','FontSize',15)
text(3.75,2000,'Cladding','Interpreter','latex','FontSize',15)
text(3.2,2000,'Gap $\rightarrow$','Interpreter','latex','FontSize',15)
text(2,2000,'Fuel','Interpreter','latex','FontSize',15)

figure(10)
plot(linspace(0,2800),polyval(P_Ik,linspace(0,2800)),'LineWidth',2)
hold on
grid on 
grid minor
plot(linspace(0,2800),polyval(P_Ik95,linspace(0,2800)),'LineWidth',2)
plot(linspace(0,2800),polyval(P_Ik98,linspace(0,2800)),'LineWidth',2)
ylabel('Integral of thermal conductivity [W/cm]','Interpreter','latex','FontSize',15)
xlabel('Temperature [$^{\circ}$C]','Interpreter','latex','FontSize',15)
xlim([0,2800])
legend('MOX, 88\%','MOX, 95\%','MOX, 98\%','Location','NorthWest','Interpreter','latex','FontSize',12)


    %% c)             
figure(2)
plot(Tfl_ax_avg, z, 'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on
grid on
grid minor
plot(Tco_ax_avg, z, 'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
plot(Tci_ax_avg, z, 'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(Tco_max,zmax,'h','color',[0.6350 0.0780 0.1840],'linewidth',2)
text(530,zmax,'$T_{max}$','Interpreter','latex','FontSize',15)
legend('Coolant T','Outer Cladding T', 'Inner Cladding T','Location','NorthWest','Interpreter','latex','FontSize',12)
ylabel('Axial position [m]','Interpreter','latex','FontSize',15)
xlabel('Temperature [$^{\circ}$C]','Interpreter','latex','FontSize',15)

    %% d)
    %% e)             

figure('Position', [1920/2, 950, 1200,500])
subplot(1,2,2)
plot(1000*r_f, sigmaR, 'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
hold on
grid on
grid minor
plot(1000*r_f, sigmaTheta, 'Color',[0 0.4470 0.7410],'LineWidth',2)
plot(1000*r_f, sigmaZ, 'Color', [0.9290 0.6940 0.1250],'LineWidth',1)
plot(1000*r_f, sigma_r(r_f), 'Color',[0.6350 0.0780 0.1840],'LineWidth',2)
title('Minimum linear heat rate','Interpreter','latex','FontSize',15)
ylabel('Stress [MPa]','Interpreter','latex','FontSize',15)
xlabel('Radial position [mm]','Interpreter','latex','FontSize',15)
legend('Radial, ${\sigma_r}$','Tangential, ${\sigma_\theta}$', 'Axial, ${\sigma_z}$', 'Fracture','Location','southeast','Interpreter','latex','FontSize',12)

subplot(1,2,1)
plot(1000*r_f, sigmaR_old, 'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
hold on
grid on
grid minor
plot(1000*r_f, sigmaTheta_old, 'Color',[0 0.4470 0.7410],'LineWidth',2)
plot(1000*r_f, sigmaZ_old, 'Color', [0.9290 0.6940 0.1250],'LineWidth',1)
plot(1000*r_f, sigma_r_old(r_f), 'Color',[0.6350 0.0780 0.1840],'LineWidth',2)
title('Average linear heat rate','Interpreter','latex','FontSize',15)
ylabel('Stress [MPa]','Interpreter','latex','FontSize',15)
xlabel('Radial position [mm]','Interpreter','latex','FontSize',15)
legend('Radial, ${\sigma_r}$','Tangential, ${\sigma_\theta}$', 'Axial, ${\sigma_z}$', 'Fracture','Location','southeast','Interpreter','latex','FontSize',12)


   %% (III)
    %% a)
    %% b)
    %% c)             
figure('Position', [1920/2, 950, 1300,500])
plot(1000*r_HC,T_HC,'LineWidth',2)
hold on
grid on 
grid minor
ylabel('Temperature [$^{\circ}$C]','Interpreter','latex','FontSize',15)
xlabel('Radius [mm]','Interpreter','latex','FontSize',15)
plot(1000*R_co*ones(100,1),linspace(Tfl_0,Tmax,100),'k--')
plot(1000*R_ox*ones(100,1),linspace(Tfl_0,Tmax,100),'k--')
plot(1000*R_ci*ones(100,1),linspace(Tfl_0,Tmax,100),'k--')
plot(1000*R_fo*ones(100,1),linspace(Tfl_0,Tmax,100),'k--')
text(4.4,2000,'Coolant','Interpreter','latex','FontSize',15)
text(3.75,2000,'Cladding','Interpreter','latex','FontSize',15)
text(3.2,2000,'Gap $\rightarrow$','Interpreter','latex','FontSize',15)
text(2,2200,'Fuel','Interpreter','latex','FontSize',15)
text(4.1,2200,'Oxide layer','Interpreter','latex','FontSize',15)
text(4.2,2125,'$\downarrow$','Interpreter','latex','FontSize',15)

    %% d)             
figure(5)
plot(Tfl_ax_HC, z, 'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on
grid on
grid minor
plot(Tco_ax_HC, z, 'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
plot(Tox_ax_HC, z, 'k','LineWidth',2)
plot(Tci_ax_HC, z, 'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(Tco_max_HC,zmax_HC,'h','color',[0.6350 0.0780 0.1840],'linewidth',2)
text(545,zmax_HC,'$T_{max}$','Interpreter','latex','FontSize',15)
legend('Coolant T', 'Outer Oxide T', 'Oxide-Cladding interface T', 'Inner Cladding T','Location','northwest','Interpreter','latex','FontSize',15)
ylabel('Axial position [m]','Interpreter','latex','FontSize',15)
xlabel('Temperature [$^{\circ}$C]','Interpreter','latex','FontSize',15)

    %% e)
    %% f)             
figure(6)
rCerchio        = R_melt;            
tCerchio        = 0:pi/24:2*pi;
xCerchio        = rCerchio*cos(tCerchio);
yCerchio        = rCerchio*sin(tCerchio);
fill(xCerchio,yCerchio,'')				
hold on
grid on
grid minor
rCerchio        = R_fo;            
tCerchio        = 0:pi/24:2*pi;
xCerchio        = rCerchio*cos(tCerchio);
yCerchio        = rCerchio*sin(tCerchio);
plot(xCerchio,yCerchio,'linewidth',2)
legend('Molten area',' Fuel area')
legend('location','northwest','Interpreter','latex','FontSize',15)
%title('Molten area section view')
xlabel('Radial position [m]','Interpreter','latex','FontSize',15)
ylabel('Radial position [m]','Interpreter','latex','FontSize',15)
axis ([-rCerchio rCerchio -rCerchio rCerchio]*1.3)
axis equal

    %% g)
    %% h)
    %% i)
   %% (IV)
    %% a)             
% Restructured fuel in HC conditions
figure('Position', [1920/2, 950, 1920,500])
hold on
plot(1000*r_z2, T_z2,'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
plot(1000*r_z3, T_z3,'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(1000*r_z1,T_z1,'Color',[0 0.4470 0.7410],'LineWidth',2)
plot(1000*r_f,Tf_noRes,'--','LineWidth',2)% Profile of unrestructured fuel
grid on 
grid minor
ylabel('Temperature [$^{\circ}$C]','Interpreter','latex','FontSize',15)
xlabel('Radius [mm]','Interpreter','latex','FontSize',15)
plot(1000*R_v*ones(100,1),linspace(1000,Tmax_HC,100),'k--')
plot(1000*R_eqax*ones(100,1),linspace(1000,Tmax_HC,100),'k--')
plot(1000*R_clm*ones(100,1),linspace(1000,Tmax_HC,100),'k--')
plot(1000*R_fo*ones(100,1),linspace(1000,Tmax_HC,100),'k--')
text(3,Tmax_res,'As-fabricated fuel','Interpreter','latex','FontSize',15)
text(2.56,2410,'Equiaxed grains','Interpreter','latex','FontSize',15)
text(2,Tmax_res+20,'Columnar grains','Interpreter','latex','FontSize',15)
text(0.3,Tmax_res,'Central void','Interpreter','latex','FontSize',15)

% Restructured fuel in variable conditions
% Variation of q1_max
figure('Position', [1920/2, 950, 1920,500])
subplot(1,2,1)
plot(q1_vect, (Tmax_res_val(1,:)-Tmax_HC_vect(1,:))./Tmax_HC_vect(1,:)*100,'Color',[0 0.4470 0.7410],'LineWidth',2) % rho = 80
hold on 
grid on 
grid minor
plot(q1_vect, (Tmax_res_val(9,:)-Tmax_HC_vect(9,:))./Tmax_HC_vect(9,:)*100,'Color',[0.9290 0.6940 0.1250],'LineWidth',2) % rho = 88
plot(q1_vect, (Tmax_res_val(end,:)-Tmax_HC_vect(end,:))./Tmax_HC_vect(11,:)*100,'Color',[0.8500 0.3250 0.0980],'LineWidth',2) % rho = 90
ylabel('\% Variation of $T_{max}$','Interpreter','latex','FontSize',15)
xlabel('Linear power rate [W/m]','Interpreter','latex','FontSize',15)
legend('${\rho_0 = 80 \%}$','${\rho_0 = 88 \%}$', '${\rho_0 = 94 \%}$','Location','northeast','Interpreter','latex','FontSize',15)

% Variation of P0
subplot(1,2,2)
plot(rho0_vect, (Tmax_res_val(:,1)-Tmax_HC_vect(:,1))./Tmax_HC_vect(:,1)*100,'Color',[0 0.4470 0.7410],'LineWidth',2) % rho = 80
hold on 
grid on 
grid minor
plot(rho0_vect, (Tmax_res_val(:,5)-Tmax_HC_vect(:,5))./Tmax_HC_vect(:,5)*100,'Color',[0.9290 0.6940 0.1250],'LineWidth',2) % rho = 88
plot(rho0_vect, (Tmax_res_val(:,11)-Tmax_HC_vect(:,11))./Tmax_HC_vect(:,11)*100,'Color',[0.8500 0.3250 0.0980],'LineWidth',2) % rho = 90
ylabel('\% Variation of $T_{max}$','Interpreter','latex','FontSize',15)
xlabel('Density [\%TD]','Interpreter','latex','FontSize',15)
legend('${\dot{q}_{max} = 26.03}$','${\dot{q}_{max} = 31.2}$','${\dot{q}_{max} = 38.96  }$', 'Location','northwest','Interpreter','latex','FontSize',15)

    %% b)             

figure('Position', [1920/2, 950, 1920/2,500])
plot(1000*R_fo_exp,z,'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on 
grid on 
grid minor
plot(1000*R_ci_exp,z,'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(1000*R_ox_exp,z,'k','LineWidth',2)
plot(1000*R_co_exp,z,'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
plot(1000*R_fo*ones(N,1),z,'--','Color',[0 0.4470 0.7410],'LineWidth',2)
plot(1000*R_ci*ones(N,1),z,'--','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(1000*R_ox_HC,z,'k--','LineWidth',2)
plot(1000*R_co*ones(N,1),z,'--','Color',[0.9290 0.6940 0.1250],'LineWidth',2)
xlabel('Radial position [mm]','Interpreter','latex','FontSize',15)
ylabel('Axial position [m]','Interpreter','latex','FontSize',15)
legend('Outer fuel radius', 'Inner cladding radius','Cladding-Oxide interface','Outer cladding radius','Location','northeastoutside','Interpreter','latex','FontSize',12)

    %% c)             
 % Task : Compare the temperature profiles of IVa and IVb
 % Restructured fuel in HC conditions
r_res = [r_z3,r_z2,r_z1];
T_res_prof = [T_z3,T_z2,T_z1];
figure('Position', [1920/2, 950, 1920,500])
plot(1000*r_res, T_res_prof,'LineWidth',2)
hold on
plot(1000*rf_0_exp, Tf_0_exp,'LineWidth',2)
plot(1000*r_f,Tf_noRes,'LineWidth',2)% Profile of unrestructured fuel
grid on 
grid minor
legend('Restructured fuel','Expanded fuel', 'As-fabricated fuel','Location','northeast','Interpreter','latex','FontSize',15)
ylabel('Temperature [$^{\circ}$C]','Interpreter','latex','FontSize',15)
xlabel('Radius [mm]','Interpreter','latex','FontSize',15)
plot(1000*R_fo*ones(100,1),linspace(800,2600,100),'k--','HandleVisibility','off')
plot(1000*rf_0_exp(end).*ones(100,1),linspace(800,2600,100),'k--','HandleVisibility','off')
plot(1000*R_v*ones(100,1),linspace(800,2600,100),'k--','HandleVisibility','off')
plot(1000*R_eqax*ones(100,1),linspace(800,2600,100),'k--','HandleVisibility','off')
plot(1000*R_clm*ones(100,1),linspace(800,2600,100),'k--','HandleVisibility','off')

   %% (V)
    %% b)             
figure('Position', [1920/2, 950, 1920/2,500])
semilogy(PLM,sigmaPLM,'LineWidth',2)
grid on
grid minor
ylabel('Applied stress [MPa]','Interpreter','latex','FontSize',15)
xlabel('Larson-Miller Parameter','Interpreter','latex','FontSize',15)
hold on
plot([PLM_noOx,PLM_noOx], [10, sigma_vM_noOx],':','Color',[0.4660 0.6740 0.1880],'linewidth',2)
plot([PLM_Ox,PLM_Ox], [10, sigma_vM_Ox],':','Color',[0.4940 0.1840 0.5560],'linewidth',2)
plot([23000,PLM_noOx], [sigma_vM_noOx, sigma_vM_noOx],':','Color',[0.4660 0.6740 0.1880],'linewidth',2)
plot([23000,PLM_Ox], [sigma_vM_Ox, sigma_vM_Ox],':','Color',[0.4940 0.1840 0.5560],'linewidth',2)
legend(' Kawasaki et al.,1991',' Without Oxide',' With Oxide')
legend('Location','northeast','Interpreter','latex','FontSize',15)

toc
% close all

 %% Functions
   %% (III) e)        
function [k,Tmelt] = Magni(T,x,wtPu,P,bu)
% Thermal Conductivity
% Validity range: 
%   Temperature, T: [500, 2700] K.
%   Deviation from stoichiometry, x: [0, 0.04] (hypo-stoichiometry). 
%   Plutonium content, [Pu]: [0, 45] at.%.
%   Porosity, p: [0, 7] %.
%   Burn-up, bu: [0, 130] GWd/tHM.

A0  = 0.01926;
Ax  = 1.06e-6;
APu = 2.63e-8;
B0  = 2.39e-4;
BPu = 1.37e-13;
D   = 5.27e9;
E   = 17109.5;

k_inf = 1.755;
phi   = 128.75;

% Conversion to atomic percentage
Z_Pu = 239.0521634;
Z_U  = 238.05078826;
atPu = (wtPu/Z_Pu)/(wtPu/Z_Pu+(1-wtPu)/Z_U);

k0 = (1./(A0+Ax*x+APu*atPu+(B0+BPu*atPu).*T)+(D./T.^2).*exp(-E./T)).*(1-P).^2.5;
k  = k_inf + (k0-k_inf)*exp(-bu/phi);

% Melting Temperature
% Validity range:
%   Deviation from stoichiometry, x: [0, 0.06] (hypo-stoichiometry). 
%   Plutonium content, [Pu]: [0, 50] at.%.
%   Burn-up, bu: [0, 110] GWd/tHM.

Tmelt_UO2 = 3147;
gamma_Pu  = 364.85;
gamma_x   = 1014.15;

Tmelt_inf = 2964.92;
delta     = 40.43;

Tmelt0 = Tmelt_UO2 - gamma_Pu*atPu - gamma_x*x;
Tmelt  = Tmelt_inf + (Tmelt0-Tmelt_inf)*exp(-bu/delta);

end

function [k,Tmelt] = Even(T,x,wtPu,P,bu)

% Thermal Conductivity
A0  = 1.81e-2;
APu = 7.92e-2;
B0  = 2.39e-4;
Bx  = 4.59e-4;
C   = 2.34e-13;
D   = 1.97e9;
E   = 1.28e4;

k_inf = 1.58;
phi   = 72.29;

k0 = (1./(A0+APu*wtPu+(B0+Bx*x).*T)+C.*T.^3+(D./T.^2).*exp(-E./T)).*(1-P).^2.5;
k  = k_inf + (k0-k_inf)*exp(-bu/phi);

% Melting Temperature
Tmelt_UO2 = 3150;
c_Pu      = 332.25;
c_x       = 995.88;

Tmelt_inf = 2912.98;
delta     = 133.27;

Tmelt0 = Tmelt_UO2 - c_Pu*wtPu - c_x*x;
Tmelt  = Tmelt_inf + (Tmelt0-Tmelt_inf)*exp(-bu/delta);

end

   %% (IV)  a)        
function [Tf_HC, Tmax] = T_NoRes(q1_max, P)
load('Results.mat', 'H_e','m_dot_ch','C_fl','H_a', 'Tfl_in', 'h','R_co','R_ci','R_fo','k_c','deltaG','N','Tmax_HC')
% Fluid - Cladding:
z = 0; % Midplane axial coordinate 
deltaT_HC = 2*H_e*q1_max/pi/m_dot_ch/C_fl*sin(pi*H_a/2/H_e);
Tfl_0_HC = Tfl_in + 0.5*deltaT_HC*(1+sin(pi*z./H_a)/sin(pi*H_a/2/H_e)); % Coolant T [°C]
Tco_0_HC = Tfl_0_HC + q1_max/h/2/pi./R_co; % Coolant radial T profile [°C]
% Oxide layer:
deltaC_0 = (0.26*(Tco_0_HC-380)+1)*10^(-6); % Thickness of the oxide layer [m]
R_ox_0 = R_co-deltaC_0; k_ox = 1; % Radial position of the cladding-oxide interface
Tox_0_HC = Tco_0_HC + q1_max*log(R_co./R_ox_0)./2./pi./k_ox; % T profile in the oxide layer [°C]
% Cladding:
Tci_0_HC = Tox_0_HC + q1_max*log(R_ox_0./R_ci)./2./pi./k_c; % Cladding radial T profile [°C]
% Gap:
theta = Tci_0_HC+273.15; % Initialization of the average gap T [K]
e = 1; % Initialization of the error
while e>1e-5
    k_g_HC    = 15.8e-4*theta^0.79; % Gap thermal conductivity [W/mK]
    Tfo_0_HC  = Tci_0_HC+q1_max*deltaG/(2*pi*R_fo*k_g_HC); % Fuel outer T [°C]
    theta_new = (Tfo_0_HC+Tci_0_HC)/2+273.15; % Update of iteration variable
    e         = abs(theta_new-theta); % Error
    theta     = theta_new;
end
Tfo_0_HC = Tci_0_HC+q1_max*(R_ci-R_fo)./(2*pi*R_fo*k_g_HC); % Gap radial T profile [°C]
% Fuel
A0 = 3.08e-2; B0 = 2.516e-4; C = 4.715e9; D = 16361;
k_fuel = @(T,P) (1./(A0+B0.*(T+273.15)) + C./(T+273.15).^2.*exp(-D./(T+273.15))).*(1-P).^2.5;
q3 = q1_max/pi/R_fo^2;

e=1;
Tmax = Tmax_HC;
while e>1e-5
    k = 1/(Tmax-Tfo_0_HC)*integral(@(T) k_fuel(T,P),Tfo_0_HC,Tmax);
    Tmax_new = Tfo_0_HC + q3*R_fo^2/4/k;
    e = abs(Tmax_new-Tmax);
    Tmax = Tmax_new;
end
r_f = linspace(0,R_fo,N);
Tf_HC = Tmax - q3.*r_f.^2/4/k;
end
