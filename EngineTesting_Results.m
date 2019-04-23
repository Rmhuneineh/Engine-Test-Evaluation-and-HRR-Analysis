%% Engine
n_engine = 1600*ones(4, 1)*2*pi/60; % [rad/s] Engine Speed
P_dyno = [3.35; 8.38; 16.75; 21.78]; % [kW] Measured Power
M_dyno = [20; 50.02; 99.99; 129.99]; % [Nm] Measured Torque
qm_fuel = [1.723; 2.485; 4.072; 5.063]/3600; % [kg/s]
qm_air = [134.75; 81.96; 97.71; 111.36]/3600; % [kg/s]
qm_EGR = [2.781567151; 8.476373487; 9.539912943; 10.30479605]; % [g/s]

%% Ambient
p_a = [1000.12; 1000.08; 1000.03; 1000]*1e-3; % [bar] Barometric Pressure
T_a = [25.2; 25.3; 25.3; 25.5]; % [C] Temperature
Rel_Hum = [41.35; 40.3; 42.24; 48.35]/100; % [%] Relative Humidity

%% Intake Manifold
T_i_MF = [27.82; 64.21; 77.7; 78.14]; % [C]
p_i_MF = [131.1; 82.3; 200.9; 308.5]*1e-3; % [bar] Relative Manifold Pressure

%% Water Saturation Coefficient
a0 = 0.611511;
a1 = 0.0437968;
a2 = 0.00155832;
a3 = 1.7572*1e-5;
a4 = 5.4683*1e-7;

%% Fuel Characteristics
Hip = 44.2*1e6; % [J/kg]
R_f = 417; % [J/kg]

%% Air
R_a = 287.05; % [J/kg]

%% Engine Characteristics
cr = 52*1e-3; % [m]
crl = 0.158; % [m]
bore = 95.8*1e-3; % [m]
epsilon = 14.6; % [-]
nc = 4; % [-]