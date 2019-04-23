%% Clear ALL
clear all
close all
clc

%% Load Data
EngineTesting_Results

%% Correction Factor
T0 = 298; % [K] Reference Temperature
p0_dry = 99; % [kPa] Reference Dry Pressure
p_sat_H2O = a0 + a1*T_a + a2*T_a.^2 + a3*T_a.^3 + a4*T_a.^4; % [kPa]
p_hum = Rel_Hum.*p_sat_H2O; % [kPa]
p_a_dry = p_a*100 - p_hum; % [kPa]
fa = (p0_dry./p_a_dry).^0.7.*((T_a+273)./T0).^1.5; % [-]
Vd = pi*bore^2/4*(2*cr); % [m^3]
Vc = Vd/(epsilon-1); % [m^3]
Vt = nc*(Vd + Vc)*1e3; % [L]
qf = qm_fuel*1e6/(1600/60/2)/Vt; % [mg/(L.cycle)]
qc = qf./((p_i_MF+p_a)./p_a); % [mg/(L.cycle)]
fm = 0.3*ones(4, 1);% [mg/(L.cycle)]
muC = fa.^fm; % Correction Factor
P0 = muC.*P_dyno; % [kW] Corrected Power
M0 = muC.*M_dyno; % [Nm] Corrected Torque