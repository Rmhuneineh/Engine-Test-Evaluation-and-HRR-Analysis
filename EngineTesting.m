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
qf = qm_fuel*1e6/(n_engine/2)/nc/Vd/3600; % [mg/(L.cycle)]
qc = qf./((p_i_MF+p_a)./p_a); % [mg/(L.cycle)]
fm = 0.3*ones(4, 1);% [mg/(L.cycle)]
muC = fa.^fm; % Correction Factor
P0 = muC.*P_dyno; % [kW] Corrected Power
M0 = muC.*M_dyno; % [Nm] Corrected Torque

%% Fuel Conversion Efficiency
md_EGR = qm_EGR*1e-3; % [kg/s] EGR mass flow rate
m_EGR = md_EGR/(n_engine); % [kg] EGR mass
md_a = qm_air/3600; % [kg/s] Air mass flow rate
m_a = md_a/(n_engine); % [kg] Air mass
md_f = qm_fuel/3600; % [kg/s] Fuel mass flow rate
m_f = md_f/(n_engine); % [kg] Fuel mass
eta_f = P0*1e3./(md_f*Hip); % [-] Fuel Conversion Efficiency

%% Brake Specific Fuel Consumption (bsfc)
bsfc = qm_fuel*1e3./P0; % [g/kW.h]

%% Volumetric Efficiency
ro = (p_i_MF+p_a)*1e5./(R_a*(T_a+273)); % [kg/m^3] Density
lambdav = (md_a+md_EGR)./(ro*nc*Vd*n_engine/2); % [-]

