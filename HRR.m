%% Clear ALL
clear all
close all
clc

%% Load Data
[NUM, ~, ~] = xlsread('Report2_datasheet_79', 'Data');
theta = NUM(:, 1)+360; % [deg]
p_cyl1 = NUM(:, 2); % [bar]
inj1 = NUM(:, 3); % [-]

%% Engine Characteristics
cr = 52*1e-3; % [m] Crank Radius
crl = 0.158; % [m] Connecting-rod Length
bore = 95.8*1e-3; % [m] Bore
epsilon = 14.6; % [-] Compression Ratio
nc = 4; % [-] Number of Cylinders

%% Working Point
rpm = 1600; % [rpm] Engine Speed
Torque = 99.99; % [Nm]
qm_a = 97.71; % [kg/h] Air mass flow rate
qm_f = 4.072; % [kg/h] Fuel mass flow rate
qm_egr = 9.533991; % [g/s]
cp_egr = 1.038*1e3; % [J/(kg.K)]
R_a = 287.05; % [J/kg]
Hip = 44.5*1e6; % [J/kg]
%% Set Parameters
lambda = cr/crl;
Vd = pi*bore^2/4*(2*cr); % [m^3]
Vc = Vd/(epsilon-1); % [m^3]
area = pi*bore^2/4;
x = cr*((1-cosd(theta)) + 1/lambda*(1-sqrt(1-lambda^2*sind(theta).^2))) + Vc/area; % [m] piston position
V = area*x; % [m^3] cylinder volume vs theta

%% Calculation
dp(1) = (p_cyl1(1) - p_cyl1(end))*1e5; % [Pa] delta_P
dV(1) = V(1) - V(end); % [m^3] delta_V
for i = 2:length(p_cyl1)
   dp(i) = (p_cyl1(i) - p_cyl1(i-1))*1e5; % [Pa] 
   dV(i) = V(i) - V(i-1); % [m^3]
end
dp = dp';
dV = dV';
m_a = qm_a/60/(nc*rpm/2); % [kg/(cycle.cylinder)] Air mass flow rate
m_egr = qm_egr*60/1000/(nc*rpm/2); % [kg/cycle.cylinder)] Fuel mass flow rate
R_egr = cp_egr*(1-1/1.4); % [J/(kg.K)]
T_cyl = p_cyl1*1e5.*V./(m_a*R_a + m_egr*R_egr); 