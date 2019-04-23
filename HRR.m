%% Clear ALL
clear all
close all
clc

%% Load Data
[NUM, ~, ~] = xlsread('Report2_datasheet_79', 'Data');
crank_angle = NUM(:, 1); % [deg]
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
md_a = 97.71; % [kg/h] Air mass flow rate
md_f = 4.072; % [kg/h] Fuel mass flow rate
qm_egr = 9.533991; % [g/s]
cp_egr = 1.038; % [J/(g.K)]
Hip = 44.5*1e6; % [J/kg]
%% Set Parameters
theta = (0:0.1:719.9)'; % [deg]
lambda = cr/crl;
Vd = pi*bore^2/4*(2*cr); % [m^3]
Vc = Vd/(epsilon-1); % [m^3]
area = pi*bore^2/4;
x = cr*((1-cosd(theta)) + 1/lambda*(1-sqrt(1-lambda^2*sind(theta).^2))) + Vc/area; % [m] piston position
V = area*x; % [m^3] cylinder volume vs theta