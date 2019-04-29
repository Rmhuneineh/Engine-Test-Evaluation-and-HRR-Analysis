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
Hip = 44.5; % [J/kg]
%% Set Parameters
lambda = cr/crl;
Vd = pi*bore^2/4*(2*cr); % [m^3]
Vc = Vd/(epsilon-1); % [m^3]
area = pi*bore^2/4;
x = cr*((1-cosd(theta)) + 1/lambda*(1-sqrt(1-lambda^2*sind(theta).^2))) + Vc/area; % [m] piston position
V = area*x; % [m^3] cylinder volume vs theta

%% Temperature and Elastic Constant
m_a = qm_a/60/nc/rpm; % [kg/(cycle.cylinder)] Air mass
m_f = 2*qm_f/60/nc/(rpm*2*pi/60); % [kg/(cycle.cylinder)] Fuel mass flow rate
m_egr = qm_egr*60/1000/rpm; % [kg/rad.cylinder)] EGR mass
R_egr = cp_egr*(1-1/1.4); % [J/(kg.K)] EGR elastic constant
T_cyl = p_cyl1*1e5.*V./(m_a*R_a + m_egr*R_egr) + 273; % [K] in-cylinder temperature
K = 1.338 - 6*1e-5*T_cyl + 1e-8*T_cyl.^2; % [-] Polytropic Coefficient

%% Graphically
figure(1)
hold on
plot(theta(theta>=320 & theta<=430), p_cyl1(theta>=320 & theta<=430))
plot(theta(theta>=320 & theta<=430), inj1(theta>=320 & theta<=430))
ylim([0 max(p_cyl1)])
SOI = 348.7; % [deg] Start Of Injection
SOC = 358; % [deg] Start Of Combustion
ID = SOC - SOI; % [deg] Ignition Delay

%% Pressure
m = 1.45; % [-] Compression Polytropic Coefficient
mp = 1.3; % [-] Expansion Polytropic Coficient
Vt = Vd + Vc;
% Intake (7 ->1)
P(1:1801) = min(p_cyl1);

% Compression (1 -> 2)
V12 = V(1802:3601);
P(1802:3601) =  P(1801)*(Vt./V12).^m;

% Expansion (3 -> 4)
V34 = V(3602:5402);
P(3602:5402) = P(3601)*(Vc./V34).^mp;

P(5403:7200) = P(5402);
P(theta<=SOC) = p_cyl1(theta<=SOC);
P(theta>SOC) = P(theta==SOC)*(V(theta==SOC)./V(theta>SOC)).^mp;

%% Mass Fraction Burned
for i = 2:length(theta)
       if(theta(i) <= SOI)
          Vdp(i) = 0;
          pdV(i) = 0;
       else
          Vdp(i) = trapz([p_cyl1(i-1) p_cyl1(i)], [V(i-1) V(i)]/(K(i)-1)) + Vdp(i-1);
          pdV(i) = trapz([V(i-1) V(i)], [p_cyl1(i-1) p_cyl1(i)]*K(i)/(K(i)-1)) + pdV(i-1);
       end
end
dQ_ch = Vdp + pdV;
xb = dQ_ch/m_f/Hip;


plot(theta(theta>=320 & theta<=430), P(theta>=320 & theta<=430))
yyaxis right
plot(theta(theta>=320 & theta<=430), xb(theta>=320 & theta<=430), ':', 'LineWidth', 1.5)
ylim([0 1])
plot([SOI SOI], [0 1], 'm--')
plot([SOC SOC], [0 1], 'k--')
plot([362.18 362.18], [0 1], 'g--')
plot([360.25 360.25], [0 1], 'c--')
plot([376.8 376.8], [0 1], 'c--')
legend('PCYL1', 'p_m_o_t_o_r_e_d', 'INJ1', 'x_b')
grid on