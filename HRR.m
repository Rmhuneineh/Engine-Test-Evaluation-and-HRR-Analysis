%% Clear ALL
clear all
close all
clc

%% Load Data
[NUM, ~, ~] = xlsread('Report2_datasheet_79', 'Data');
crank_angle = NUM(:, 1); % [deg]
p_cyl1 = NUM(:, 2); % [bar]
inj1 = NUM(:, 3); % [-]