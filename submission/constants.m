clear all; clc; 
%% DEVELOPMENT
%% Link Lengths
crank_angle = 1:720;            % Range of crank angles (degrees)
crank_length = 46.5e-3;             % Length of the crank (m)
coupler_length = 148.8e-3;          % Length of the coupler (m)

%% Crankshaft Offsets
r0 = 0;                         % Offset between crankshaft and piston centerlines (m)
r1 = -18e-3;                      % Offset between crankshaft and piston centerlines (m)

%% Link Masses
m2 = 0.81;
m3 = 0.66;
m4 = 0.71;

%% Links Parameters
ground_angle = 0;               % Angle of the ground (degrees)
rpm = 2000;
w2 = rpm * pi/30;              % Angular velocity of the crank (ms^-1)
alpha2 = 0;                     % Initial angular acceleration (ms^-2)
OG2 = crank_length/10;
AG3 = coupler_length/3;
BG3 = 2*AG3;
J3 = (1/12)*m3*coupler_length^2;

%% Engine Parameters
rc0 = 9.8;                    % Compression ratio
rc1 = 10.1;                     % Compression ratio for offset
L = 2 * crank_length;           % Total stroke length (m)
R = coupler_length / crank_length; % Ratio of coupler to driver length
bore = 83e-3;                        % Bore radius (m)
g = 9.81;

%% General Parameters
p_atm = 101.325 * 10^3; %% pascal

%% Constants
Taf = 300; T_initial = 300; %% kelvin
yi =  1.3;
yc = 1.25;
ye = 1.48;
R = 287; %% gas constant of air (J/kg * Kelvin)
phi = 1;
CBF = 0.96; 
AFR = 14.7; 
Qhv = 43.5e6; %% J/kg
cv = R/(yc - 1); %% J/kg*k
a = 5;
m = 2;
 
%% Per Engine
d_valve = 0.4*bore;
dt = (1/(rpm*2*pi/60))*(pi/180); %% seconds
Cd = 0.7; %% check unit
At = (pi/4)*(d_valve)^2; 
Pe = 1 * p_atm; %% pascal
Pi = 1 * p_atm; %% pascal
A_cyl = (pi/4) * bore^2;

%% Valve Timings
combustion_shift = 30;
intake_start = 1;
intake_end = 222;
compression_start = 223;
compression_end = 329;
combustion_start = 360 - combustion_shift;
combustion_end = 360 + combustion_shift -1;
expansion_start = 390;
expansion_end = 487;
exhaust_start = 488;
exhaust_end = 720;

%% Convergence Variables
error = inf;
tolerence= 1e-2;
counter = 0;

flywheel_cst = 0.01;
flywheel_thickness = 50e-3;
flywheel_density = 7800;

radius_dynamic = 0.38;
rpm_max = 4000;

vel_max_car = 245;
car_mass = 2600;
friction_coef = 0.16;
incl_coef = 7/100;
incl_angle = atan2d(incl_coef, 1);
lambda = 1.3;
car_acc = 4.55;
total_efficiency = 0.9;
car_max_torque = 450;
iD = 3.636;
direct_drive = 6;
gears_total = 8;

Pf = 0.2;

nr = 2;
air_density = 1.12;

save('constants.mat')
disp('done')