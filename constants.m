clear all; clc;
%% Engine and vehicle parameters, saved to constants.mat for part1.m and part2.m
% Lengths in m, angles in degrees, pressures in Pa, temperatures in K.

%% Crank angle grid
crank_angle = 1:720;                % one four-stroke cycle, 1 degree steps (720 = 0 = TDC)

%% Slider-crank geometry
crank_length = 46.5e-3;             % crank radius (m)
coupler_length = 148.8e-3;          % connecting-rod length (m)
r0 = 0;                             % cylinder axis through the crankshaft (m)
r1 = -18e-3;                        % offset cylinder axis, for comparison (m)
bore = 83e-3;                       % bore diameter (m)
n_cyl = 4;                          % cylinders, firing every 180 degrees (1-3-4-2)

%% Link masses and inertia
m2 = 0.81;                          % crank (kg)
m3 = 0.66;                          % connecting rod (kg)
m4 = 0.71;                          % piston assembly (kg)
OG2 = crank_length/10;              % crank centre of mass from the main journal (m)
AG3 = coupler_length/3;             % rod centre of mass from the crank pin (m)
BG3 = 2*AG3;                        % rod centre of mass from the piston pin (m)
J3 = (1/12)*m3*coupler_length^2;    % rod moment of inertia about its centre of mass (kg m^2)
g = 9.81;

%% Operating point
rpm = 2000;                         % design speed for the force analysis and flywheel
w2 = rpm*pi/30;                     % crank angular velocity (rad/s)
alpha2 = 0;                         % constant crank speed
rpm_sweep = 1000:500:6000;          % speeds for the full-load torque and power curves

%% Compression ratios
rc0 = 9.8;                          % cylinder axis through the crankshaft
rc1 = 10.1;                         % offset cylinder axis

%% Gas properties
p_atm = 101325;                     % Pa
R = 287;                            % gas constant (J/kg K)
yi = 1.3;                           % ratio of specific heats, fresh mixture
yc = 1.25;                          % ratio of specific heats, burned gas
Taf = 300;                          % intake temperature (K)
Pi = p_atm;                         % intake manifold pressure (Pa)
Pe = p_atm;                         % exhaust manifold pressure (Pa)

%% Fuel and combustion
phi = 1;                            % equivalence ratio
AFR = 14.7;                         % stoichiometric air-fuel ratio
Qhv = 43.5e6;                       % fuel heating value (J/kg)
CBF = 0.96;                         % combustion efficiency
a = 5;                              % Wiebe efficiency factor
m = 2;                              % Wiebe form factor

%% Valve and combustion timing, in crank degrees from TDC at the start of intake
intake_start = 1;                   % intake valve opens (TDC)
intake_end = 222;                   % intake valve closes (42 degrees after BDC)
compression_start = 223;
compression_end = 329;
combustion_shift = 30;              % spark advance (degrees before TDC)
combustion_start = 360 - combustion_shift;      % spark, start of the burn
combustion_end = 360 + combustion_shift - 1;    % 60 degree burn duration
expansion_start = 390;
expansion_end = 487;
exhaust_start = 488;                % exhaust valve opens (52 degrees before BDC)
exhaust_end = 720;                  % exhaust valve closes (TDC)

%% Valves
d_valve = 0.4*bore;                 % valve head diameter (m)
lift_max = 0.25*d_valve;            % peak lift: the curtain area reaches the port area
Cd = 0.7;                           % discharge coefficient

%% Heat transfer to the walls (Woschni)
T_wall = 400;                       % mean wall temperature (K)

%% Friction: Heywood's correlation for four-cylinder SI engines,
% fmep (bar) = 0.97 + 0.15 (N/1000) + 0.05 (N/1000)^2
fmep_coef = [0.97 0.15 0.05];

%% Cycle convergence
tolerance = 1e-5;                   % relative change of the state from one cycle to the next
max_cycles = 40;

%% Flywheel
flywheel_cst = 0.04;                % coefficient of speed fluctuation
flywheel_thickness = 50e-3;         % steel disc (m)
flywheel_density = 7800;            % kg/m^3

%% Vehicle and transmission (the vehicle the gearbox is designed for)
radius_dynamic = 0.38;              % tyre dynamic radius (m)
rpm_max = 4000;                     % engine speed at top speed
vel_max_car = 245;                  % top speed (km/h)
car_mass = 2600;                    % kg
friction_coef = 0.16;               % rolling resistance coefficient
incl_coef = 7/100;                  % 7% grade
incl_angle = atan(incl_coef)*180/pi;
lambda = 1.3;                       % rotating mass factor
car_acc = 4.55;                     % acceleration demanded in first gear (m/s^2)
total_efficiency = 0.9;             % driveline efficiency
car_max_torque = 450;               % engine torque the gearbox is designed for (N m)
iD = 3.636;                         % final drive ratio
direct_drive = 6;                   % gear with ratio 1
gears_total = 8;

save('constants.mat', '-v7')
disp('done')
