clear all; clc;

load('part4.mat')
Pf = 0.2;
%%%%%%%%%%%%%%%%%%%%%%% Work/Cycle %%%%%%%%%%%%%%%%%%%%%%%%%
% Define your data and calculations
engine_work = trapz(vol, pressure); % work/cycle
power_ind = ((engine_work*(rpm/60))/nr); % indicated power in W
MEP = (power_ind*nr)/(vdisp0*1000*rpm/60); % mean effective pressure in KPa
power_fr = Pf*power_ind; % Friction power
mech_efficiency = (1-(power_fr/power_ind)); % mechanical efficiency
power_brake = mech_efficiency*power_ind; % Brake power in W
BMEP = (power_brake*nr)/(vdisp0*1000*rpm/60); % brake mean effective pressure in KPa
thermo_torque = (BMEP*vdisp0*1000)/(6.28*nr); % thermo torque in Nm
fuel_mf = (mf*(rpm/60))/nr; % fuel mass flow rate in kg/s
fuel_conv_efficiency = (engine_work)/((mf)*(Qhv)); % Fuel conversion efficiency
SFC = (fuel_mf*10^(-3))/(power_ind*10^(-3)); % specific fuel consumption in mg/J
vol_efficiency = (max(mass)/(air_density*vdisp0)); % Volumetric efficiency
piston_speed = 2*L*rpm/60; % piston speed in m/sec
Engine_power = (power_ind*4)*(10^-3); % kW

% Create a UI figure window
fig = uifigure('Name', 'Engine Performance Parameters', 'Position', [100, 100, 500, 600]);

% Add labels to display the calculated values
lbl_engine_work = uilabel(fig, 'Position', [20, 550, 460, 22], 'Text', ['Engine Work: ', num2str(engine_work), ' J/cycle']);
lbl_power_ind = uilabel(fig, 'Position', [20, 510, 460, 22], 'Text', ['Indicated Power: ', num2str(power_ind), ' W']);
lbl_MEP = uilabel(fig, 'Position', [20, 470, 460, 22], 'Text', ['Mean Effective Pressure (MEP): ', num2str(MEP), ' KPa']);
lbl_power_fr = uilabel(fig, 'Position', [20, 430, 460, 22], 'Text', ['Friction Power: ', num2str(power_fr), ' W']);
lbl_mech_efficiency = uilabel(fig, 'Position', [20, 390, 460, 22], 'Text', ['Mechanical Efficiency: ', num2str(mech_efficiency)]);
lbl_power_brake = uilabel(fig, 'Position', [20, 350, 460, 22], 'Text', ['Brake Power: ', num2str(power_brake), ' W']);
lbl_BMEP = uilabel(fig, 'Position', [20, 310, 460, 22], 'Text', ['Brake Mean Effective Pressure (BMEP): ', num2str(BMEP), ' KPa']);
lbl_thermo_torque = uilabel(fig, 'Position', [20, 270, 460, 22], 'Text', ['Thermo Torque: ', num2str(thermo_torque), ' Nm']);
lbl_fuel_mf = uilabel(fig, 'Position', [20, 230, 460, 22], 'Text', ['Fuel Mass Flow Rate: ', num2str(fuel_mf), ' kg/s']);
lbl_fuel_conv_eff = uilabel(fig, 'Position', [20, 190, 460, 22], 'Text', ['Fuel Conversion Efficiency: ', num2str(fuel_conv_efficiency)]);
lbl_SFC = uilabel(fig, 'Position', [20, 150, 460, 22], 'Text', ['Specific Fuel Consumption (SFC): ', num2str(SFC), ' mg/J']);
lbl_vol_efficiency = uilabel(fig, 'Position', [20, 110, 460, 22], 'Text', ['Volumetric Efficiency: ', num2str(vol_efficiency)]);
lbl_piston_speed = uilabel(fig, 'Position', [20, 70, 460, 22], 'Text', ['Piston Speed: ', num2str(piston_speed), ' m/s']);
lbl_Engine_power = uilabel(fig, 'Position', [20, 30, 460, 22], 'Text', ['Engine Power: ', num2str(Engine_power), ' kW']);

