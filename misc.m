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

fprintf('Engine Work per Cycle: %.4f J\n', engine_work);
fprintf('Indicated Power: %.4f W\n', power_ind);
fprintf('Mean Effective Pressure (MEP): %.4f kPa\n', MEP);
fprintf('Friction Power: %.4f W\n', power_fr);
fprintf('Mechanical Efficiency: %.4f\n', mech_efficiency);
fprintf('Brake Power: %.4f W\n', power_brake);
fprintf('Brake Mean Effective Pressure (BMEP): %.4f kPa\n', BMEP);
fprintf('Thermodynamic Torque: %.4f Nm\n', thermo_torque);
fprintf('Fuel Mass Flow Rate: %.4f kg/s\n', fuel_mf);
fprintf('Fuel Conversion Efficiency: %.4f\n', fuel_conv_efficiency);
fprintf('Specific Fuel Consumption (SFC): %.4f mg/J\n', SFC);
fprintf('Volumetric Efficiency: %.4f\n', vol_efficiency);
fprintf('Piston Speed: %.4f m/s\n', piston_speed);
fprintf('Engine Power: %.4f kW\n', Engine_power);

save('part5.mat')

