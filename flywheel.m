%set K =  0.01
% Is = E/K*wavg^2
% E = integral
% we need w given
% circular steel disk
% mass density = 7800 kg/m3 and 50 mm thick

%% Steps
% get mean torque
% get angles for which torque equal to avg torque
% get E from integrating T - Tavg from angle min to angle max
% calculate Is = E/K*wavg^2
% get radius from derived equation ( radius = 2*Is/(pi * thickness * density)

clear all; clc;

%% Load Necessary Data
load('part3.mat')

%% Calculating Mean Torque
T_avg = mean(T0); %per cylinder
T_avg_total = mean(T_total); %total for engine

%% Calculating the crank angles at which the avg and instantaneous torque intersect
theta_int = find(abs(diff(sign(T0 - T_avg))) > 0);
theta_int = theta_int(theta_int >= combustion_start & theta_int <= expansion_end);

%% Calculating the Energy using trapz
%energy = trapz(crank_angle(min(theta_int):max(theta_int)), T0(min(theta_int):max(theta_int)) - T_avg);  % Integrate over theta_range and T_diff_range
T_avg1 = 14250;
% 
% % Loop through crank angles with a step of 0.1
% for i = 1:length(crank_angle)
%     Titi(i) = 2200 * sind(2 * i) - 1800 * cosd(2 * i);
% end
% 
% integral_Torque = trapz(20:110, Titi(20:110));
% 
% % Calculate E in J
% energy = (integral_Torque);
% 
% disp(['Energy: ', num2str(energy)]);
% plot(crank_angle(20:110), Titi(20:110))

integ_range = 20:110;

fun = @(x) 2200*sind(2*x) - 1800*cosd(2*x);
energy = integral(fun, 20, 110);
disp(['Integral: ', num2str(energy)]);



%% Calculating Interia Using Energy
w21 = 150 * pi/30;
Is = energy/(flywheel_cst * w21^2);

%% Calculating The Flywheel Radius
flywheel_radius = ((2 * Is) / (pi * flywheel_thickness * flywheel_density))^(1/4);
disp(['Radius: ', num2str(flywheel_radius)]);

save('part4.mat')
