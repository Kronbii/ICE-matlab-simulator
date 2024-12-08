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
% T_avg1 = 14250;
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

% integ_range = 20:110;
% 
% fun = @(x) 2200*sind(2*x) - 1800*cosd(2*x);
% energy = integral(fun, 20, 110);
% disp(['Integral: ', num2str(energy)]);

integ_range = min(theta_int):max(theta_int);

energy = trapz(integ_range,T0(integ_range) - T_avg);
disp(['Energy: ', num2str(energy)]);

%% Calculating Interia Using Energy
Is = energy/(flywheel_cst * w2^2);

%% Calculating The Flywheel Radius
flywheel_radius = ((2 * Is) / (pi * flywheel_thickness * flywheel_density))^(1/4);
disp(['Radius: ', num2str(flywheel_radius)]);
disp(['Average Torque: ', num2str(T_avg)]);

% Create a new figure window for the plot
figure;

% Plot the torque before adding the flywheel (T0) in red color
plot(crank_angle, T0, 'r', 'LineWidth', 1);  % Red line for torque before flywheel
hold on;  % Keep the current plot and add the next one

% Plot the average torque (T_avg) in blue with circle markers
plot(crank_angle, T_avg, '-ob', 'MarkerSize', 1, 'LineWidth', 1);  % Blue line with circle markers

% Add a title to the plot to describe what is being displayed
title('Torque Before and After Adding Flywheel', 'FontSize', 14);

% Label the x-axis to indicate crank angle in degrees
xlabel('Crank Angle (degrees)', 'FontSize', 12);

% Label the y-axis to indicate torque in Newton-meters (Nm)
ylabel('Torque (Nm)', 'FontSize', 12);

% Add grid lines to the plot for better readability
grid on;

% Add a legend to differentiate between the torque before adding the flywheel and the average torque
legend('Torque Before Flywheel', 'Torque With Flywheel', 'Location', 'Best');

% Set axis limits to ensure both curves are fully visible
axis([min(crank_angle) max(crank_angle) min([T0, T_avg]) max([T0, T_avg])]);

% Optionally add a horizontal line at y=0 to visualize the baseline (zero torque level)
yline(0, 'k--');  % Horizontal dashed line at y=0

% Display the plot
hold off;  % Release the current plot hold

save('part4.mat')
