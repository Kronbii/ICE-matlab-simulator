clear all; clc; % Clear command window and workspace

%% Load Necessary Data
load('constants.mat')

%% Declaring iterating variables
s_length0 = zeros(1, length(crank_angle)); % Preallocate s_length for efficiency
s_length1 = zeros(1, length(crank_angle)); % Preallocate s_length for efficiency
coupler_angle0 = zeros(1, length(crank_angle)); % Preallocate for coupler angles
coupler_angle1 = zeros(1, length(crank_angle)); % Preallocate for coupler angles
vel_piston0 = zeros(1, length(crank_angle)); % Preallocate for velocities
vel_piston1 = zeros(1, length(crank_angle)); % Preallocate for velocities
w3 = zeros(1, length(crank_angle)); % Preallocate for velocities
alpha3 = zeros(1, length(crank_angle));
a2x = zeros(1, length(crank_angle));
a2y = zeros(1, length(crank_angle));
a3x = zeros(1, length(crank_angle));
a3y = zeros(1, length(crank_angle));
a4y = zeros(1, length(crank_angle));

%% Numerical Analysis
for i = 1:length(crank_angle)
    %%% Position Analysis %%%
    % Solve for s_length (length on follower track)
    tempa = 1;
    tempb = -2 * crank_length * cosd(crank_angle(i));
    tempc0 = crank_length^2 + r0^2 - coupler_length^2 + 2 * r0 * crank_length * sind(crank_angle(i));
    tempc1 = crank_length^2 + r1^2 - coupler_length^2 + 2 * r1 * crank_length * sind(crank_angle(i));

    tempbdot = 2 * crank_length * w2 * sind(crank_angle(i));
    tempc0dot = 0; % No derivative for tempc0
    tempc1dot = 2 * r1 * crank_length * w2 * cosd(crank_angle(i));

    % Calculate s_length for both configurations
    s_length0(i) = (-tempb + sqrt(tempb^2 - 4 * tempa * tempc0)) / (2 * tempa);
    s_length1(i) = (-tempb + sqrt(tempb^2 - 4 * tempa * tempc1)) / (2 * tempa);

    coupler_angle0(i) = atan2d((-r0-crank_length*sind(crank_angle(i))),(-crank_length*cosd(crank_angle(i))+s_length0(i)));
    coupler_angle1(i) = atan2d((-r1-crank_length*sind(crank_angle(i))),(-crank_length*cosd(crank_angle(i))+s_length1(i)));
                                
    % Calculate piston velocities for both configurations
    vel_piston0(i) = 0.5 * (-tempbdot + (0.5 / sqrt(tempb^2 - 4 * tempa * tempc0)) * (2 * tempb * tempbdot - 4 * tempc0dot));
    vel_piston1(i) = 0.5 * (-tempbdot + (0.5 / sqrt(tempb^2 - 4 * tempa * tempc1)) * (2 * tempb * tempbdot - 4 * tempc1dot));
    
    % solve for w3                          
    w3(i) = -(crank_length*w2*cosd(crank_angle(i)))/(coupler_length*cosd(coupler_angle0(i)));
    
    %solve for alpha3
    % Calculate alpha3 using the formula for the current step
    alpha3(i) = (crank_length * alpha2 * cosd(crank_angle(i)) - ...
                        crank_length * w2^2 * sind(crank_angle(i)) + ...
                        coupler_length * w3(i)^2 * sind(coupler_angle0(i))) / ...
                        (coupler_length * cosd(coupler_angle0(i)));
    
    % Calculate the x and y components of the acceleration a2 for the current step
    a2x(i) =  w2^2*OG2*sind(crank_angle(i));
    a2y(i) =  -w2^2*OG2*cosd(crank_angle(i));
    
    % Calculate the x and y components of acceleration a3 for the current step
    a3x(i) = w2^2 * crank_length * sind(crank_angle(i)) + ...
                    w3(i)^2 * AG3 * sind(coupler_angle0(i)) - ...
                    alpha3(i) * AG3 * cosd(coupler_angle0(i));

    a3y(i) = -w2^2 * crank_length * cosd(crank_angle(i)) - ...
                    w3(i)^2 * AG3 * cosd(coupler_angle0(i)) - ...
                    alpha3(i) * AG3 * sind(coupler_angle0(i));
    
                % Define the equation for d_double_dot
    a4y(i) = -crank_length * alpha2 * sind(crank_angle(i)) - ...
                   crank_length * w2^2 * cosd(crank_angle(i)) + ...
                   coupler_length * alpha3(i) * sind(coupler_angle0(i)) + ...
                   coupler_length * w3(i)^2 * cosd(coupler_angle0(i));
               

end

%% Displacement Calculation
c_length0 = max(s_length0) - s_length0;
c_length1 = max(s_length1) - s_length1;

vdisp0 = ((pi*bore^2)/4)*max(s_length0);
vdisp1 = ((pi*bore^2)/4)*max(s_length1);

%% Displaced Volume Calculation
total_volume0 = vdisp0/(rc0 -1) + ((pi*bore^2)/4)*c_length0;
total_volume1 = vdisp1/(rc1 -1) + ((pi*bore^2)/4)*c_length1;

%% Plots
% Plot Velocities vs Driver Angle
figure('Name', 'Velocity vs Driver Angle');
hold on; % Hold on to add more elements to the plot
grid on;
% Plot Title
title('Piston Velocity vs Driver Angle', 'FontWeight', 'bold', 'FontSize', 14);
% Plot the velocities
plot(crank_angle, vel_piston0, 'LineWidth', 2, 'DisplayName', 'No Offset');
plot(crank_angle, vel_piston1, 'LineWidth', 2, 'DisplayName', 'Offseted');
% Axis labels
xlabel('Crank Angle (degrees)', 'FontWeight', 'bold');
ylabel('Velocity (cm/s)', 'FontWeight', 'bold');
% Set limits
xlim([0 length(crank_angle)]); % Assuming crank_angle ranges from 0 to 360 degrees
% Legend
legend('Location', 'best', 'FontSize', 10);
xline(90, 'k--', 'DisplayName', '90 degrees')
xline(180, 'k--', 'DisplayName', '180 degrees')
xline(270, 'k--', 'DisplayName', '270 degrees')
% Set axis properties
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
hold off;
 
 
 
% Plot Displacement vs Driver Angle
figure('Name', 'Displacement vs Driver Angle');
hold on; % Hold on to add more elements to the plot
grid on;
% Plot Title
title('Piston Displacement vs Driver Angle', 'FontWeight', 'bold', 'FontSize', 14);
% Plot the velocities
plot(crank_angle, 100*c_length0, 'LineWidth', 2, 'DisplayName', 'No Offset');
plot(crank_angle, 100*c_length1, 'LineWidth', 2, 'DisplayName', 'Offseted');
% Axis labels
xlabel('Crank Angle (degrees)', 'FontWeight', 'bold');
ylabel('Displacement (cm)', 'FontWeight', 'bold');
% Set limits
xlim([0 length(crank_angle)]); % Assuming crank_angle ranges from 0 to 360 degrees
% Legend
legend('Location', 'best', 'FontSize', 10);
xline(90, 'k--', 'DisplayName', '90 degrees')
xline(180, 'k--', 'DisplayName', '180 degrees')
xline(270, 'k--', 'DisplayName', '270 degrees')
% Set axis properties
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
hold off;
 
 
 
% Plot Cylinder Volume vs Driver Angle
figure('Name', 'Piston Volume vs Driver Angle');
hold on; % Hold on to add more elements to the plot
grid on;
% Plot Title
title('Cylinder Volume vs Driver Angle', 'FontWeight', 'bold', 'FontSize', 14);
% Plot the velocities
plot(crank_angle, total_volume0*10^6, 'LineWidth', 2, 'DisplayName', 'No Offset');
plot(crank_angle, total_volume1*10^6, 'LineWidth', 2, 'DisplayName', 'Offseted');
% Axis labels
xlabel('Crank Angle (degrees)', 'FontWeight', 'bold');
ylabel('Piston Volume (cc)', 'FontWeight', 'bold');
% Set limits
xlim([0 length(crank_angle)]); % Assuming crank_angle ranges from 0 to 360 degrees
% Legend
legend('Location', 'best', 'FontSize', 10);
xline(90, 'k--', 'DisplayName', '90 degrees')
xline(180, 'k--', 'DisplayName', '180 degrees')
xline(270, 'k--', 'DisplayName', '270 degrees')
% Set axis properties
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
hold off;

%% Data save
save('part1.mat')
            
disp('done')