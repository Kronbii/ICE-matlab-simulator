clc; clear; % Clear command window and workspace

%% Variables and Declarations
%%%%% Link Lengths %%%%
% Inputs
crank_angle = 1:720;            % Range of crank angles (degrees)
crank_length = 3.75;             % Length of the crank (cm) %%%%%changed
coupler_length = 11.25;          % Length of the coupler (cm) %%%%%changed
r0 = 0;                         % Offset between crankshaft and piston centerlines (cm)
r1 = -1.8;                      % Offset between crankshaft and piston centerlines (cm)
ground_angle = 0;               % Angle of the ground (degrees)
w2 = 2200 * pi/30;              % Angular velocity of the crank (ms^-1) %%%%%changed
alpha2 = 0;                     % Initial angular acceleration (ms^-2)
L = 2 * crank_length;           % Total stroke length (cm)
R = coupler_length / crank_length; % Ratio of coupler to driver length
rc0 = 10;                    % Compression ratio %%%%%changed
rc1 = 10.1;                     % Compression ratio for offset
B = 8.4;                        % Bore radius %%%%%changed
vdisp0 = 415.67;                 % total displaced volume in (cc)
vdisp1 = 697.9;                 % total displaced volume for offset in (cc)
 
s_length0 = zeros(1, length(crank_angle)); % Preallocate s_length for efficiency
s_length1 = zeros(1, length(crank_angle)); % Preallocate s_length for efficiency
coupler_angle0 = zeros(1, length(crank_angle)); % Preallocate for coupler angles
coupler_angle1 = zeros(1, length(crank_angle)); % Preallocate for coupler angles
vel_piston0 = zeros(1, length(crank_angle)); % Preallocate for velocities
vel_piston1 = zeros(1, length(crank_angle)); % Preallocate for velocities
 
%% Numerical Analysis
for i = 1:length(crank_angle)
    %%% Position Analysis %%%
    % Solve for s_length (length on follower track)
    tempa = 1;
    tempb = -2*crank_length*cosd(crank_angle(i));
    tempc0 = crank_length^2 + r0^2 - coupler_length^2 +2*r0*crank_length*sind(crank_angle(i));
    tempc1 = crank_length^2 + r1^2 - coupler_length^2 +2*r1*crank_length*sind(crank_angle(i));

    tempbdot = 2*crank_length*w2*sind(crank_angle(i));
    tempc0dot = 0;
    tempc1dot = 2*r1*crank_length*w2*cosd(crank_angle(i));

    s_length0(i) = (-tempb + sqrt(tempb^2 - 4*tempa*tempc0))/(2*tempa);
    s_length1(i) = (-tempb + sqrt(tempb^2 - 4*tempa*tempc1))/(2*tempa);


    % Calculate coupler angle using atan2
    coupler_angle0(i) = atan2d(sqrt((s_length0(i) * sind(ground_angle) - r0 - crank_length * sind(crank_angle(i)))^2 / coupler_length^2), ...
                               sqrt((s_length0(i) * cosd(ground_angle) - crank_length * cosd(crank_angle(i)))^2 / coupler_length^2));
    
 
    % Calculate coupler angle using atan2
    coupler_angle1(i) = atan2d(sqrt((s_length1(i) * sind(ground_angle) - r1 - crank_length * sind(crank_angle(i)))^2 / coupler_length^2), ...
                               sqrt((s_length1(i) * cosd(ground_angle) - crank_length * cosd(crank_angle(i)))^2 / coupler_length^2));
                           
    %%% Velocity Analysis %%%
    vel_piston0(i) = 0.5*(-tempbdot + (0.5/(sqrt(tempb^2 - 4*tempa*tempc0)))*(2*tempb*tempbdot - 4*tempc0dot));
    vel_piston1(i) = 0.5*(-tempbdot + (0.5/(sqrt(tempb^2 - 4*tempa*tempc1)))*(2*tempb*tempbdot - 4*tempc1dot));
end

%% Plotting Live Curves
%{
for i=1:length(crank_angle)
    % Joint positions
    o1_x = 0;
    o1_y = 0;
    o2_x = 0; o2_y = r0; % Origin of the coupler
    A_x = o2_x + crank_length * cosd(crank_angle(i)); % Position of point A
    A_y = o2_y + crank_length * sind(crank_angle(i));
    o4_x = s_length0(i) * cosd(ground_angle); % Position of point O4
    o4_y = s_length0(i) * sind(ground_angle);
    
    plot([o1_x o4_x],[o1_y,o4_y],[o1_x o2_x],[o1_y,o2_y],[o2_x A_x],[o2_y,A_y],[A_x, o4_x],[A_y,o4_y],'LineWidth',3)
    grid on;
    rectangle('Position',[o4_x-0.35 o4_y-0.35 0.5 0.5],'FaceColor',[0 1 1],'LineWidth',3)
    axis([-10 30 -10 20])
    pause(0.01);
end

for i=1:length(crank_angle)
    o1_x = 0;
    o1_y = 0;
    o2_x = 0; o2_y = r1; % Origin of the coupler
    A_x = o2_x + crank_length * cosd(crank_angle(i)); % Position of point A
    A_y = o2_y + crank_length * sind(crank_angle(i));
    o4_x = s_length1(i) * cosd(ground_angle); % Position of point O4
    o4_y = s_length1(i) * sind(ground_angle);
    
    plot([o1_x o4_x],[o1_y,o4_y],[o1_x o2_x],[o1_y,o2_y],[o2_x A_x],[o2_y,A_y],[A_x, o4_x],[A_y,o4_y],'LineWidth',3)
    grid on;
    rectangle('Position',[o4_x-0.35 o4_y-0.35 0.5 0.5],'FaceColor',[0 1 1],'LineWidth',3)
    axis([-10 30 -10 20])
    pause(0.01);
end
%}

%% Volume Calculations
%%% Displacement Calculation %%%
c_length0 = max(s_length0) - s_length0;
%%% Displaced Volume Calculation %%%
total_volume0 = vdisp0/(rc0 -1) + ((pi*B^2)/4)*c_length0;
 
%%% Displacement Calculation %%%
c_length1 = max(s_length1) - s_length1;
%%% Displaced Volume Calculation %%%
total_volume1 = vdisp1/(rc1 -1) + ((pi*B^2)/4)*c_length1;
 

%% Plotting Results
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
%calculate max values
[piston0_max, piston0_index] = max(vel_piston0);
[piston1_max, piston1_index] = max(vel_piston1);
plot(crank_angle(piston0_index), piston0_max, 'kx', 'MarkerSize', 8, 'DisplayName', sprintf('Max baseline (%.1f, %.2f)', crank_angle(piston0_index), piston0_max)); 
plot(crank_angle(piston1_index), piston1_max, 'kx', 'MarkerSize', 8, 'DisplayName', sprintf('Max offset (%.1f, %.2f)', crank_angle(piston1_index), piston1_max));
%calculate minimum values
[piston0_min, piston0_index] = min(vel_piston0);
[piston1_min, piston1_index] = min(vel_piston1);
plot(crank_angle(piston0_index), piston0_min, 'kx', 'MarkerSize', 8, 'DisplayName', sprintf('Min baseline (%.1f, %.2f)', crank_angle(piston0_index), piston0_min)); 
plot(crank_angle(piston1_index), piston1_min, 'kx', 'MarkerSize', 8, 'DisplayName', sprintf('Min offset (%.1f, %.2f)', crank_angle(piston1_index), piston1_min));
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
plot(crank_angle, c_length0, 'LineWidth', 2, 'DisplayName', 'No Offset');
plot(crank_angle, c_length1, 'LineWidth', 2, 'DisplayName', 'Offseted');
% Axis labels
xlabel('Crank Angle (degrees)', 'FontWeight', 'bold');
ylabel('Displacement (cm)', 'FontWeight', 'bold');
% Set limits
xlim([0 length(crank_angle)]); % Assuming crank_angle ranges from 0 to 360 degrees
% Legend
legend('Location', 'best', 'FontSize', 10);
%calculate max values
[piston0_max, piston0_index] = max(c_length0);
[piston1_max, piston1_index] = max(c_length1);
plot(crank_angle(piston0_index), piston0_max, 'kx', 'MarkerSize', 8, 'DisplayName', sprintf('Max baseline (%.1f, %.2f)', crank_angle(piston0_index), piston0_max)); 
plot(crank_angle(piston1_index), piston1_max, 'kx', 'MarkerSize', 8, 'DisplayName', sprintf('Max offset(%.1f, %.2f)', crank_angle(piston1_index), piston1_max));
%calculate minimum values
[piston0_min, piston0_index] = min(c_length0);
[piston1_min, piston1_index] = min(c_length1);
plot(crank_angle(piston0_index), piston0_min, 'kx', 'MarkerSize', 8, 'DisplayName', sprintf('Min baseline (%.1f, %.2f)', crank_angle(piston0_index), piston0_min)); 
plot(crank_angle(piston1_index), piston1_min, 'kx', 'MarkerSize', 8, 'DisplayName', sprintf('Min offset (%.1f, %.2f)', crank_angle(piston1_index), piston1_min));
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
plot(crank_angle, total_volume0, 'LineWidth', 2, 'DisplayName', 'No Offset');
plot(crank_angle, total_volume1, 'LineWidth', 2, 'DisplayName', 'Offseted');
% Axis labels
xlabel('Crank Angle (degrees)', 'FontWeight', 'bold');
ylabel('Piston Volume (cc)', 'FontWeight', 'bold');
% Set limits
xlim([0 length(crank_angle)]); % Assuming crank_angle ranges from 0 to 360 degrees
% Legend
legend('Location', 'best', 'FontSize', 10);
%calculate max values 
[volume0_max, volume0_index] = max(abs(total_volume0));
[volume1_max, volume1_index] = max(abs(total_volume1));
xline(crank_angle(volume0_index), 'b--', 'LineWidth', 1.5, 'DisplayName', sprintf(' (%.2f, %.2f)', volume0_index, volume0_max));
xline(crank_angle(volume1_index), 'r--', 'LineWidth', 1.5, 'DisplayName', sprintf(' (%.2f, %.2f)', volume1_index, volume1_max));
xline(crank_angle(piston1_index), 'r--', 'LineWidth', 1.5, 'DisplayName', sprintf(' (%.2f, %.2f)', piston1_index, piston1_max));
xline(crank_angle(90), 'k--', 'LineWidth', 1.5, 'DisplayName', "Crank Angle: 90");
xline(crank_angle(180), 'k--', 'LineWidth', 1.5, 'DisplayName', "Crank Angle: 180");
xline(crank_angle(270), 'k--', 'LineWidth', 1.5, 'DisplayName', "Crank Angle: 270");
% Set axis properties
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
hold off;

save('data.mat' , 'total_volume0', 'total_volume1', 'crank_angle');