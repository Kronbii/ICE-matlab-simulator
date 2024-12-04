clc; clear; % Clear command window and workspace
 
%% Link Lengths
crank_angle = 1:720;            % Range of crank angles (degrees)
crank_length = 4.65;             % Length of the crank (cm)
coupler_length = 14.88;          % Length of the coupler (cm)

%% Crankshaft Offsets
r0 = 0;                         % Offset between crankshaft and piston centerlines (cm)
r1 = -1.8;                      % Offset between crankshaft and piston centerlines (cm)

%% Link Masses
m2 = 0.81;
m3 = 0.66;
m4 = 0.71;

%% Links Parameters
ground_angle = 0;               % Angle of the ground (degrees)
rpm = 1000;
w2 = rpm * pi/30;              % Angular velocity of the crank (ms^-1)
alpha2 = 0;                     % Initial angular acceleration (ms^-2)
OG2 = crank_length/10;
AG3 = coupler_length/3;
BG3 = 2*AG3;
J3 = (1/12)*m3*coupler_length^2;

%% Engine Parameters
rc0 = 9.8;                    % Compression ratio
rc1 = 10.1;                     % Compression ratio for offset
L = 2 * crank_length;           % Total stroke length (cm)
R = coupler_length / crank_length; % Ratio of coupler to driver length
B = 8.3;                        % Bore radius
g = 9.81;

%% General Parameters
p_atm = 101.325 * 10^3; %% pascal

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

    % Calculate coupler angles using atan2
    coupler_angle0(i) = atan2d(sqrt((s_length0(i) * sind(ground_angle) - r0 - crank_length * sind(crank_angle(i)))^2 / coupler_length^2), ...
                                   sqrt((s_length0(i) * cosd(ground_angle) - crank_length * cosd(crank_angle(i)))^2 / coupler_length^2));

    coupler_angle1(i) = atan2d(sqrt((s_length1(i) * sind(ground_angle) - r1 - crank_length * sind(crank_angle(i)))^2 / coupler_length^2), ...
                                   sqrt((s_length1(i) * cosd(ground_angle) - crank_length * cosd(crank_angle(i)))^2 / coupler_length^2));
                               
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

vdisp0 = ((pi*B^2)/4)*max(s_length0);
vdisp1 = ((pi*B^2)/4)*max(s_length1);

%% Displaced Volume Calculation
total_volume0 = vdisp0/(rc0 -1) + ((pi*B^2)/4)*c_length0;
total_volume1 = vdisp1/(rc1 -1) + ((pi*B^2)/4)*c_length1;

save('R:\University\Courses\Fall 24\Automotive\assignment-3\codes\part1.mat','a2x', 'a2y', 'a3x', 'a3y', 'w3', 'w2', 'alpha2',...
                'alpha3', 'g', 'm2', 'm3', 'm4', 'OG2', 'AG3', 'BG3', 'J3', 'coupler_angle0', ...
                'crank_angle', 'crank_length', 'coupler_length', 'total_volume0', 'B', 'rpm', 'p_atm', 'a4y')
            
disp('done')
