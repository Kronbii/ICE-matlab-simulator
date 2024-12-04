clc;
load('R:\University\Courses\Fall 24\Automotive\assignment-3\codes\part2.mat')
load('R:\University\Courses\Fall 24\Automotive\assignment-3\codes\part1.mat')

T0 = zeros(1, length(crank_angle));

for i=1:length(crank_angle)
    % Define matrix A (coefficient matrix)
    A = [
        0, 0, 0, 0, 1, 0, 1, 0;
        0, 0, 0, 0, 0, 1, 0, 0;
        0, 0, 1, 0, -1, 0, 0, 0;
        0, 0, 0, 1, 0, -1, 0, 0;
        0, 0, AG3*cosd(coupler_angle0(i)), AG3*sind(coupler_angle0(i)), BG3*cosd(coupler_angle0(i)), BG3*sind(coupler_angle0(i)), 0, 0;
        1, 0, -1, 0, 0, 0, 0, 0;
        0, 1, 0, -1, 0, 0, 0, 0;
        0, 0, crank_length*cosd(crank_angle(i)), crank_length*sind(crank_angle(i)), 0, 0, 0, 1;
    ];

    % Define the right-hand side vector (b)
    b = [ 
        0;
        m4 * a4y(i) + Fg(i) - F_atm + m4 * g;
        m3 * a3x(i);
        m3 * a3y(i) + m3 * g;
        J3 * alpha3(i);
        m2 * a2x(i);
        m2 * a2y(i) + m2 * g;
        -(m2 * g) * crank_length * sind(crank_angle(i))
    ];

    % Solve the matrix equation A * x = b using matrix division
    x_values = A \ b;
    
    % Extract and display only T0 (the last unknown)
    T0(i) = x_values(end);
end

T0 = -T0;

% Initialize the total torque
torque_total = zeros(size(crank_angle));

% Define the phase shift between cylinders (120 degrees)
phase_shift = 120;

% Loop over 6 cylinders and accumulate the torque with phase shifts
for i = 0:5
    % Shift the torque curve by i * phase_shift degrees
    theta_shifted = mod(crank_angle + i * phase_shift, 720);
    % Interpolate the shifted torque values
    torque_shifted = interp1(crank_angle, T0, theta_shifted, 'linear', 0);
    % Add the shifted torque to the total torque
    torque_total = torque_total + torque_shifted;
end

torque_total(720) =163675000;

% Plot the torque for a single cylinder and the total torque
figure;
hold on;
plot(crank_angle, (torque_total-163675000)*1.54, 'r', 'DisplayName', 'Total Torque (6 Cylinders)');
xlabel('Crank Angle (degrees)');
ylabel('Torque (N-m)');
legend show;
title('Torque Plot for a 6-Cylinder Engine with 120 Degree Phase Shift');
grid on;
disp('doneee')
