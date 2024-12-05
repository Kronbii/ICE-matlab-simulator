clc;
load('part2.mat')
load('part1.mat')

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

% Plot the torque for a single cylinder and the total torque
figure;
hold on;
plot(crank_angle, T0, 'r');
xlabel('Crank Angle (degrees)');
ylabel('Torque (N-m)');
legend show;
title('Torque Plot');
grid on;
disp('doneee')
