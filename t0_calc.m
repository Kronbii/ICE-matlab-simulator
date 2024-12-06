clc;
load('part2.mat')

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
T0(1:300) =  0;
T0(540:720) = 0;
T1 = circshift(T0, 180);
T2 = circshift(T0, 360);
T3 = circshift(T0, 540);

T_total = T0 + T1 + T2 + T3;

% Plot the torque for a single cylinder and the total torque
figure;
hold on;
plot(crank_angle, T0, 'r');
% plot(crank_angle, T1, 'r');
% plot(crank_angle, T2, 'r');
% plot(crank_angle, T3, 'r');
xlabel('Crank Angle (degrees)');
ylabel('Torque (N-m)');
legend show;
title('Torque Plot');
grid on;

figure;
hold on;
plot(crank_angle, T_total, 'r');
xlabel('Crank Angle (degrees)');
ylabel('Torque total (N-m)');
legend show;
title('Torque Plot');
grid on;
disp('doneee')

% Plot coupler_angle0
figure;
plot(crank_angle, coupler_angle0);
title('Coupler Angle0');
xlabel('Crank Angle (degrees)');
ylabel('Coupler Angle0');

% Plot a4y
figure;
plot(crank_angle, a4y);
title('a4y');
xlabel('Crank Angle (degrees)');
ylabel('a4y (m/s^2)');

% Plot Fg
figure;
plot(crank_angle, Fg);
title('Fg');
xlabel('Crank Angle (degrees)');
ylabel('Fg (N)');

% Plot alpha3
figure;
plot(crank_angle, alpha3);
title('Alpha3');
xlabel('Crank Angle (degrees)');
ylabel('Alpha3 (rad/s^2)');

% Plot a3x
figure;
plot(crank_angle, a3x);
title('a3x');
xlabel('Crank Angle (degrees)');
ylabel('a3x (m/s^2)');

% Plot a3y
figure;
plot(crank_angle, a3y);
title('a3y');
xlabel('Crank Angle (degrees)');
ylabel('a3y (m/s^2)');

% Plot m3 (constant value)
figure;
plot(crank_angle, m3 * ones(size(crank_angle)));
title('m3 (Mass of Link 3)');
xlabel('Crank Angle (degrees)');
ylabel('m3 (kg)');

% Plot m2 (constant value)
figure;
plot(crank_angle, m2 * ones(size(crank_angle)));
title('m2 (Mass of Link 2)');
xlabel('Crank Angle (degrees)');
ylabel('m2 (kg)');

% Plot m4 (constant value)
figure;
plot(crank_angle, m4 * ones(size(crank_angle)));
title('m4 (Mass of Link 4)');
xlabel('Crank Angle (degrees)');
ylabel('m4 (kg)');

% Plot J3 (constant value)
figure;
plot(crank_angle, J3 * ones(size(crank_angle)));
title('J3 (Moment of Inertia of Link 3)');
xlabel('Crank Angle (degrees)');
ylabel('J3 (kg*m^2)');

% Plot gravitational constant g (constant value)
figure;
plot(crank_angle, g * ones(size(crank_angle)));
title('Gravitational Constant (g)');
xlabel('Crank Angle (degrees)');
ylabel('g (m/s^2)');

% Plot crank_length (constant value)
figure;
plot(crank_angle, crank_length * ones(size(crank_angle)));
title('Crank Length');
xlabel('Crank Angle (degrees)');
ylabel('Crank Length (m)');