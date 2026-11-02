clear all; clc;
%% Slider-crank kinematics and cylinder volume, saved to part1.mat
load('constants.mat')

th = crank_angle;                   % degrees
r = crank_length;
l = coupler_length;
A_cyl = pi/4*bore^2;                % piston area (m^2)

%% Position, velocity and acceleration for both cylinder offsets
% The crank pin sits at (r cos th, r sin th); the piston pin slides on the line
% y = -e at x = s, so th = 0 is top dead centre. phi3 is the rod angle.
[s_length0, coupler_angle0, w3, alpha3, vel_piston0, acc_piston0] = slider_crank(th, r, l, r0, w2, alpha2);
[s_length1, coupler_angle1, w3_1, alpha3_1, vel_piston1, acc_piston1] = slider_crank(th, r, l, r1, w2, alpha2);

%% Accelerations of the centres of mass for the force analysis (part2.m)
% The force analysis stands the engine upright: Y runs up the cylinder axis
% (the x above) and X = -y, so gravity acts along -Y.
a2x = w2^2*OG2*sind(th) - alpha2*OG2*cosd(th);
a2y = -w2^2*OG2*cosd(th) - alpha2*OG2*sind(th);
a3x = w2^2*r*sind(th) + w3.^2*AG3.*sind(coupler_angle0) - alpha3*AG3.*cosd(coupler_angle0);
a3y = -w2^2*r*cosd(th) - w3.^2*AG3.*cosd(coupler_angle0) - alpha3*AG3.*sind(coupler_angle0);
a4y = acc_piston0;                  % piston acceleration along +Y

%% Check the analytic derivatives against numerical ones
dt = (pi/180)/w2;                   % seconds per crank degree
ddt = @(x) (circshift(x, [0 -1]) - circshift(x, [0 1]))/(2*dt);
err_v = max(abs(vel_piston0 - ddt(s_length0)))/max(abs(vel_piston0));
err_a = max(abs(acc_piston0 - ddt(vel_piston0)))/max(abs(acc_piston0));
err_w = max(abs(alpha3 - ddt(w3)))/max(abs(alpha3));
clear ddt
fprintf('Kinematics check against numerical derivatives (max relative error):\n');
fprintf('  piston velocity %.1e, piston acceleration %.1e, rod angular acceleration %.1e\n', err_v, err_a, err_w);

%% Stroke, swept volume and cylinder volume
s_max0 = max(s_length0); s_min0 = min(s_length0);
s_max1 = max(s_length1); s_min1 = min(s_length1);
stroke0 = s_max0 - s_min0;
stroke1 = s_max1 - s_min1;
vdisp0 = A_cyl*stroke0;             % swept volume of one cylinder (m^3)
vdisp1 = A_cyl*stroke1;
c_length0 = s_max0 - s_length0;     % piston distance from top dead centre (m)
c_length1 = s_max1 - s_length1;
total_volume0 = vdisp0/(rc0 - 1) + A_cyl*c_length0;
total_volume1 = vdisp1/(rc1 - 1) + A_cyl*c_length1;
fprintf('Stroke %.1f mm, swept volume %.1f cc per cylinder, %.2f L engine\n', ...
    stroke0*1e3, vdisp0*1e6, n_cyl*vdisp0*1e3);
fprintf('Compression ratio %.2f (offset cylinder %.2f)\n', ...
    max(total_volume0)/min(total_volume0), max(total_volume1)/min(total_volume1));

%% Plots
figure('Name', 'Piston velocity');
plot(th, vel_piston0, 'LineWidth', 2); hold on;
plot(th, vel_piston1, 'LineWidth', 2); hold off; grid on;
title('Piston velocity'); xlabel('Crank angle (degrees)'); ylabel('Velocity (m/s)');
legend('No offset', '18 mm offset'); xlim([0 720]); set(gca, 'XTick', 0:90:720);

figure('Name', 'Piston acceleration');
plot(th, acc_piston0, 'LineWidth', 2); hold on;
plot(th, acc_piston1, 'LineWidth', 2); hold off; grid on;
title('Piston acceleration'); xlabel('Crank angle (degrees)'); ylabel('Acceleration (m/s^2)');
legend('No offset', '18 mm offset'); xlim([0 720]); set(gca, 'XTick', 0:90:720);

figure('Name', 'Piston displacement');
plot(th, 100*c_length0, 'LineWidth', 2); hold on;
plot(th, 100*c_length1, 'LineWidth', 2); hold off; grid on;
title('Piston distance from top dead centre'); xlabel('Crank angle (degrees)'); ylabel('Distance (cm)');
legend('No offset', '18 mm offset'); xlim([0 720]); set(gca, 'XTick', 0:90:720);

figure('Name', 'Cylinder volume');
plot(th, total_volume0*1e6, 'LineWidth', 2); hold on;
plot(th, total_volume1*1e6, 'LineWidth', 2); hold off; grid on;
title('Cylinder volume'); xlabel('Crank angle (degrees)'); ylabel('Volume (cc)');
legend('No offset', '18 mm offset'); xlim([0 720]); set(gca, 'XTick', 0:90:720);

%% Data save
save('part1.mat', '-v7')
disp('done')
