clear all; clc;
%% Thermodynamic cycle, crank torque, flywheel, performance, speed sweep and gearbox
load('part1.mat')
c = load('part1.mat');                  % the same workspace as a struct, for the functions

%% Cycle at the design speed
cyc = engine_cycle(c, rpm);
pressure = cyc.pressure; temp = cyc.temp; mass = cyc.mass; xb = cyc.xb;
lift_in = cyc.lift_in; lift_ex = cyc.lift_ex; mdot_in = cyc.mdot_in; mdot_ex = cyc.mdot_ex;
dQ_comb = cyc.dQ_comb; dQ_wall = cyc.dQ_wall; burned = cyc.burned;
vol = total_volume0;
cycle_count = cyc.cycles;
fprintf('Cycle converged after %d cycles (state change %.1e)\n', cyc.cycles, cyc.change);
[p_peak, i_peak] = max(pressure);
fprintf('Peak pressure %.1f bar at %d degrees after TDC, peak temperature %.0f K\n', ...
    p_peak/1e5, i_peak - 360, max(temp));
fprintf('Residual gas fraction %.1f%%, exhaust temperature %.0f K\n', 100*cyc.residual, cyc.T_exhaust);

%% Crank torque: four cylinders 180 degrees apart, firing 1-3-4-2
Fg = pressure*A_cyl;
T0 = crank_torque(c, Fg, m2, m3, m4, J3);           % gas, inertia and gravity
T0_gas = crank_torque(c, Fg, 0, 0, 0, 0);           % gas alone
T0_inertia = T0 - T0_gas;                           % inertia and gravity
T1 = circshift(T0, [0 180]);
T2 = circshift(T0, [0 360]);
T3 = circshift(T0, [0 540]);
T_total = T0 + T1 + T2 + T3;
T_avg = mean(T0);
T_avg_total = mean(T_total);
W_loop = trapz([vol vol(1)], [pressure pressure(1)]);     % area of the sampled P-V loop
fprintf('Torque check: mean crank torque %.2f N m, P-V loop work x %d cylinders / 4 pi = %.2f N m\n', ...
    T_avg_total, n_cyl, n_cyl*W_loop/(4*pi));

%% Flywheel: the largest swing in shaft energy over the cycle sets its size
E_fly = cumsum(T_total - T_avg_total)*pi/180;       % energy above the mean load (J)
energy = max(E_fly) - min(E_fly);
Is = energy/(flywheel_cst*w2^2);                    % flywheel inertia (kg m^2)
flywheel_radius = ((2*Is)/(pi*flywheel_thickness*flywheel_density))^(1/4);
% the shaft speed this flywheel allows, against a steady load at the mean torque
w_fluc = sqrt(w2^2 + 2*(E_fly - mean(E_fly))/Is);
w_max = max(w_fluc); w_min = min(w_fluc);
cs_check = (w_max - w_min)/mean(w_fluc);

%% Performance at the design speed
Vd = vdisp0;
IMEP_g = cyc.W_gross/Vd; IMEP_n = cyc.W_net/Vd; PMEP = IMEP_g - IMEP_n;
FMEP = 1e5*polyval(fliplr(fmep_coef), rpm/1000);
BMEP = IMEP_n - FMEP;
engine_work = cyc.W_net;                            % net work per cylinder per cycle (J)
power_ind = n_cyl*IMEP_n*Vd*rpm/120;                % W
power_brake = n_cyl*BMEP*Vd*rpm/120;                % W
torque_brake = n_cyl*BMEP*Vd/(4*pi);                % N m
mech_efficiency = BMEP/IMEP_n;
fuel_conv_efficiency = cyc.W_gross/(cyc.m_fuel*Qhv);        % gross indicated
brake_efficiency = BMEP*Vd/(cyc.m_fuel*Qhv);
fuel_mf = n_cyl*cyc.m_fuel*rpm/120;                 % kg/s
BSFC = fuel_mf/power_brake*3.6e9;                   % g/kWh
air_density = Pi/(R*Taf);
vol_efficiency = cyc.m_air/(air_density*Vd);
piston_speed = 2*stroke0*rpm/60;
heat_to_walls = cyc.Q_wall/cyc.Q_comb;
Engine_power = power_brake/1e3;                     % kW

fprintf('\nPerformance at %d rpm (whole engine unless noted)\n', rpm);
fprintf('Net work per cylinder per cycle: %.1f J\n', engine_work);
fprintf('IMEP gross %.2f bar, net %.2f bar, pumping %.2f bar\n', IMEP_g/1e5, IMEP_n/1e5, PMEP/1e5);
fprintf('FMEP %.2f bar, BMEP %.2f bar\n', FMEP/1e5, BMEP/1e5);
fprintf('Indicated power %.1f kW, brake power %.1f kW\n', power_ind/1e3, power_brake/1e3);
fprintf('Brake torque %.1f N m, mechanical efficiency %.1f%%\n', torque_brake, 100*mech_efficiency);
fprintf('Fuel conversion efficiency: indicated %.1f%%, brake %.1f%%\n', ...
    100*fuel_conv_efficiency, 100*brake_efficiency);
fprintf('Heat to the walls: %.1f%% of the fuel energy released\n', 100*heat_to_walls);
fprintf('BSFC %.0f g/kWh, fuel flow %.2f g/s\n', BSFC, 1e3*fuel_mf);
fprintf('Volumetric efficiency %.1f%%, mean piston speed %.1f m/s\n', 100*vol_efficiency, piston_speed);
fprintf('Flywheel: energy swing %.1f J, inertia %.4f kg m^2, radius %.1f cm (speed swing %.1f%%)\n', ...
    energy, Is, 100*flywheel_radius, 100*cs_check);

%% Full-load torque and power from 1000 to 6000 rpm
n_sw = numel(rpm_sweep);
sweep_torque = zeros(1, n_sw); sweep_power = sweep_torque; sweep_bmep = sweep_torque;
sweep_eta_v = sweep_torque; sweep_eta_b = sweep_torque; sweep_bsfc = sweep_torque;
fprintf('\n  rpm   torque   power   BMEP   vol.eff  BSFC\n');
for j = 1:n_sw
  sw = engine_cycle(c, rpm_sweep(j));
  fmep_j = 1e5*polyval(fliplr(fmep_coef), rpm_sweep(j)/1000);
  bmep_j = sw.W_net/Vd - fmep_j;
  sweep_bmep(j) = bmep_j;
  sweep_torque(j) = n_cyl*bmep_j*Vd/(4*pi);
  sweep_power(j) = n_cyl*bmep_j*Vd*rpm_sweep(j)/120;
  sweep_eta_v(j) = sw.m_air/(air_density*Vd);
  sweep_eta_b(j) = bmep_j*Vd/(sw.m_fuel*Qhv);
  sweep_bsfc(j) = 3.6e9*sw.m_fuel/(bmep_j*Vd);
  fprintf('%5d  %6.1f  %6.1f  %5.2f   %5.1f%%  %4.0f\n', rpm_sweep(j), sweep_torque(j), ...
      sweep_power(j)/1e3, bmep_j/1e5, 100*sweep_eta_v(j), sweep_bsfc(j));
end
[torque_peak, j] = max(sweep_torque); rpm_torque_peak = rpm_sweep(j);
[power_peak, j] = max(sweep_power); rpm_power_peak = rpm_sweep(j);
fprintf('Peak torque %.1f N m at %d rpm, peak power %.1f kW at %d rpm\n', ...
    torque_peak, rpm_torque_peak, power_peak/1e3, rpm_power_peak);

%% Eight-speed gearbox for the vehicle in constants.m
% Top gear from the top speed, first gear from the steepest demand; the ratios in
% between follow i_n = i_top phi1^(z-n) phi2^((z-n)(z-n-1)/2) with gear 6 direct.
iAmin = (3.6*(pi/30)*radius_dynamic*rpm_max)/vel_max_car/iD;
iAmax = ((car_mass*g*(friction_coef*cosd(incl_angle) + sind(incl_angle)) + ...
    car_mass*lambda*car_acc)*radius_dynamic)/(total_efficiency*car_max_torque)/iD;
z = gears_total; d = direct_drive;
% the two conditions are linear in log(phi1) and log(phi2)
M = [z - d, (z - d)*(z - d - 1)/2; z - 1, (z - 1)*(z - 2)/2];
x = M\[-log(iAmin); log(iAmax/iAmin)];
phi1 = exp(x(1)); phi2 = exp(x(2));
k = 1:z;
iA = iAmin*phi1.^(z - k).*phi2.^((z - k).*(z - k - 1)/2);
v_gear = 3.6*(rpm_max*pi/30)*radius_dynamic./(iA*iD);    % km/h at rpm_max
fprintf('\nGearbox for a %d kg vehicle and a %d N m engine: phi1 = %.4f, phi2 = %.4f\n', ...
    car_mass, car_max_torque, phi1, phi2);
fprintf('Gear ratios: %s\n', sprintf('%.2f ', iA));
fprintf('Speed at %d rpm (km/h): %s\n', rpm_max, sprintf('%.0f ', v_gear));

%% Plots
figure('Name', 'Pressure');
plot(crank_angle, pressure/1e5, 'LineWidth', 2); grid on;
title('Cylinder pressure'); xlabel('Crank angle (degrees)'); ylabel('Pressure (bar)');
xlim([0 720]); set(gca, 'XTick', 0:90:720);

figure('Name', 'P-V diagram');
loglog(vol*1e6, pressure/1e5, 'LineWidth', 2); grid on;
title('Pressure-volume diagram (log scales)'); xlabel('Volume (cc)'); ylabel('Pressure (bar)');

figure('Name', 'Temperature');
plot(crank_angle, temp, 'LineWidth', 2); grid on;
title('Gas temperature'); xlabel('Crank angle (degrees)'); ylabel('Temperature (K)');
xlim([0 720]); set(gca, 'XTick', 0:90:720);

figure('Name', 'Heat release');
plot(crank_angle, dQ_comb, 'LineWidth', 2); hold on;
plot(crank_angle, dQ_wall, 'LineWidth', 2); hold off; grid on;
title('Heat release and heat lost to the walls'); xlabel('Crank angle (degrees)'); ylabel('J per degree');
legend('Heat release', 'Heat to the walls'); xlim([300 540]);

figure('Name', 'Gas exchange');
plot(crank_angle, 1e3*mdot_in, 'LineWidth', 2); hold on;
plot(crank_angle, -1e3*mdot_ex, 'LineWidth', 2); hold off; grid on;
title('Mass flow through the valves'); xlabel('Crank angle (degrees)'); ylabel('Mass flow (g/s)');
legend('Into the cylinder through the intake', 'Out through the exhaust'); xlim([0 720]); set(gca, 'XTick', 0:90:720);

figure('Name', 'Torque, one cylinder');
plot(crank_angle, T0, 'LineWidth', 2); hold on;
plot(crank_angle, T0_gas, '--', 'LineWidth', 1.5);
plot(crank_angle, T0_inertia, ':', 'LineWidth', 1.5); hold off; grid on;
title('Crank torque of one cylinder'); xlabel('Crank angle (degrees)'); ylabel('Torque (N m)');
legend('Total', 'Gas', 'Inertia and gravity'); xlim([0 720]); set(gca, 'XTick', 0:90:720);

figure('Name', 'Torque, four cylinders');
plot(crank_angle, T_total, 'LineWidth', 2); hold on;
plot([0 720], [T_avg_total T_avg_total], '--', 'LineWidth', 1.5); hold off; grid on;
title('Crank torque of the engine'); xlabel('Crank angle (degrees)'); ylabel('Torque (N m)');
legend('Four cylinders', 'Mean'); xlim([0 720]); set(gca, 'XTick', 0:90:720);

figure('Name', 'Flywheel');
plot(crank_angle, E_fly, 'LineWidth', 2); grid on;
title('Shaft energy above the mean load'); xlabel('Crank angle (degrees)'); ylabel('Energy (J)');
xlim([0 720]); set(gca, 'XTick', 0:90:720);

figure('Name', 'Crank speed');
plot(crank_angle, w_fluc*30/pi, 'LineWidth', 2); hold on;
plot([0 720], [rpm rpm], '--'); hold off; grid on;
title('Crank speed with the flywheel'); xlabel('Crank angle (degrees)'); ylabel('Speed (rpm)');
xlim([0 720]); set(gca, 'XTick', 0:90:720);

figure('Name', 'Torque and power');
plot(rpm_sweep, sweep_torque, '-o', 'LineWidth', 2); hold on;
plot(rpm_sweep, sweep_power/1e3, '-s', 'LineWidth', 2); hold off; grid on;
title('Full-load torque and power'); xlabel('Engine speed (rpm)'); ylabel('Torque (N m) and power (kW)');
legend('Brake torque', 'Brake power');

figure('Name', 'Gear ratios');
bar(iA); grid on;
title('Gear ratios'); xlabel('Gear'); ylabel('Ratio');

%% Data save
clear c cyc sw
save('part2.mat', '-v7')
disp('done')
