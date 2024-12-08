clear all; clc;

%% Load Necessary Data
load('part5.mat')

car_max_torque = max(T0);
iA = zeros(1, length(gears_total));

%% Calculate iAmin
iAmin = (3.6 * (pi / 30) * radius_dynamic * rpm_max)/(vel_max_car);

%% Calculate iAmax
iAmax = ((car_mass*g*(friction_coef*cosd(incl_angle) + sind(incl_angle)) + ...
          car_mass*lambda*car_acc) * radius_dynamic) / (total_efficiency * car_max_torque);
      
%% Divide by the differential ratio
iAmin = iAmin / iD;
iAmax = iAmax / iD;

%% Calculate phi1 and phi2 
equations = @(phi) [(iAmin * (phi(1)^(gears_total - direct_drive)) * phi(2)^(0.5 * (gears_total - direct_drive) * (gears_total - direct_drive - 1))) - 1;
                    (iAmin * (phi(1)^(gears_total - 1)) * phi(2)^(0.5 * (gears_total - 1) * (gears_total - 2))) - iAmax];


%% Solve for phi1 and phi2
solution = fsolve(equations, [1 ; 1], optimset('Display', 'iter'));

phi1=solution(1);
phi2=solution(2);
 
disp(['phi1 = ', num2str(phi1)]);
disp(['phi2 = ', num2str(phi2)]);
disp(['iAmin = ', num2str(iAmin)]);
disp(['iAmax = ', num2str(iAmax)]);

%gear ratios 
for i=1:gears_total
    iA(i)=((iAmin * phi1^(gears_total-i))*phi2^(0.5*(gears_total-i)*(gears_total-i-1)));
end
disp('Gear Ratios (iA):');
disp(iA);