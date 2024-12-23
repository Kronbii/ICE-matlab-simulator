clear all; clc

%% Load Necessary Data
load('part1.mat')

%% Variables
vol = total_volume0; 
temp = zeros(1, length(crank_angle)); 
pressure = zeros(1, length(crank_angle));
mass = zeros(1, length(crank_angle));
xb = zeros(1, combustion_end - combustion_start);

temp(compression_start) = T_initial; 

%% Numerical Analysis
while error > tolerence
    counter = counter + 1;
    %% Compression Numerical Analysis
    pressure(compression_start) = Pi;
    mi = (pressure(compression_start) * vol(compression_start)) / (R * temp(compression_start));
    mf = mi / (1 + AFR/phi);
    mass(compression_start:expansion_end) = mi;
    for i = (compression_start + 1):compression_end
        pressure(i) = (pressure(compression_start) * (vol(compression_start))^yi) / ((vol(i))^yi);
        temp(i) = (pressure(i) * vol(i)) / (mi * R);
    end

    %% Combustion Numerical Analysis
    pressure(combustion_start) = pressure(compression_end);
    temp(combustion_start) = temp(compression_end);
    
    temp(combustion_end) = (CBF * mf * Qhv) / (mi * cv) + temp(combustion_start);
    pressure(combustion_end) = (mi * R * temp(combustion_end)) / (vol(combustion_end));
    
    for i = (combustion_start):(combustion_end - 1)
        xb(i) = 1 - exp(-a * ((crank_angle(i) - combustion_start) / (combustion_end - combustion_start))^(m + 1));
        pressure(i) = (xb(i) * ((pressure(combustion_end) * vol(combustion_end)^yc) - ...
            (pressure(combustion_start) * vol(combustion_start)^yc)) + ...
            pressure(combustion_start) * vol(combustion_start)^yc) / (vol(i)^yc);
        temp(i) = (pressure(i)*(vol(i)))/(R*mass(i));
    end

    %% Expansion Numerical Analysis
    for i = (expansion_start):expansion_end
        pressure(i) = (pressure(combustion_end) * vol(combustion_end)^ye) / (vol(i)^ye);
        temp(i) = (pressure(i) * vol(i)) / (mi * R);
    end
    
    %% Exhaust Numerical Analysis
    mass(exhaust_start) = mi;
    
    pressure(exhaust_start) = pressure(expansion_end);
    temp(exhaust_start) = temp(expansion_end);

    for i = (exhaust_start + 1):exhaust_end
      mass(i)=real(mass(i-1)-(dt*((Cd*At*pressure(i-1))/(sqrt((pressure(i-1)*vol(i-1))/mass(i-1))))* ...
      ((Pe/pressure(i-1))^(1/ye))*((((2*ye)/(ye-1))*(1-(Pe/pressure(i-1))^((ye-1)/ye)))^0.5)));
      pressure(i) = ((((temp(i-1))*(mass(i))*R)/(vol(i)))^ye) * (1/((pressure(i-1))^(ye-1)));
      temp(i) = (pressure(i)*(vol(i)))/(R*mass(i));
    end

    %% Intake Numerical Analysis
    mass(intake_start) = mass(exhaust_end);
    pressure(intake_start) = pressure(exhaust_end);
    temp(intake_start) = temp(exhaust_end);
    
    for i = (intake_start + 1):intake_end    
        mass(i)=real(mass(i-1)+(dt*((Cd*At*pressure(i-1))/(sqrt(R*Taf)))*((pressure(i-1)/Pi)^(1/yi))* ...
        ((((2*yi)/(yi-1))*(1-(pressure(i-1)/Pi)^((yi-1)/yi))).^0.5)));
        xr = mass(exhaust_end)/mass(i);
        temp(i) = xr*temp(exhaust_end) + ((1- xr))*Taf;
        pressure(i) = (mass(i)*R*temp(i))/vol(i);
        %pressure(i) = Pi;
    end

    % Calculate error based on the temperature change
    error = abs(temp(intake_end) - temp(compression_start));
    temp(compression_start) = temp(intake_end);
    
    % Debugging output
    fprintf("Iteration: %d, Error: %.6f, Temp: %.2f\n", counter, error, temp(intake_end));
end

%% Force Calculation
fprintf("Converged with error %.6f after %d iterations\n", error, counter);

%% Torque %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Torque Calculation
F_atm = p_atm * A_cyl;
Fg = pressure * A_cyl;

T0 = zeros(1,length(crank_angle));
a4y = -a4y;

for k= 1:length(crank_angle)
A1=[0 0 0 0 1 0 1 0; 
     0 0 0 0 0 1 0 0;
     0 0 1 0 -1 0 0 0;
     0 0 0 1 0 -1 0 0;
     0 0 AG3*cosd(coupler_angle0(k)) AG3*sind(coupler_angle0(k)) BG3*cosd(coupler_angle0(k)) BG3*sind(coupler_angle0(k)) 0 0;
     1 0 -1 0 0 0 0 0;
     0 1 0 -1 0 0 0 0;
     0 0 crank_length*cosd(crank_angle(k)) crank_length*sind(crank_angle(k)) 0 0 0 1];

B1=[0;
    m4*a4y(k)+Fg(k)-F_atm+m4*g;
    m3*a3x(k);
    m3*a3y(k)+m3*g;
    J3*alpha3(k);
    m2*a2x(k);
    m2*a2y(k)+m2*g;
    -(m2*g)*OG2*sind(crank_angle(k))];
z = linsolve(A1,B1);
T0(k) = -z(8);
end

T1 = circshift(T0, 180);
T2 = circshift(T0, 360);
T3 = circshift(T0, 540);

T_total = T0 + T1 + T2 + T3;

%% Flywheel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating Mean Torque
T_avg = mean(T0); %per cylinder
T_avg_total = mean(T_total); %total for engine

%% Calculating the crank angles at which the avg and instantaneous torque intersect
theta_int = find(abs(diff(sign(T0 - T_avg))) > 0);
disp('Intersection Angles (degrees):')
disp(theta_int)
theta_int = theta_int(theta_int >= combustion_start);
disp('Intersection Angles (degrees):')
disp(theta_int)

%% Calculating the Energy using trapz
integ_range = min(theta_int):max(theta_int);
integ_range_rad = integ_range * pi/180;

energy = trapz(integ_range_rad, T0(integ_range) - T_avg);

%% Calculating Interia Using Energy
Is = energy/(flywheel_cst * w2^2);

%% Calculating The Flywheel Radius
flywheel_radius = ((2 * Is) / (pi * flywheel_thickness * flywheel_density))^(1/4);

%% Calculating minimum and maximum rotational speed
w_max = (w2/2) * (flywheel_cst+2);
w_min = w2*(2-(flywheel_cst+2)/2);

% Generate the sine wave
w_fluc = ((w_max - w_min)/2) * sind(crank_angle) + w2;

%% Perfomance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating Performance Characteristics
engine_work = trapz(vol, pressure); % work/cycle
power_ind = ((engine_work*(rpm/60))/nr); % indicated power in W
MEP = (power_ind*nr)/(vdisp0*1000*rpm/60); % mean effective pressure in KPa
power_fr = Pf*power_ind; % Friction power
mech_efficiency = (1-(power_fr/power_ind)); % mechanical efficiency
power_brake = mech_efficiency*power_ind; % Brake power in W
BMEP = (power_brake*nr)/(vdisp0*1000*rpm/60); % brake mean effective pressure in KPa
thermo_torque = (BMEP*vdisp0*1000)/(6.28*nr); % thermo torque in Nm
fuel_mf = (mf*(rpm/60))/nr; % fuel mass flow rate in kg/s
fuel_conv_efficiency = (engine_work)/((mf)*(Qhv)); % Fuel conversion efficiency
SFC = (fuel_mf*10^(-3))/(power_ind*10^(-3)); % specific fuel consumption in mg/J
vol_efficiency = (max(mass)/(air_density*vdisp0)); % Volumetric efficiency
piston_speed = 2*L*rpm/60; % piston speed in m/sec
Engine_power = (power_ind*4)*(10^-3); % kW

%% Displaying Performance Characteristics
fprintf('Engine Work per Cycle: %.4f J\n', engine_work);
fprintf('Indicated Power: %.4f W\n', power_ind);
fprintf('Mean Effective Pressure (MEP): %.4f kPa\n', MEP);
fprintf('Friction Power: %.4f W\n', power_fr);
fprintf('Mechanical Efficiency: %.4f\n', mech_efficiency);
fprintf('Brake Power: %.4f W\n', power_brake);
fprintf('Brake Mean Effective Pressure (BMEP): %.4f kPa\n', BMEP);
fprintf('Thermodynamic Torque: %.4f Nm\n', thermo_torque);
fprintf('Fuel Mass Flow Rate: %.4f kg/s\n', fuel_mf);
fprintf('Fuel Conversion Efficiency: %.4f\n', fuel_conv_efficiency);
fprintf('Specific Fuel Consumption (SFC): %.10f mg/J\n', SFC);
fprintf('Volumetric Efficiency: %.4f\n', vol_efficiency);
fprintf('Piston Speed: %.4f m/s\n', piston_speed);
fprintf('Engine Power: %.4f kW\n\n', Engine_power);
fprintf('Energy of flywheel: %.4f J\n', energy);
fprintf('Flywheel Radius: %.4f m\n', flywheel_radius);
fprintf('Average Torque: %.4f Nm\n\n', T_avg);

%% Transmission %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
equations = @(phi) [(iAmin * (phi(1)^(gears_total - direct_drive)) * phi(2)^(0.5 * (gears_total - direct_drive) * ...
(gears_total - direct_drive - 1))) - 1;
                    (iAmin * (phi(1)^(gears_total - 1)) * phi(2)^(0.5 * (gears_total - 1) * (gears_total - 2))) - iAmax];


%% Solve for phi1 and phi2
solution = fsolve(equations, [1 ; 1], optimset('Display', 'iter'));

phi1=solution(1);
phi2=solution(2);
 
disp(['phi1 = ', num2str(phi1)]);
disp(['phi2 = ', num2str(phi2)]);
disp(['iAmin = ', num2str(iAmin)]);
disp(['iAmax = ', num2str(iAmax)]);

%% Calculating gear ratios 
for i=1:gears_total
    iA(i)=((iAmin * phi1^(gears_total-i))*phi2^(0.5*(gears_total-i)*(gears_total-i-1)));
end

%% Displating Gear ratios
disp('Gear Ratios (iA):');
disp(iA);

%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting Pressure vs theta
figure;  % Create a new figure window
plot(crank_angle, pressure/p_atm);
xlabel('Crank Angle');  % Label for x-axis
ylabel('Pressure/Patm'); % Label for y-axis 
title('Pressure vs Volume'); % Title of the plot
grid on;  % Add grid for better readability

%% Plotting Pressure vs Volume
figure;  % Create a new figure window
plot(vol*10^6, pressure/p_atm);
xlabel('Volume/Vmax');  % Label for x-axis
ylabel('Pressure/Patm'); % Label for y-axis 
title('Pressure vs Volume'); % Title of the plot
grid on;  % Add grid for better readability

%% Plotting the torque for a single cylinder
figure;
hold on;
plot(crank_angle, T0, 'r');
xlabel('Crank Angle (degrees)');
ylabel('Torque (N-m)');
legend show;
title('Torque Plot');
grid on;

%% Plotting the total torque
figure;
hold on;
plot(crank_angle, T_total, 'r');
xlabel('Crank Angle (degrees)');
ylabel('Torque total (N-m)');
legend show;
title('Torque Plot');
grid on;

%% Plotting the Torque with and without a flywheel
figure;
plot(crank_angle, T0, 'r', 'LineWidth', 1);
hold on;
plot(crank_angle, (T_avg), '-ob', 'MarkerSize', 1, 'LineWidth', 1);  
title('Torque Before and After Adding Flywheel', 'FontSize', 14);
xlabel('Crank Angle (degrees)', 'FontSize', 12);
ylabel('Torque (Nm)', 'FontSize', 12);
grid on;
legend('Torque Before Flywheel', 'Torque With Flywheel', 'Location', 'Best');
yline(0, 'k--'); 
hold off;


%% Plotting the RPM wave
figure;
plot(crank_angle, w_fluc);
yline(w2);
xlabel('Crank Angle');
ylabel('Amplitude');
title('RPM fluctuation Centered Around w2');
grid on;

%% saving data
save('part2.mat')