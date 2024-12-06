clear all; clc

%% Load Necessary Data
load('part1.mat') 

%% Variables and constants
%%Constants
vol = total_volume0; 
Taf = 300; T_inital = 300; %% kelvin
yi =  1.3;
yc = 1.25;
ye = 1.48;
R = 287; %% gas constant of air (J/kg * Kelvin)
phi = 1;
CBF = 0.96; 
AFR = 14.7; 
Qhv = 43.5e6; %% J/kg
cv = R/(yc - 1); %% J/kg*k
a = 5;
m = 2;
 
%%Per Engine
d_in_valve = 0.4*bore; %% mm
d_out_valve = 0.4*bore; %% mm
dt = (1/(rpm*2*pi/60))*(pi/180); %% seconds
Cd = 0.7; %% check unit
At = (pi/4)*(d_in_valve)^2; 
Pe = 1 * p_atm; %% pascal
Pi = 1 * p_atm; %% pascal
A_cyl = (pi/4) * bore^2;

intake_start = 1;
intake_end = 222;
compression_start = 223;
compression_end = 329;
combustion_start = 330;
combustion_end = 389;
expansion_start = 390;
expansion_end = 487;
exhaust_start = 488;
exhaust_end = 720;

%%Variables
temp = zeros(1, length(crank_angle)); 
pressure = zeros(1, length(crank_angle));
mass = zeros(1, length(crank_angle));
xb = zeros(1, combustion_end - combustion_start);

error = inf;
tolerence= 1e-3;
counter = 0;

temp(compression_start) = T_inital; 

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
    
    %pressure(exhaust_start)= (((temp(expansion_end) * mi *R)/vol(exhaust_start))^(ye))*(1/(pressure(expansion_end)^(ye-1)));
    %temp(exhaust_start)= (pressure(exhaust_start)*vol(exhaust_end)/(mi*R));

    for i = (exhaust_start + 1):exhaust_end
      mass(i)=real(mass(i-1)-(dt*((Cd*At*pressure(i-1))/(sqrt((pressure(i-1)*vol(i-1))/mass(i-1))))*((Pe/pressure(i-1))^(1/ye))*((((2*ye)/(ye-1))*(1-(Pe/pressure(i-1))^((ye-1)/ye)))^0.5)));
      pressure(i) = ((((temp(i-1))*(mass(i))*R)/(vol(i)))^ye) * (1/((pressure(i-1))^(ye-1)));
      temp(i) = (pressure(i)*(vol(i)))/(R*mass(i));
    end

    %% Intake Numerical Analysis
    mass(intake_start) = mass(exhaust_end);
    
    pressure(intake_start) = pressure(exhaust_end);
    temp(intake_start) = temp(exhaust_end);
    
    for i = (intake_start + 1):intake_end    
        mass(i)=real(mass(i-1)+(dt*((Cd*At*pressure(i-1))/(sqrt(R*Taf)))*((pressure(i-1)/Pi)^(1/yi))*((((2*yi)/(yi-1))*(1-(pressure(i-1)/Pi)^((yi-1)/yi))).^0.5)));
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

F_atm = p_atm * A_cyl;
Fg = pressure * A_cyl;


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

%% saving data
save('part2.mat', 'pressure', 'F_atm', 'Fg')
disp('donee')
disp('stable1.1')