clear all; clear; clc

%% User Input
load('data.mat')
while true 
    condition_input = input('Select baseline or offset conditions: ', 's'); 
    condition_input = lower(condition_input);  
    
    if strcmp(condition_input, 'baseline') 
        vol = total_volume0 * 10^(-6); 
        disp('Using Baseline Condition'); 
        break; 
    elseif strcmp(condition_input, 'offset') 
        vol = total_volume1 * 10^(-6); 
        disp('Using Offset Condition'); 
        break; 
    else 
        disp('Invalid input. Please enter either "Baseline" or "Offset".'); 
    end 
end


%% Variables and constants
%%Constants
p_atm = 101.325 * 10^3; %% pascal
Taf = 303; %% kelvin
yi =  1.3;
yc = 1.25;
ye = 1.48;
R = 8.314; %%287.05; %% gas constant of air (J/kg * Kelvin)
phi = 1;
pressure_offset = 5325; %% pascal
CBF = 0.95; 
AFR = 14.7; 
Qhv = 44.4*10^6; %% j/kg
cv = 850; %% J/kgk
a = 5;
m = 2;
 
%%Per Engine
bore = 84/1000; %% mm
d_in_valve = 33/1000; %% mm
d_out_valve = 30.5/1000; %% mm
rpm = 2200;
dt = 1/(rpm * 6);
Cd = 0.7; %% check unit
Ati = 0.785398*(d_in_valve)^2; 
Ate = 0.785398*(d_out_valve)^2;
Pe = p_atm + pressure_offset; %% pascal
Pi = p_atm - pressure_offset; %% pascal

intake_start = 1; %%%%CHANGE
intake_end = 240; %%%%CHANGE
compression_start = 240; %%%%CHANGE
compression_end = 346; %%%%CHANGE
combustion_start = 346; %%%%CHANGE
combustion_end = 386; %%%%CHANGE
expansion_start = 386; %%%%CHANGE
expansion_end = 492; %%%%CHANGE
exhaust_start = 492; %%%%CHANGE
exhaust_end = 720; %%%%CHANGE

%%Variables
temp = zeros(1, length(crank_angle)); 
pressure = zeros(1, length(crank_angle));
mass = zeros(1, length(crank_angle));
xb = zeros(1, combustion_end-combustion_start);

pressure(compression_start) = p_atm; %pascal
temp(compression_start) = Taf; %degrees celsius %%%%CHANGE
error = inf;
titi = 0;

while error > 1
    titi = titi + 1;
    %% Compression Numerical Analysis
    mi = (pressure(compression_start) * vol(compression_start)) / (R * temp(compression_start));
    mf = (mi) / (1 + AFR/phi);
    mass(1, compression_start:expansion_end) = mi;
    for i=(compression_start + 1):compression_end
        %%Calculating instanteneous pressure
        pressure(i) = (pressure(compression_start)*vol(compression_start)^yi) / (vol(i)^yi);

        %%Calculating temperature
        temp(i) = (pressure(i) * vol(i))/ (mi * R);
    end

    %% Combustion Numerical Analysis
    temp(combustion_end) = (CBF * mf * Qhv)/(mi * cv) + temp(combustion_start);
    pressure(combustion_end) = (mi * R * temp(combustion_end)) / (vol(combustion_end));
    for i=(combustion_start + 1):(combustion_end - 1)
        xb(i) = 1- exp((-a)*((i-combustion_start)/(combustion_end - combustion_start))^(m+1));
        pressure(i) = (xb(i) * ((pressure(combustion_end)*vol(combustion_end)^yc) - ...
       (pressure(combustion_start)*vol(combustion_start)^yc)) + ...
        pressure(combustion_start)*vol(combustion_start)^yc)/(vol(i)^yc);
    end

    %% Expansion Numerical Analysis
    for i=(expansion_start + 1):expansion_end
        %%Calculating instanteneous pressure
        pressure(i) = (pressure(expansion_start) * vol(expansion_start)^ye) / (vol(i)^ye);

        %%Calculating temperature
        temp(i) = (pressure(i) * vol(i))/ (mi * R);
    end
    
%% Exhaust Numerical Analysis
for i = (exhaust_start + 1):exhaust_end
    mass(i) = mass(i - 1) - dt * ((Cd * Ate * pressure(i - 1)) / ...
        sqrt((pressure(i - 1) * vol(i - 1)) / mass(i - 1))) * ...
        ((Pe / pressure(i - 1))^(1 / ye)) * ...
        ((2 * ye) / (ye - 1) * (1 - (Pe / pressure(i - 1))^((ye - 1) / ye)))^0.5;

    pressure(i) = (((temp(i - 1) * mass(i) * R) / vol(i))^ye) * ...
        (pressure(i - 1))^(1 - ye);

    temp(i) = (pressure(i) * vol(i)) / (mass(i) * R);
end

    %% Intake Numerical Analysis
    mass(intake_start) = mass(exhaust_end);
    pressure(intake_start:intake_end) = pressure(exhaust_end);
    temp(intake_start) = temp(exhaust_end);
    
    for i=(intake_start + 1):intake_end    
       mass(i) = mass(i-1) + dt * ((Cd*Ati*Pi)/(sqrt(R*Taf))) * ((pressure(i-1)/Pi)^(1/yi)) * ...
       (((2*yi)/(yi-1)) * (1-((2*yi)/(yi-1))^((yi-1)/(yi))))^0.5;

       xr = mass(exhaust_end)/mass(i);

       temp(i) = xr * temp(exhaust_end) + (1-xr) * temp(compression_start);
    end
    error = abs(temp(intake_end) - Taf)
    Taf = temp(intake_end);
    titi
end
%% Plotting Pressure vs Volume
figure;  % Create a new figure window
plot(vol/max(vol), pressure/p_atm);
xlabel('Volume/Vmax');  % Label for x-axis
ylabel('Pressure/Patm'); % Label for y-axis 
title('Pressure vs Volume'); % Title of the plot
grid on;  % Add grid for better readability