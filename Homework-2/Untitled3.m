%%Variables
vol = total_volume0 (350:530);
mi = 0:719;
temp = 0:719;

%%Constants
p_atm = 10^5; %pascal
p_start = p_atm; %pascal
t_start = 300; %degrees celsius
yi =  1.3;
R = 8.3; %%%%CHANGE
pv_cst = p_start * vol(1)^yi;
m_start = (p_start * vol(1)) / (R * t_start);
kiki = 0;
%%Unknowns
p_end = 0;

for i=1:length(vol)
    %%Calculating instanteneous pressure
    pressure(i) = pv_cst / (vol(i)^yi);
    
    %%Calculating temperature
    temp(i) = (pressure(i) * vol(i))/ (m_start * R);
    
    %%Calculating mass of fuel air mixture
    mi(i) = (pressure(i) * vol(i))/ (R * temp(i));
    kiki = kiki + 1
end

%% Plotting Pressure vs Volume
figure;  % Create a new figure window
plot(vol, pressure);
xlabel('Volume (m³)');  % Label for x-axis
ylabel('Pressure (Pa)'); % Label for y-axis
title('Pressure vs Volume'); % Title of the plot
grid on;  % Add grid for better readability
   