%% Variables and constants
%%Constants
vol = total_volume0; 
Taf = 300; %% kelvin
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
tolerence= 1e-2;
counter = 0;