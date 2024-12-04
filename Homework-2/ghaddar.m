Bore = 0.084;  
throat_area = ((pi * (0.39* Bore)^2)/4); 
rpm = 1500; 
theta = 0:719; 
load('baselineData.mat');  
while true 
condition_input = input('Baseline or Offset condition: ', 's'); 
condition_input = lower(condition_input);  
if strcmp(condition_input, 'baseline') 
V = (v_baseline_values + vc_base) * 10^-6; 
disp('Using Baseline Condition'); 
break; 
elseif strcmp(condition_input, 'offset') 
V = (v_offsetted_values + vc_off) * 10^-6; 
disp('Using Offset Condition'); 
break; 
else 
disp('Invalid input. Please enter either "Baseline" or "Offset".'); 
end 
end 
tolerance = 1; 
error = inf; 
P=zeros(1,720); 
T=zeros(1,720); 
mass_exh =zeros(1,180); 
mass_int =zeros(1,180); 
xb = zeros(1,55); 
%Intervals 
theta_start_int = 1;    theta_end_int = 227; 
theta_start_comp = 228; theta_end_comp = 330; 
theta_start_comb = 331; theta_end_comb = 390; 
theta_start_exp = 391;  theta_end_exp = 492; 
theta_start_exh = 493;  theta_end_exh = 720; 
P_intake = 101.325 * 10^3; T_intake  = 300; 
R = 287.05; %% gas constant of air (J/kg * Kelvin) 
AFR = 14.7; 
CBF = 0.95; 
Qhv = 45 * 10^6; 
Taf  = 303; 
% gamma values  
g_comp = 1.3; 
g_comb = 1.25; 
g_exp = 1.48; 
g_exh = 1.35; 
g_int = 1.32; 
dt = (1 / (rpm*2*pi / 60)) * ( pi / 180); 
Cd = 0.7; 
c_v= R/(g_comb-1); 
% xb constant values 
a = 5; 
m = 2 ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
PV Diagram     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%% 
while error > tolerance 
%% Compression 
T(theta_start_comp)=T_intake; 
P(theta_start_comp)=P_intake; 
m_i = (P_intake * V(theta_start_comp))/(R * T_intake); 
for i=theta_start_comp+1:theta_end_comp 
P(i)=(P_intake*(V(theta_start_comp)^g_comp))/(V(i)^g_comp); 
T(i)= (P(i)* V(i)) / (m_i * R); 
end 
%% Combustion 
m_f = m_i / (1 + AFR); 
T(theta_end_comb)=((CBF*m_f*Qhv)/(m_i*c_v))+T(theta_end_comp); 
P(theta_end_comb)=m_i*R*T(theta_end_comb)/V(theta_end_comb); 
P(theta_start_comb)= P(theta_end_comp); 
for i=theta_start_comb:theta_end_comb 
xb(i-theta_start_comb+1)= 1 - exp(-a * ((i - theta_start_comb) / (theta_end_comb - theta_start_comb))^(m + 1)); 
P(i) = (((P(theta_end_comb)*V(theta_end_comb)^(g_comb))
(P(theta_start_comb)*V(theta_start_comb)^(g_comb)))*xb(i-theta_start_comb+1) + 
(P(theta_start_comb)*V(theta_start_comb)^(g_comb))) / V(i)^(g_comb); 
end 
%% Expansion 
exp_cst = P(theta_end_comb) * V(theta_end_comb)^g_exp; 
for i= theta_start_exp:theta_end_exp 
P(i) = exp_cst / V(i)^g_exp;  
T(i)=P(i)*V(i)/(m_i*R); 
end 
%% Exhaust 
mass_exh(1) = m_i; 
P(theta_start_exh) =P(theta_end_exp);  
T(theta_start_exh) =T(theta_end_exp); 
pe = 101.325*10^3; 
for i= theta_start_exh+1:theta_end_exh 
mass_exh(i - theta_start_exh + 1) = mass_exh(i - theta_start_exh) - (dt * Cd * 
throat_area * P(i-1) / sqrt(P(i-1) * V(i-1) / mass_exh(i - theta_start_exh))) * ((pe / P(i
1))^(1 / g_exh)) * sqrt((2 * g_exh / (g_exh - 1)) * real(1 - (pe / P(i-1))^((g_exh - 1) / 
g_exh))); 
mass_exh(i - theta_start_exh + 1) = real(mass_exh(i - theta_start_exh + 1)); 
P(i) = (((T(i-1)*mass_exh(i-theta_start_exh+1)*R)/V(i))^g_exh)*(1/(P(i-1)^(g_exh-1))); 
T(i) = (P(i)*V(i)/(mass_exh(i-theta_start_exh+1)*R)); 
end 
%% Intake 
mass_end_exh = mass_exh(theta_end_exh - theta_start_exh + 1); 
mass_int(1) = mass_exh(theta_end_exh - theta_start_exh + 1); 
P(theta_start_int) = P(theta_end_exh);  
T(theta_start_int) = T(theta_end_exh); 
P_int = 101.325*10^3; 
for i= theta_start_int+1:theta_end_int 
mass_int(i - theta_start_int + 1) = mass_int(i - theta_start_int) + (dt * Cd * 
throat_area * P(i-1) / sqrt(R * T(theta_start_int))) * ((P(i-1) / P_int )^(1 / g_int)) * sqrt((2 
* g_int / (g_int - 1)) * (1 - (P(i-1) / P_int )^((g_int - 1) / g_int))); 
mass_int(i - theta_start_int + 1) =  real(mass_int(i - theta_start_int + 1)); 
P(i) = (((T(i-1)*mass_int(i-theta_start_int+1)*R) / V(i))^g_int)*(1/(P(i-1)^(g_int-1))); 
Xr = mass_end_exh/ mass_int(i - theta_start_int + 1); 
T(i) = Xr *T(theta_end_exh) + (1 - Xr) * Taf; 
end 
error = abs(T(223) - T_intake); 
T_intake = T_intake + error; 
disp(['Current error in temperature: ', num2str(error)]); 
end 
f
 igure; 
plot(V(1:720) * 10^6 , P(1:720) /101325, 'b-', 'LineWidth', 1.5); 
xlabel('Volume (cubic centimeters)'); 
ylabel('Pressure (atm)'); 
title('Pressure vs. Volume'); 
grid on; 
f
 igure; 
plot(theta , P(1:720) /101325, 'b-', 'LineWidth', 1.5); 
xlabel('theta (deg)'); 
ylabel('Pressure (atm)'); 
title('Pressure vs. theta'); 
grid on;