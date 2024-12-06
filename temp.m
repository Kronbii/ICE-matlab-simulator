clear all; clc;
%%%%%%%% VOLUME %%%%%%%%%

%given: 
offset =0; % offset in m
a_radius = 50e-3; % crank radius in m
l = 150e-3; % connecting rod in m
theta = 1:720;
B = 80e-3; % bore in m
L =2*a_radius; % stroke in m
rc=9.8; % compression ratio
CBF=0.95; % coefficient of burned fuel
g_comp=1.3; % gamma for compression
g_exp=1.48; % gamma for expansion
g_comb=1.25; % gamma for combustion
g_exh=1.48; % gamma for exhaust
g_int=1.2; % gamma for intake
R=287;  % universal gas constant J/(kg.k)
T_Scomp = 300; % initial temperature at compression(K)
P_in = 0.9*101325; % starting pressure
P_atm=0.9*101325; % atmospheric pressure
rpm=2000; % N rev/min
cd=0.7; % discharge coefficient
Q_hv=44*10^6; % heating value (J/kg)
sAFR = 14.7; % stoichiometric air to fuel ratio
a1 = 5; % constant in Weibe function
m1 = 2; % constant in Weibe function
Phi = 0.98; % fuel/air equivalence ratio
Vc=(pi*(B^2)*L/4)/(rc-1); %clearance volume
Rf=l/a_radius; % Ratio of connecting rod length to crank radius
P=zeros(1,720); % pressure array 
T=zeros(1,720); % temperature array
m=zeros(1,720); % mass array
o2 = 1:720;
nr = 2; % for 4 stroke
air_density = 1.204; % kg/m^3
w2 = rpm*2*pi/60; % angular velocity
gravity = 9.81;
m3 = 0;%2.7215; % coupler mass
m2 = 0;%1.8143; % crank mass
m4 = 0;%0.765; % piston mass
F_atm = 0; % P_in*pi*(B^2/4);
OG2 = (1/3)*a_radius;
OA = a_radius;
AG2 = (2/3)*a_radius;
AG3 = (1/3)*l;
BG3 = (2/3)*l;
alpha2 =0;
ground_angle = 0;
J3 = (1/12)*m3*l^2 + m3*((l/2)-AG3)^2;
Pf=0; % friction power
K=0.01; % coefficient of fluctuation
density_steel=7800;


%%%%%%%%%%%%%%%%%%%%%%%%% Volume %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S is the distance between piston center and crank center
S = ((2*a_radius.*cosd(theta))+sqrt(((2*a_radius.*cosd(theta)).^2)-4*((a_radius^2)+2*offset*a_radius.*sind(theta)+(offset^2)-(l^2))))/2;
    
% b is the angle between l+a_radius and r1_max
b = asind(offset/(l+a_radius));
r1_max = (l+a_radius)*cosd(b);
    
% c is clearance 
c = r1_max - S;

% Vtotal
Vt_1 = (pi*B*B*c/4)+Vc;

% Flipping V:
V = flip(Vt_1);

%Finding theta at volume peaks

x1 = find (V == max(V));
BDC = [x1(1) x1(2)];
x2 = find (V == min(V)); 
TDC = [x2(1) x2(2)];

%%%%%%%%%%%%%%%%%% Valve timing angles %%%%%%%%%%%%

theta_Sint= TDC(1); %1;
theta_Eint=42 + BDC(1); %223;%180;
theta_Scomp=43 + BDC(1); %224;%181;
theta_Ecomp=TDC(2)-37; %330;%339;
theta_Scomb=TDC(2)-37; %331;%340;
theta_Ecomb=TDC(2)+28;%390;%380;
theta_Sexp=TDC(2)+29;%391;%381;
theta_Eexp=BDC(2)-52;%488;%540;
theta_Sexh=BDC(2)-53;%489;%541;
theta_Eexh= 720;%720;%720;

error = 10;
N=0;
P(theta_Scomp) = P_in;
T(theta_Scomp) = T_Scomp;
mi = (P(theta_Scomp).*(V(theta_Scomp)))/(R.*T(theta_Scomp)); %mass of air fuel mixture
m(theta_Scomp)=mi;
while error>5
%%%%%%%%%%%%%%%%%%%%% Compression Stroke%%%%%%%%%%%%%%%%%%%

%calculate compression stroke pressure and temperature%
for i=theta_Scomp+1:theta_Ecomp
P(i) = (P(i-1).*(V(i-1)^g_comp)./(V(i))^g_comp);
T(i)= (P(i).* (V(i)))/(mi*R);
m(i)=mi;
end 
 
%%%%%%%%%%%%%%%%%%%%% Combustion Stroke%%%%%%%%%%%%%%%%%%%
c_v = R/(g_comb-1);
T(theta_Scomb) = T(theta_Ecomp);
P(theta_Scomb) = P(theta_Ecomp);
m(theta_Scomb)=mi;

% Calculating P_Ecomb
m_f = mi /(1+(sAFR/Phi));
T(theta_Ecomb) = ((CBF*m_f*Q_hv)/(mi * c_v)) + T(theta_Scomb); 
P(theta_Ecomb) = (mi * R .* T(theta_Ecomb) )./ V(theta_Ecomb);

% Calculating P
for i = theta_Scomb:theta_Ecomb
    x_b(i) = 1 - exp(-a1*((theta(i) - theta_Scomb)/(theta_Ecomb - theta_Scomb)).^(m1+1));
    P(i) = ((x_b(i).*((P(theta_Ecomb)*(V(theta_Ecomb)^(g_comb)))- (P(theta_Scomb).*V(theta_Scomb).^(g_comb))) + (P(theta_Scomb).*(V(theta_Scomb).^(g_comb)))))./(V(i).^(g_comb));
    T(i)= (P(i).* (V(i)))./(mi*R);
    m(i)=mi;
end

%%%%%%%%%%%%%%%%%%%%% Expansion Stroke%%%%%%%%%%%%%%%%%%%
for i=theta_Sexp:theta_Eexp
P(i) = P(theta_Ecomb)*(V(theta_Ecomb)^g_exp)/(V(i)^g_exp);
T(i)=P(i)*V(i)/(mi*R);
m(i)=mi;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Exhaust %%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants and initial conditions
At = (pi/4)*(0.4*B)^2; % Throat area
m(theta_Sexh) = mi;
P(theta_Sexh)= (((T(theta_Eexp).*mi.*R)/V(theta_Sexh)).^(g_exh)).*(1/(P(theta_Eexp).^(g_exh-1)));
T(theta_Sexh)= (P(theta_Sexh).*V(theta_Sexh)./(mi.*R));

%calculation of dt
dt = (1/(rpm*2*pi/60))*(pi/180);

%calculation of mass theta
  for i = theta_Sexh+1: theta_Eexh
      m(i)=real(m(i-1)-(dt.*((cd*At.*P(i-1))./(sqrt((P(i-1).*V(i-1))./m(i-1)))).*((P_atm./P(i-1)).^(1/1.48)).*((((2*1.48)/(1.48-1)).*(1-(P_atm./P(i-1)).^((1.48-1)/1.48))).^0.5)));
      P(i) = ((((T(i-1))*(m(i))*R)/(V(i)))^g_exh) * (1/((P(i-1))^(g_exh-1)));
      T(i) = (P(i)*(V(i)))./(R*m(i));
      
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Intake %%%%%%%%%%%%%%%%%%%%%%%%%%
Ti = T_Scomp; % atmospheric temperature
m(theta_Sint) = m(theta_Eexh);
x(1)=m(theta_Eexh)/m(theta_Sint);
T(theta_Sint)=(x(1)*T(theta_Eexh))+((1-x(1))*Ti);
P(theta_Sint)=(m(theta_Sint)*R*T(theta_Sint))/(V(theta_Sint));


for i =theta_Sint+1:theta_Eint
    m(i)=real(m(i-1)+(dt.*((cd*At.*P(i-1))./(sqrt(R.*Ti))).*((P(i-1)./P_atm).^(1/1.3)).*((((2*1.3)/(1.3-1)).*(1-(P(i-1)./P_atm).^((1.3-1)/1.3))).^0.5)));
    x(i) = m(theta_Eexh)/m(i);
    T(i) = x(i)*T(theta_Eexh) + (1- (x(i)))*Ti;
    P(i) = (m(i)*R*T(i))/V(i);

end


 error=abs(T(theta_Eint)-T(theta_Scomp));
 T(theta_Scomp)=T(theta_Eint);  
 P(theta_Scomp)=P(theta_Eint);  
 mi=(P(theta_Scomp)*V(theta_Eint)/(R*T(theta_Scomp)));
 m(theta_Scomp)=mi;
 N=N+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%Finding torque%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Shifting angles by 180 degrees for each cylinder
o22 = circshift(o2,180);
o23 = circshift(o2,360);
o24 = circshift(o2,540);

P2 = circshift(P,180);
P3 = circshift(P,360);
P4 = circshift(P,540);

F_gas = P*pi*(B^2/4);
F_gas2 = P2*pi*(B^2/4);
F_gas3 = P3*pi*(B^2/4);
F_gas4 = P4*pi*(B^2/4);


%position analysis%
a = 1;
b = 2.*a_radius.*cosd(o2);
c = a_radius.^2 + 2.*offset.*a_radius.*sind(o2)+ offset.^2 - l.^2;

r1 = -b - sqrt(b.^2 -4.*a.*c)/(2.*a);
r1 = abs(r1);
r12 = circshift(r1,180);
r13 = circshift(r1,360);
r14 = circshift(r1,540);

o3 = atand((-offset-a_radius.*sind(o2))./(-a_radius.*cosd(o2) + r1));
o32 = circshift(o3,180);
o33 = circshift(o3,360);
o34 = circshift(o3,540);

%velocity analysis%
w3 = -(a_radius*w2.*cosd(o2))./(l.*cosd(o3));
w32 = circshift(w3,180);
w33 = circshift(w3,360);
w34 = circshift(w3,540);

vel = -a_radius*w2.*sind(o2) +a_radius*w2.*cosd(o2).*tand(o3);

%acceleration analysis%

a3 = (+w2^2 *a_radius.*sind(o2) -w3.^2 *l.*sind(o3))./(l.*cosd(o3));
a32 = circshift(a3,180);
a33 = circshift(a3,360);
a34 = circshift(a3,540);

acc = -w2.^2.*a_radius.*cosd(o2) -a3.*l.*sind(o3) +w3.^2 .*l.*cosd(o3);
acc2 = circshift(acc,180);
acc3 = circshift(acc,360);
acc4 = circshift(acc,540);

a2x = w2.^2 .*OG2 .*sind(o2);
a2x2 = circshift(a2x,180);
a2x3 = circshift(a2x,360);
a2x4 = circshift(a2x,540);

a2y = -w2.^2 .*OG2 .*cosd(o2);
a2y2 = circshift(a2y,180);
a2y3 = circshift(a2y,360);
a2y4 = circshift(a2y,540);

a3x = w2.^2.*OA.*sind(o2) + w3.^2 .*AG3.*sind(o3) -a3.*AG3.*cosd(o3);
a3x2 = circshift(a3x,180);
a3x3 = circshift(a3x,360);
a3x4 = circshift(a3x,540);

a3y = -w2.^2.*OA.*cosd(o2) - w3.^2 .*AG3.*cosd(o3) -a3.*AG3.*sind(o3);
a3y2 = circshift(a3y,180);
a3y3 = circshift(a3y,360);
a3y4 = circshift(a3y,540);

% solving for torque in cylinder 1
for k=1:length(o2) 
    
A1=[0 0 0 0 1 0 1 0; 
     0 0 0 0 0 1 0 0;
     0 0 1 0 -1 0 0 0;
     0 0 0 1 0 -1 0 0;
     0 0 AG3.*cosd(o3(k)) AG3.*sind(o3(k)) BG3.*cosd(o3(k)) BG3.*sind(o3(k)) 0 0;
     1 0 -1 0 0 0 0 0;
     0 1 0 -1 0 0 0 0;
     0 0 OA.*cosd(o2(k)) OA.*sind(o2(k)) 0 0 0 1];

B1=[0;
    m4.*acc(k)+F_gas(k)-F_atm+m4*gravity;
    m3.*a3x(k);
    m3.*a3y(k)+m3.*gravity;
    J3.*a3(k);
    m2.*a2x(k);
    m2.*a2y(k)+m2.*gravity;
    -(m2.*gravity).*OG2.*sind(o2(k))];
z = linsolve(A1,B1);
Torque(k) = -z(8);


% solving for torque in cylinder 2   
C1=[0 0 0 0 1 0 1 0; 
     0 0 0 0 0 1 0 0;
     0 0 1 0 -1 0 0 0;
     0 0 0 1 0 -1 0 0;
     0 0 AG3.*cosd(o32(k)) AG3.*sind(o32(k)) BG3.*cosd(o32(k)) BG3.*sind(o32(k)) 0 0;
     1 0 -1 0 0 0 0 0;
     0 1 0 -1 0 0 0 0;
     0 0 OA.*cosd(o22(k)) OA.*sind(o22(k)) 0 0 0 1];

D1=[0;
    m4.*acc2(k)+F_gas2(k)-F_atm+m4*gravity;
    m3.*a3x2(k);
    m3.*a3y2(k)+m3.*gravity;
    J3.*a32(k);
    m2.*a2x2(k);
    m2.*a2y2(k)+m2.*gravity;
    -(m2.*gravity).*OG2.*sind(o22(k))];
z2 = linsolve(C1,D1);
Torque2(k) = -z2(8);

% solving for torque in cylinder 3
E1  =[0 0 0 0 1 0 1 0; 
     0 0 0 0 0 1 0 0;
     0 0 1 0 -1 0 0 0;
     0 0 0 1 0 -1 0 0;
     0 0 AG3.*cosd(o33(k)) AG3.*sind(o33(k)) BG3.*cosd(o33(k)) BG3.*sind(o33(k)) 0 0;
     1 0 -1 0 0 0 0 0;
     0 1 0 -1 0 0 0 0;
     0 0 OA.*cosd(o23(k)) OA.*sind(o23(k)) 0 0 0 1];

F1 =[0;
    m4.*acc3(k)+F_gas3(k)-F_atm+m4*gravity;
    m3.*a3x3(k);
    m3.*a3y3(k)+m3.*gravity;
    J3.*a33(k);
    m2.*a2x3(k);
    m2.*a2y3(k)+m2.*gravity;
    -(m2.*gravity).*OG2.*sind(o23(k))];
z3 = linsolve(E1,F1);
Torque3(k) = -z3(8);

% solving for torque in cylinder 4

G1=[0 0 0 0 1 0 1 0; 
     0 0 0 0 0 1 0 0;
     0 0 1 0 -1 0 0 0;
     0 0 0 1 0 -1 0 0;
     0 0 AG3.*cosd(o34(k)) AG3.*sind(o34(k)) BG3.*cosd(o34(k)) BG3.*sind(o34(k)) 0 0;
     1 0 -1 0 0 0 0 0;
     0 1 0 -1 0 0 0 0;
     0 0 OA.*cosd(o24(k)) OA.*sind(o24(k)) 0 0 0 1];

H1=[0;
    m4.*acc4(k)+F_gas4(k)-F_atm+m4*gravity;
    m3.*a3x4(k);
    m3.*a3y4(k)+m3.*gravity;
    J3.*a34(k);
    m2.*a2x4(k);
    m2.*a2y4(k)+m2.*gravity;
    -(m2.*gravity).*OG2.*sind(o24(k))];
z4 = linsolve(G1,H1);
Torque4(k) = -z4(8);

end

Engine_torque = Torque + Torque2 + Torque3 + Torque4;

% Plot the torque for a single cylinder and the total torque
figure;
hold on;
plot(o2, Torque, 'r');
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
plot(o2, Engine_torque, 'r');
xlabel('Crank Angle (degrees)');
ylabel('Torque total (N-m)');
legend show;
title('Torque Plot');
grid on;

% Plot coupler_angle0
figure;
plot(o2, o3);
title('Coupler Angle0');
xlabel('Crank Angle (degrees)');
ylabel('Coupler Angle0');

% Plot a4y
figure;
plot(o2, acc4);
title('a4y');
xlabel('Crank Angle (degrees)');
ylabel('a4y (m/s^2)');

% Plot Fg
figure;
plot(o2, F_gas);
title('Fg');
xlabel('Crank Angle (degrees)');
ylabel('Fg (N)');

% Plot alpha3
figure;
plot(o2, a3);
title('Alpha3');
xlabel('Crank Angle (degrees)');
ylabel('Alpha3 (rad/s^2)');

% Plot a3x
figure;
plot(o2, a3x);
title('a3x');
xlabel('Crank Angle (degrees)');
ylabel('a3x (m/s^2)');

% Plot a3y
figure;
plot(o2, a3y);
title('a3y');
xlabel('Crank Angle (degrees)');
ylabel('a3y (m/s^2)');



%%%%%%%%%%%%%%%%%%%%%%% Work/Cycle %%%%%%%%%%%%%%%%%%%%%%%%%
Wi = trapz(V,P); % work/cycle

Vd = max(V)-Vc; % displaced volume per cylinder in m3

Pi = ((Wi*(rpm/60))/nr); % indicated power in W

mep = (Pi*nr)/(Vd*1000*rpm/60); % mean effective pressure in KPa

Pf=Pf*Pi; % Friction power

Mechanical_efficiency=(1-(Pf/Pi)); % mechanical efficiency

Pb=Mechanical_efficiency*Pi; % Brake power in W

bmep=(Pb*nr)/(Vd*1000*rpm/60); % brake mean effective pressure in KPa

thermo_torque=(bmep*Vd*1000)/(6.28*nr); % thermo torque in Nm

m_f_dot=(m_f*(rpm/60))/nr; % fuel mass flow rate in kg/s

Fuel_conv_efficiency=(Wi)/((m_f)*(Q_hv)); % Fuel conversion efficiency

sfc=(m_f_dot*10^-3)/(Pi*10^-3); % specific fuel consumption in mg/J

%sfc_3=(m_f_dot*10e-3*3600)/(Pi*10e-3); % specific fuel consumption in g/kwh

Volumetric_efficiency = (max(m)/(air_density*Vd)); % Volumetric efficiency

Sp_bar = 2*L*rpm/60; % piston speed in m/sec

Engine_power = (Pi*4)*(10^-3);% kW
Engine_horse_power = (Pi*4)*(10^-3)/(0.746);% kW
Total_torque = mean(Engine_torque);

Engine_volume =  ((max(V))*4)* 10^6; % cubic centimeter
Engine_volume_liters=Engine_volume/1000;


%%%%%%%%%%%%%%%%%%%%%%%% Flywheel Design %%%%%%%%%%%%%%%%%%%%%%%%%

% Average Torque
avg_Torque=mean(Torque);
avg_Torque=linspace(avg_Torque,avg_Torque,720);
avg_Torque_total= mean(Engine_torque);
avg_Torque_total=linspace(avg_Torque_total,avg_Torque_total,720);

%Calculating w_min and w_max from Torque
x3 = find(abs(diff(sign(Torque - avg_Torque))) > 0);

% Extract the x and y values at the points of intersection of Torque
% and average torque
x_intersect = theta(x3);
y_intersect=avg_Torque(x_intersect-1);

% Display the intersection point(s)
for i = 1:length(x3)
    if y_intersect(i)>Torque(x_intersect(i)-1)
        disp(['Intersection point theta_w_min is (' num2str(x_intersect(i)) ', ' num2str(y_intersect(i)) ')']);
        theta_w_min=x_intersect(i);
    else
        disp(['Intersection point theta_w_max is (' num2str(x_intersect(i)) ', ' num2str(y_intersect(i)) ')']);
        theta_w_max=x_intersect(i);
    end
    
end

%Calculating E
integral_Torque = trapz(theta(theta_w_min:theta_w_max), Torque(theta_w_min:theta_w_max));
integral_avg_Torque = trapz(theta(theta_w_min:theta_w_max), avg_Torque(theta_w_min:theta_w_max));

% Calculate E in J
E = 4*(integral_Torque - integral_avg_Torque);

%Calculating w_min and w_max
w_avg=(2*pi*rpm)/60;
syms w_min w_max
eqn1 = (w_max - w_min)/w_avg == K;
eqn2 = (w_max + w_min)/2 == w_avg;
sol = solve([eqn1, eqn2], [w_min, w_max]);
disp("w_min = " + string(vpa(sol.w_min, 6)));
disp("w_max = " + string(vpa(sol.w_max, 6)));

%Calculating I 
I=0.9*E/(K*((w_avg*(180/pi))^2)); %kg.m2

%Flywheel Design
%Our flywheel is a solid circular steel disk
%The cross-section of the rim is square
R_flywheel=30/w_avg; %m
m_flywheel=(2*I)/(R_flywheel^2); %kg
t_flywheel=sqrt(m_flywheel/(2*pi*R_flywheel*density_steel)); %m
b_flywheel=t_flywheel; %m

%%%%%%%%%%%%%%%%%%%%%%%% Transmission System %%%%%%%%%%%%%%%%%%%%%%%%%

mv=2320; %mass of vehicle in kg
Trans_torque= avg_Torque_total(1); %from code
fr=0.16;
lambda=1.3;
rdyn=0.392;
eff=0.885;
q=0.07;
nmax=4000;
n=8; 
d=6;
acceleration=4.54; %from code 
iE=3.636; 
vmax=245; 
iS=1;
g=9.81;

%Calculation of inclination angle alphast
alphast=atand(q);

%calculation of iamax
iamax= ((mv*g*(0.16*cosd(alphast)+sind(alphast)))+(mv*lambda*acceleration))/(Trans_torque*eff*1/rdyn);

%calculation of igmax
igmax=iamax/iE*iS;

%Calculation of iamin
iamin=3.6*pi/30*nmax*rdyn/vmax;

%calculation of igmin
igmin=iamin/iE*iS;

%Solve for phi1 and phi2 
eq1 = @(phi1,phi2) igmin*phi1.^(n-d).*phi2^(n-d-1) - 1;
eq2 = @(phi1,phi2) igmin*phi1.^(n-1).*phi2.^(0.5*(n-1)*(n-1-1)) - igmax;

% Use fsolve to solve the system of equations
sol = fsolve(@(phi1phi2) [eq1(phi1phi2(1),phi1phi2(2)); eq2(phi1phi2(1),phi1phi2(2))], [1;1]);

phi1=sol(1);
phi2=sol(2);

%gear ratios 
for x=1:n
    z=n;
    i(x)=(igmin*phi1^(z-x))*phi2^(0.5*(z-x)*(z-x-1)); 
end



% figure(1)
% plot(V,P),xlabel('Volume (m3)'), ylabel('Pressure (Pa)'), title('PV Diagram')
% figure(2)
% plot(o2,P),xlabel('theta (degree)'), ylabel('Cylinder Pressure (Pa)'), title('P - theta Diagram')
% figure(3)
% plot(o2,avg_Torque_total), xlabel('theta (degree)'), ylabel('Torque (Nm)'), title('Flywheel Output Torque vs. Crank Angle')
% figure(4)
% hold on
% plot(o2,Torque),xlabel('theta (degree)'), ylabel('Cylinder Torque (Nm)'), title('Cylinder Torque vs. Crank Angle')
% plot(o2,avg_Torque)
% figure(5)
% hold on
% plot(o2,Engine_torque),xlabel('theta (degree)'), ylabel('Total engine Torque (Nm)'), title('Engine Torque vs. Crank Angle')
% plot(o2,avg_Torque_total)
% figure(6)
% plot(o2, V), xlabel('theta (degree)'), ylabel('Cylinder Volume (m3)'), title('Volume vs. Crank Angle')
% 
