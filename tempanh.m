clear all; clc;
load('part1.mat')
load('part2.mat')

T0 = zeros(1,length(crank_angle));
a4y = -a4y;
alpha3 = -alpha3;

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

% Plot the torque for a single cylinder and the total torque
figure;
hold on;
plot(crank_angle, T0, 'r');
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
plot(crank_angle, T_total, 'r');
xlabel('Crank Angle (degrees)');
ylabel('Torque total (N-m)');
legend show;
title('Torque Plot');
grid on;
disp('doneee')

% % Plot coupler_angle0
% figure;
% plot(crank_angle, coupler_angle0);
% title('Coupler Angle0');
% xlabel('Crank Angle (degrees)');
% ylabel('Coupler Angle0');
% 
% % Plot a4y
% figure;
% plot(crank_angle, a4y);
% title('a4y');
% xlabel('Crank Angle (degrees)');
% ylabel('a4y (m/s^2)');
% 
% % Plot Fg
% figure;
% plot(crank_angle, Fg);
% title('Fg');
% xlabel('Crank Angle (degrees)');
% ylabel('Fg (N)');
% 
% % Plot alpha3
% figure;
% plot(crank_angle, alpha3);
% title('Alpha3');
% xlabel('Crank Angle (degrees)');
% ylabel('Alpha3 (rad/s^2)');
% 
% % Plot a3x
% figure;
% plot(crank_angle, a3x);
% title('a3x');
% xlabel('Crank Angle (degrees)');
% ylabel('a3x (m/s^2)');
% 
% % Plot a3y
% figure;
% plot(crank_angle, a3y);
% title('a3y');
% xlabel('Crank Angle (degrees)');
% ylabel('a3y (m/s^2)');

save('part3.mat')