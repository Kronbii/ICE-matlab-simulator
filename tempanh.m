
load('part1.mat')
load('part2.mat')

Torque = zeros(1,length(crank_angle));

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
Torque(k) = -z(8);
end

plot(crank_angle, Torque)