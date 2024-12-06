%set K =  0.01
% Is = E/K*wavg^2
% E = integral
% we need w given
% circular steel disk
% mass density = 7800 kg/m3 and 50 mm thick

%% Steps
% get mean torque
% get angles for which torque equal to avg torque
% get E from integrating T - Tavg from angle min to angle max
% calculate Is = E/K*wavg^2
% get radius from derived equation ( radius = 2*Is/(pi * thickness * density)

load('part3.mat')
