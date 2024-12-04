% Given values
r = 3.75; % crank length
s = 7.5; % slider length
theta = 1:180; % angle in radians

% Calculate coupler length
l = sqrt(r^2 + s^2 - 2 * r * s * cos(theta));

% Display the result
max(l)