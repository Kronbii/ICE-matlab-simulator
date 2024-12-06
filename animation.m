clc;

load('part2.mat')

coupler_length = coupler_length *100;
crank_length = crank_length * 100;
s_length0 = 100*s_length0;
bore = 100*bore;

o1_x = 0; o1_y = 0;

o2_x = 1.5*bore; o2_y = 0;

o3_x = 3*bore; o3_y = 0;

o4_x = 4.5*bore; o4_y = 0;

figure;

crank_angle_shifted= circshift(crank_angle, 180);
s_length_shifted = circshift(s_length0, 180);

for i=1:length(crank_angle)
    
crank_x1 = o1_x + crank_length * sind(crank_angle(i));
crank_x2 = o2_x + crank_length * sind(crank_angle_shifted(i));
crank_x3 = o3_x + crank_length * sind(crank_angle_shifted(i));
crank_x4 = o4_x + crank_length * sind(crank_angle(i));

crank_y1 = o1_y + crank_length * cosd(crank_angle(i));
crank_y2 = o2_y + crank_length * cosd(crank_angle_shifted(i));
crank_y3 = o3_y + crank_length * cosd(crank_angle_shifted(i));
crank_y4 = o4_y + crank_length * cosd(crank_angle(i));

coupler_x1 = o1_x;
coupler_x2 = o2_x;
coupler_x3 = o3_x;
coupler_x4 = o4_x;

coupler_y1 = o1_y + s_length0(i);
coupler_y2 = o2_y + s_length_shifted(i);
coupler_y3 = o3_y + s_length_shifted(i);
coupler_y4 = o4_y + s_length0(i);

plot([o1_x crank_x1], [o1_y crank_y1], [crank_x1 coupler_x1], [crank_y1 coupler_y1], ...
     [o2_x crank_x2], [o2_y crank_y2], [crank_x2 coupler_x2], [crank_y2 coupler_y2], ...
     [o3_x crank_x3], [o3_y crank_y3], [crank_x3 coupler_x3], [crank_y3 coupler_y3], ...
     [o4_x crank_x4], [o4_y crank_y4], [crank_x4 coupler_x4], [crank_y4 coupler_y4], ...
     [o1_x-bore/4 o1_x-bore/4],[min(s_length0)-bore/4 max(s_length0)+bore/4], ...
     [o2_x-bore/4 o2_x-bore/4],[min(s_length_shifted)-bore/4 max(s_length_shifted)+bore/4], ...
     [o3_x-bore/4 o3_x-bore/4],[min(s_length_shifted)-bore/4 max(s_length_shifted)+bore/4], ...
     [o4_x-bore/4 o4_x-bore/4],[min(s_length0)-bore/4 max(s_length0)+bore/4], ...
     [o1_x+bore/4 o1_x+bore/4],[min(s_length0)-bore/4 max(s_length0)+bore/4], ...
     [o2_x+bore/4 o2_x+bore/4],[min(s_length_shifted)-bore/4 max(s_length_shifted)+bore/4], ...
     [o3_x+bore/4 o3_x+bore/4],[min(s_length_shifted)-bore/4 max(s_length_shifted)+bore/4], ...
     [o4_x+bore/4 o4_x+bore/4],[min(s_length0)-bore/4 max(s_length0)+bore/4])
 
 
rectangle('Position', [coupler_x1-bore/4 coupler_y1-bore/4 bore/2 bore/2])
rectangle('Position', [coupler_x2-bore/4 coupler_y2-bore/4 bore/2 bore/2])
rectangle('Position', [coupler_x3-bore/4 coupler_y3-bore/4 bore/2 bore/2])
rectangle('Position', [coupler_x4-bore/4 coupler_y4-bore/4 bore/2 bore/2])

axis equal;
axis([-10 45 -10 25])
pause(0.0001);

end