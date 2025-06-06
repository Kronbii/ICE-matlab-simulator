clear al; clc;

load('part2.mat')

coupler_length = coupler_length *100;
crank_length = crank_length * 100;
s_length0 = 100*s_length0;
c_length0 = 100*s_length0;

bore = 100*bore;

combustion_shift = combustion_shift - 5;

o1_x = 0; o1_y = 0;
o2_x = (2.2*crank_length); o2_y = 0;
o3_x = 2*(2.2*crank_length); o3_y = 0;
o4_x = 3*(2.2*crank_length); o4_y = 0;

figure;
vclear= vdisp0/(rc0 -1);
dclear = 400*vclear/(pi*(bore/100)^2);

plot_width = (o4_x - o1_x);
plot_offset = crank_length + 3;
plot_length = max(s_length0) + bore/4 +  dclear;

crank_angle_shifted= circshift(crank_angle, 180);
s_length_shifted = circshift(s_length0, 180);

for i=1:length(crank_angle)
set(gcf, 'Position', get(0, 'ScreenSize'));  % Maximizes the window to full screen

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

% Plotting the lines with color attributes
plot([o1_x crank_x1], [o1_y crank_y1], 'b', 'LineWidth', 1); % Crank 1 (red)
hold on;

% Plotting the joints
plot(o1_x, o1_y, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
plot(o2_x, o2_y, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
plot(o3_x, o3_y, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
plot(o4_x, o4_y, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')

plot(crank_x1, crank_y1, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
plot(crank_x2, crank_y2, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
plot(crank_x3, crank_y3, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
plot(crank_x4, crank_y4, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')

plot(coupler_x1, coupler_y1, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
plot(coupler_x2, coupler_y2, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
plot(coupler_x3, coupler_y3, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
plot(coupler_x4, coupler_y4, 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')

plot([crank_x1 coupler_x1], [crank_y1 coupler_y1], 'b', 'LineWidth', 1); % Coupler 1 (green)
plot([o2_x crank_x2], [o2_y crank_y2], 'b', 'LineWidth', 1); % Crank 2 (red)
plot([crank_x2 coupler_x2], [crank_y2 coupler_y2], 'b', 'LineWidth', 1); % Coupler 2 (green)
plot([o3_x crank_x3], [o3_y crank_y3], 'b', 'LineWidth', 1); % Crank 3 (red)
plot([crank_x3 coupler_x3], [crank_y3 coupler_y3], 'b', 'LineWidth', 1); % Coupler 3 (green)
plot([o4_x crank_x4], [o4_y crank_y4], 'b', 'LineWidth', 1); % Crank 4 (red)
plot([crank_x4 coupler_x4], [crank_y4 coupler_y4], 'b', 'LineWidth', 1); % Coupler 4 (green)

% Plotting vertical lines (for each mechanism)
plot([o1_x-bore/4 o1_x-bore/4], [min(s_length0)-bore/4 max(s_length0)+bore/4+dclear], 'k', 'LineWidth', 1);
plot([o2_x-bore/4 o2_x-bore/4], [min(s_length_shifted)-bore/4 max(s_length_shifted)+bore/4+dclear], 'k', 'LineWidth', 1);
plot([o3_x-bore/4 o3_x-bore/4], [min(s_length_shifted)-bore/4 max(s_length_shifted)+bore/4+dclear], 'k', 'LineWidth', 1);
plot([o4_x-bore/4 o4_x-bore/4], [min(s_length0)-bore/4 max(s_length0)+bore/4+dclear], 'k', 'LineWidth', 1);

% Plotting horizontal lines
plot([o1_x+bore/4 o1_x+bore/4], [min(s_length0)-bore/4 max(s_length0)+bore/4+dclear], 'k', 'LineWidth', 1);
plot([o2_x+bore/4 o2_x+bore/4], [min(s_length_shifted)-bore/4 max(s_length_shifted)+bore/4+dclear], 'k', 'LineWidth', 1);
plot([o3_x+bore/4 o3_x+bore/4], [min(s_length_shifted)-bore/4 max(s_length_shifted)+bore/4+dclear], 'k', 'LineWidth', 1);
plot([o4_x+bore/4 o4_x+bore/4], [min(s_length0)-bore/4 max(s_length0)+bore/4+dclear], 'k', 'LineWidth', 1);

% Plotting top lines
plot([o1_x-bore/4 o1_x+bore/4], [max(s_length0)+bore/4+dclear max(s_length0)+bore/4+dclear], 'k', 'LineWidth', 1);
plot([o2_x-bore/4 o2_x+bore/4], [max(s_length0)+bore/4+dclear max(s_length0)+bore/4+dclear], 'k', 'LineWidth', 1);
plot([o3_x-bore/4 o3_x+bore/4], [max(s_length0)+bore/4+dclear max(s_length0)+bore/4+dclear], 'k', 'LineWidth', 1);
plot([o4_x-bore/4 o4_x+bore/4], [max(s_length0)+bore/4+dclear max(s_length0)+bore/4+dclear], 'k', 'LineWidth', 1);

if crank_angle(i) >= (combustion_start-360 + combustion_shift) && crank_angle(i) <= (combustion_end-360-combustion_shift)
    % Draw the rectangle only if the angle is between alpha and beta
    rectangle('Position', [o1_x-bore/4, max(s_length0)+bore/4, bore/2, dclear], 'FaceColor', '#FF4500');
end

if crank_angle(i) >= (combustion_start-180 + combustion_shift) && crank_angle(i) <= (combustion_end-180-combustion_shift)
    % Draw the rectangle only if the angle is between alpha and beta
    rectangle('Position', [o3_x-bore/4, max(s_length0)+bore/4, bore/2, dclear], 'FaceColor', '#FF4500');
end

if crank_angle(i) >= (combustion_start + combustion_shift) && crank_angle(i) <= combustion_end-combustion_shift
    % Draw the rectangle only if the angle is between alpha and beta
    rectangle('Position', [o4_x-bore/4, max(s_length_shifted)+bore/4, bore/2, dclear], 'FaceColor', '#FF4500');
end

if crank_angle(i) >= (combustion_start+combustion_shift+180) && crank_angle(i) <= combustion_end+180-combustion_shift
    % Draw the rectangle only if the angle is between alpha and beta
    rectangle('Position', [o2_x-bore/4, max(s_length_shifted)+bore/4, bore/2, dclear], 'FaceColor', '#FF4500');
end

% Set axis equal to maintain aspect ratio
axis equal;
axis([o1_x-plot_offset, plot_width+plot_offset, o1_y-plot_offset, plot_length+plot_offset])
 
rectangle('Position', [coupler_x1-bore/4 coupler_y1-bore/4 bore/2 bore/2])
rectangle('Position', [coupler_x2-bore/4 coupler_y2-bore/4 bore/2 bore/2])
rectangle('Position', [coupler_x3-bore/4 coupler_y3-bore/4 bore/2 bore/2])
rectangle('Position', [coupler_x4-bore/4 coupler_y4-bore/4 bore/2 bore/2])

pause(0.000001);
hold off;

end