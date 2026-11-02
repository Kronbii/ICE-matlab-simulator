clear all; clc;
%% The four cylinders running, drawn from the results of part1.m and part2.m
load('part2.mat')

offset = [0 540 180 360];           % where each cylinder is in its cycle: firing 1-3-4-2
cycles = 2;                         % cycles to play
step = 2;                           % crank degrees per frame

% geometry in mm, with the engine upright: Y up the cylinder axis
r = 1e3*crank_length; Bh = 1e3*bore/2; pitch = 2.2*r;
pin_to_crown = 30;                  % piston pin to crown (mm)
clearance = 1e3*vdisp0/(rc0 - 1)/A_cyl;
deck = 1e3*max(s_length0) + pin_to_crown + clearance;
head = 30;
phase_name = {'intake', 'compression', 'combustion', 'expansion', 'exhaust'};
phase_end = [intake_end compression_end combustion_end expansion_end exhaust_end];
visible = {'off', 'on'};

fig = figure('Name', 'Four-cylinder engine', 'Color', 'w');
axes('Position', [0.03 0.08 0.94 0.84]); hold on; axis equal; axis off;
xlim([-Bh - 40, 3*pitch + Bh + 40]); ylim([-r - 60, deck + head + 20]);
for k = 1:4
  xc = (k - 1)*pitch;
  % cylinder walls and head
  plot([xc - Bh, xc - Bh, xc + Bh, xc + Bh], [1e3*min(s_length0) - 30, deck, deck, 1e3*min(s_length0) - 30], ...
       'k', 'LineWidth', 1.5);
  patch([xc - Bh - 6, xc + Bh + 6, xc + Bh + 6, xc - Bh - 6], [deck, deck, deck + head, deck + head], ...
        [0.75 0.75 0.78], 'EdgeColor', 'k');
  h(k).gas = patch(xc + [-Bh Bh Bh -Bh], [0 0 1 1], [0.2 0.3 0.6], 'EdgeColor', 'none');
  h(k).piston = patch(xc + [-Bh Bh Bh -Bh] * 0.98, [0 0 1 1], [0.6 0.62 0.66], 'EdgeColor', 'k');
  h(k).rod = plot([0 0], [0 0], 'Color', [0.35 0.35 0.4], 'LineWidth', 5);
  h(k).crank = plot([0 0], [0 0], 'Color', [0.2 0.2 0.25], 'LineWidth', 7);
  h(k).valve_in = plot([0 0], [0 0], 'b', 'LineWidth', 3);
  h(k).valve_ex = plot([0 0], [0 0], 'r', 'LineWidth', 3);
  h(k).spark = plot(xc, deck - 3, 'p', 'MarkerSize', 16, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', [1 0.5 0]);
  plot(xc, 0, 'ko', 'MarkerFaceColor', 'k');
  h(k).label = text(xc, -r - 35, '', 'HorizontalAlignment', 'center', 'FontSize', 11);
  text(xc, -r - 15, sprintf('Cylinder %d', k), 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end
ttl = title('', 'FontSize', 13);

for cyc_n = 1:cycles
  for i = 1:step:720
    if ~ishandle(fig)
      break
    end
    for k = 1:4
      j = mod(i - offset(k) - 1, 720) + 1;          % this cylinder's own crank angle
      xc = (k - 1)*pitch;
      s = 1e3*s_length0(j);
      pin = [xc - r*sind(j), r*cosd(j)];
      crown = s + pin_to_crown;
      t = min(max((temp(j) - 300)/2200, 0), 1);     % cool charge blue, burned gas red to white
      colour = [min(1, 0.15 + 1.6*t), min(1, max(0, 2*t - 0.6)), max(0, 0.6 - 0.9*t)];
      set(h(k).gas, 'YData', [crown; crown; deck; deck], 'FaceColor', colour);
      set(h(k).piston, 'YData', [s - 25; s - 25; crown; crown]);
      set(h(k).rod, 'XData', [pin(1) xc], 'YData', [pin(2) s]);
      set(h(k).crank, 'XData', [xc pin(1)], 'YData', [0 pin(2)]);
      li = 1e3*lift_in(j); le = 1e3*lift_ex(j);
      set(h(k).valve_in, 'XData', xc - 20 + [-12 12], 'YData', deck - [li li]);
      set(h(k).valve_ex, 'XData', xc + 20 + [-12 12], 'YData', deck - [le le]);
      set(h(k).spark, 'Visible', visible{1 + (abs(j - combustion_start) <= 4)});
      set(h(k).label, 'String', phase_name{find(j <= phase_end, 1)});
    end
    set(ttl, 'String', sprintf('Crank angle %3d deg   |   %d rpm   |   firing order 1-3-4-2', i, rpm));
    drawnow;
  end
end

