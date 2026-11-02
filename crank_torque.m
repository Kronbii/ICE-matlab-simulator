function T = crank_torque(c, Fg, m2, m3, m4, J3)
% Crank torque of one cylinder at every crank angle, from the force and moment
% balances of the crank, the connecting rod and the piston (an 8x8 linear system
% per angle). The engine stands upright: Y runs up the cylinder axis, gravity acts
% along -Y, the gas force Fg pushes the piston down and the crankcase pressure
% pushes it up. The unknowns are the pin forces F12, F23, F34, the wall force on
% the piston and the torque applied to the crank; the torque the cylinder delivers
% to the shaft is minus that. Pass zero masses for the torque of the gas alone.
  n = numel(c.crank_angle);
  T = zeros(1, n);
  F_atm = c.p_atm*c.A_cyl;
  r = c.crank_length; g = c.g;
  for k = 1:n
    th = c.crank_angle(k); ph = c.coupler_angle0(k);
    A = [0 0 0 0 1 0 1 0;
         0 0 0 0 0 1 0 0;
         0 0 1 0 -1 0 0 0;
         0 0 0 1 0 -1 0 0;
         0 0 c.AG3*cosd(ph) c.AG3*sind(ph) c.BG3*cosd(ph) c.BG3*sind(ph) 0 0;
         1 0 -1 0 0 0 0 0;
         0 1 0 -1 0 0 0 0;
         0 0 r*cosd(th) r*sind(th) 0 0 0 1];
    b = [0;
         m4*c.a4y(k) + Fg(k) - F_atm + m4*g;
         m3*c.a3x(k);
         m3*c.a3y(k) + m3*g;
         J3*c.alpha3(k);
         m2*c.a2x(k);
         m2*c.a2y(k) + m2*g;
         -m2*g*c.OG2*sind(th)];
    z = A\b;
    T(k) = -z(8);
  end
end
