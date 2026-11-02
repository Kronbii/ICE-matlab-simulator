function out = engine_cycle(c, rpm)
% Single-zone, open-system model of one cylinder over a full 720 degree cycle.
% The first law is integrated in crank angle for the gas in the cylinder:
%   dU = -p dV + dQ_comb - dQ_wall + sum(h dm) over both valves,
% integrated in steps of a tenth of a degree, each volume change taken as isentropic,
% with Wiebe heat release, Woschni heat transfer to the walls, and compressible
% flow (choked or subsonic, either way) through valves that lift and fall. Cycles
% repeat from the previous cycle's end state until the state at 720 degrees
% matches the one the cycle started from.
%
% c is the workspace of part1.m as a struct; rpm is the engine speed.
% Every per-degree output is sampled at crank angles 1..720.

  nsub = 10;                                  % integration steps per crank degree
  dth = 1/nsub;
  th = 0:dth:720;
  nstep = numel(th) - 1;
  r = c.crank_length; l = c.coupler_length; B = c.bore; A_p = c.A_cyl;
  Vd = c.vdisp0; Vc = Vd/(c.rc0 - 1);
  s = slider_crank(th, r, l, c.r0, 1, 0);
  V = Vc + A_p*(max(s) - s);                  % volume at every step (m^3)
  R = c.R; gu = c.yi; gb = c.yc;
  dt = dth/(6*rpm);                           % seconds per step
  Sp = 2*c.stroke0*rpm/60;                    % mean piston speed (m/s)

  % valve lift (m): a smooth bump over each valve's open window
  lift = @(x, a, b) c.lift_max*(x >= a & x <= b).*sin(pi*(x - a)/(b - a)).^2;
  ivo = c.intake_start - 1; ivc = c.intake_end;
  evo = c.exhaust_start; evc = c.exhaust_end;
  thm = th(1:end-1) + dth/2;                  % mid-step angles
  A_in = min(pi*c.d_valve*lift(thm, ivo, ivc), pi/4*c.d_valve^2);
  A_ex = min(pi*c.d_valve*lift(thm, evo, evc), pi/4*c.d_valve^2);
  open_any = (A_in > 0) | (A_ex > 0);

  % Wiebe burned fraction of the charge at every step
  th_s = c.combustion_start; dur = c.combustion_end - c.combustion_start + 1;
  xs = max(th - th_s, 0)/dur;
  xb_th = 1 - exp(-c.a*xs.^(c.m + 1));
  burning = thm >= th_s & thm < evo;          % Woschni's combustion term applies

  AF = c.AFR/c.phi;
  Pi = c.Pi; Pe = c.Pe; Ti = c.Taf; Tw = c.T_wall;
  cpu = gu*R/(gu - 1); cpb = gb*R/(gb - 1);

  % start of the first cycle: residual gas at exhaust pressure, TDC
  T = 1000; p = Pe; mass = p*V(1)/(R*T); mb = mass; T_exh = 900;
  U = mass*R/(gb - 1)*T;

  per = round(nsub);                          % steps per degree
  for cyc = 1:c.max_cycles
    start = [mass, U, mb];
    P = zeros(1, 720); TT = P; MM = P; FB = P; XB = P; MIN = P; MEX = P; QC = P; QW = P; HW = P;
    W_net = 0; W_gross = 0; Qc_tot = 0; Qw_tot = 0; m_fresh = 0; m_res = 0; m_fuel = 0;
    Qtot = 0; pr_ = p; Vr = V(1); Tr = T;
    hout = 0; mout = 0;
    for k = 1:nstep
      f = mb/mass;
      g = gu + (gb - gu)*f;
      cv = R/(g - 1);
      T = U/(mass*cv);
      p = mass*R*T/V(k);
      % reference state for Woschni and the fuel, at intake valve closing
      if abs(th(k) - ivc) < dth/2
        pr_ = p; Vr = V(k); Tr = T;
        m_res = mb; m_fresh = mass - mb;
        m_fuel = m_fresh/(1 + AF);
        Qtot = c.CBF*m_fuel*c.Qhv;
      end
      % flow through the intake valve (positive into the cylinder)
      dm_in = 0; dH_in = 0; dmb_in = 0;
      if A_in(k) > 0
        if Pi >= p
          q = c.Cd*A_in(k)*Pi/sqrt(R*Ti)*psi(p/Pi, gu);
          dm_in = q*dt; dH_in = cpu*Ti*dm_in;
        else
          q = c.Cd*A_in(k)*p/sqrt(R*T)*psi(Pi/p, g);
          dm_in = -q*dt; dH_in = g*cv*T*dm_in; dmb_in = f*dm_in;
        end
      end
      % flow through the exhaust valve (positive into the cylinder)
      dm_ex = 0; dH_ex = 0; dmb_ex = 0;
      if A_ex(k) > 0
        if p >= Pe
          q = c.Cd*A_ex(k)*p/sqrt(R*T)*psi(Pe/p, g);
          dm_ex = -q*dt; dH_ex = g*cv*T*dm_ex; dmb_ex = f*dm_ex;
          hout = hout + T*q*dt; mout = mout + q*dt;
        else
          q = c.Cd*A_ex(k)*Pe/sqrt(R*T_exh)*psi(p/Pe, gb);
          dm_ex = q*dt; dH_ex = cpb*T_exh*dm_ex; dmb_ex = dm_ex;
        end
      end
      % heat release and the unburned charge it turns into burned gas
      dxb = xb_th(k+1) - xb_th(k);
      dQc = Qtot*dxb;
      dmb_c = m_fresh*dxb;
      % Woschni heat transfer to the walls
      if open_any(k)
        C1 = 6.18; wgas = C1*Sp;
      else
        C1 = 2.28; wgas = C1*Sp;
        if burning(k)
          p_mot = pr_*(Vr/V(k))^gu;
          wgas = wgas + 3.24e-3*Vd*Tr/(pr_*Vr)*max(p - p_mot, 0);
        end
      end
      hw = 3.26*B^-0.2*(p/1000)^0.8*T^-0.55*wgas^0.8;
      Aw = 2*A_p + pi*B*V(k)/A_p;
      dQw = hw*Aw*(T - Tw)*dt;
      % first law: the volume change taken as isentropic over the step (exact for
      % the work term), then the heat and the flows through the valves
      U_comp = U*(V(k)/V(k+1))^(g - 1);
      dW = U - U_comp;                        % work done by the gas
      U = U_comp + dQc - dQw + dH_in + dH_ex;
      mass = mass + dm_in + dm_ex;
      mb = min(max(mb + dmb_in + dmb_ex + dmb_c, 0), mass);
      W_net = W_net + dW;
      if th(k) >= 180 && th(k) < 540
        W_gross = W_gross + dW;
      end
      Qc_tot = Qc_tot + dQc; Qw_tot = Qw_tot + dQw;
      % record per degree, at the end of each whole degree
      if mod(k, per) == 0
        i = k/per;
        f = mb/mass; g = gu + (gb - gu)*f; cv = R/(g - 1);
        P(i) = mass*R*(U/(mass*cv))/V(k+1); TT(i) = U/(mass*cv); MM(i) = mass; FB(i) = f;
        XB(i) = xb_th(k+1);
      end
      i = ceil(k/per);
      MIN(i) = MIN(i) + dm_in/dt/per; MEX(i) = MEX(i) + dm_ex/dt/per;
      QC(i) = QC(i) + dQc; QW(i) = QW(i) + dQw; HW(i) = HW(i) + hw/per;
    end
    if mout > 0
      T_exh = hout/mout;
    end
    change = max(abs([mass, U, mb] - start)./abs(start));
    if change < c.tolerance
      break
    end
  end

  out.pressure = P; out.temp = TT; out.mass = MM; out.burned = FB; out.xb = XB;
  out.mdot_in = MIN; out.mdot_ex = MEX; out.dQ_comb = QC; out.dQ_wall = QW; out.h_wall = HW;
  idx = (1:720)*per + 1;
  out.lift_in = lift(th(idx), ivo, ivc); out.lift_ex = lift(th(idx), evo, evc);
  out.W_net = W_net; out.W_gross = W_gross; out.W_pump = W_net - W_gross;
  out.Q_comb = Qc_tot; out.Q_wall = Qw_tot;
  out.m_fuel = m_fuel; out.m_fresh = m_fresh; out.m_air = m_fresh - m_fuel;
  out.residual = m_res/(m_res + m_fresh);
  out.T_exhaust = T_exh;
  out.cycles = cyc; out.change = change; out.rpm = rpm;
end

function y = psi(pr, g)
% Compressible flow function of an orifice: mass flow = Cd A p0/sqrt(R T0) psi,
% for the ratio pr of downstream to upstream pressure. Choked below the critical
% ratio; linear in the last 1% before the pressures equalise, which keeps the
% explicit integration steady as a valve flow dies out.
  crit = (2/(g + 1))^(g/(g - 1));
  if pr <= crit
    y = sqrt(g)*(2/(g + 1))^((g + 1)/(2*(g - 1)));
  elseif pr < 0.99
    y = pr^(1/g)*sqrt(2*g/(g - 1)*(1 - pr^((g - 1)/g)));
  else
    y0 = 0.99^(1/g)*sqrt(2*g/(g - 1)*(1 - 0.99^((g - 1)/g)));
    y = y0*(1 - pr)/0.01;
  end
end
