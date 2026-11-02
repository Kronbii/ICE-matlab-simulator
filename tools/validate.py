"""Independent checks of the saved results, using none of the scripts' own formulas.

    python3 tools/validate.py [path/to/part2.mat]

- Mass: what flows in through the intake equals what leaves through the exhaust.
- Work: the cycle integrator's net work against the area of the sampled P–V loop.
- Gas torque: the 8x8 force solve against virtual work, -(p - p_atm) A ds/dtheta.
- Inertia torque: the force solve against the energy method, -d(KE + PE)/dt / omega,
  with every velocity taken by differentiating positions numerically.
"""

import sys

import numpy as np

from simdata import load


def main():
    d = load(sys.argv[1] if len(sys.argv) > 1 else None)
    f = lambda k: float(d[k])
    rpm, w = f('rpm'), f('w2')
    dt = 1 / (6 * rpm)                                      # seconds per crank degree
    m_in = np.sum(d['mdot_in']) * dt
    m_out = -np.sum(d['mdot_ex']) * dt
    print(f'mass per cycle: in {m_in * 1e6:.2f} mg, out {m_out * 1e6:.2f} mg '
          f'(difference {abs(m_in - m_out) / m_in:.1e})')

    V, p = np.asarray(d['vol']), np.asarray(d['pressure'])
    W_loop = np.trapz(np.append(p, p[0]), np.append(V, V[0]))
    print(f"net work: integrator {f('engine_work'):.2f} J, P-V loop {W_loop:.2f} J, "
          f"one cylinder's mean torque x 4 pi {f('T_avg') * 4 * np.pi:.2f} J")

    th = np.radians(1.0)
    diff = lambda x: (np.roll(x, -1) - np.roll(x, 1)) / (2 * th)
    s = np.asarray(d['s_length0'])
    T_vw = -(p - f('p_atm')) * f('A_cyl') * diff(s)
    gas = np.asarray(d['T0_gas'])
    print(f'gas torque: force solve vs virtual work, largest difference '
          f'{np.max(np.abs(T_vw - gas)):.3f} N m (peak {gas.max():.0f} N m)')

    t = np.radians(np.asarray(d['crank_angle'], dtype=float))
    ph = np.radians(np.asarray(d['coupler_angle0']))
    r, OG2, AG3 = f('crank_length'), f('OG2'), f('AG3')
    m2, m3, m4, J3, g = f('m2'), f('m3'), f('m4'), f('J3'), f('g')
    X2, Y2 = -OG2 * np.sin(t), OG2 * np.cos(t)              # engine upright: Y up the bore
    X3 = -r * np.sin(t) - AG3 * np.sin(ph)
    Y3 = r * np.cos(t) + AG3 * np.cos(ph)
    vel = lambda x: diff(x) * w
    KE = 0.5 * m2 * (vel(X2) ** 2 + vel(Y2) ** 2) + 0.5 * m3 * (vel(X3) ** 2 + vel(Y3) ** 2) \
        + 0.5 * J3 * vel(np.unwrap(ph)) ** 2 + 0.5 * m4 * vel(s) ** 2
    PE = g * (m2 * Y2 + m3 * Y3 + m4 * s)
    T_energy = -(vel(KE) + vel(PE)) / w
    inertia = np.asarray(d['T0_inertia'])
    print(f'inertia torque: force solve vs energy method, largest difference '
          f'{np.max(np.abs(T_energy - inertia)):.3f} N m (peak {np.max(np.abs(inertia)):.0f} N m)')


if __name__ == '__main__':
    main()
