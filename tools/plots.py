"""Every figure the scripts produce, drawn from part2.mat.

Each draw function takes an axes, the loaded workspace, a theme from simdata and a
progress between 0 and 1: at 1 the figure is complete (the static images); below 1 the
curves are drawn only that far, which is how the run video draws them in.
"""

import numpy as np

from simdata import CYLINDER_COLOR, phases


def _style(ax, th, title, xlabel, ylabel):
    ax.set_facecolor(th['bg'])
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(th['axis'])
    ax.tick_params(colors=th['muted'], labelsize=13)
    ax.grid(True, color=th['grid'], linewidth=1)
    ax.set_axisbelow(True)
    ax.set_title(title, loc='left', fontsize=19, color=th['text'], pad=14)
    ax.set_xlabel(xlabel, fontsize=14, color=th['text'])
    ax.set_ylabel(ylabel, fontsize=14, color=th['text'])


def _twin(ax, th, label):
    ax2 = ax.twinx()
    ax2.set_ylabel(label, fontsize=13, color=th['muted'])
    ax2.tick_params(colors=th['muted'], labelsize=12)
    for side in ('top', 'left', 'bottom'):
        ax2.spines[side].set_visible(False)
    ax2.spines['right'].set_color(th['axis'])
    return ax2


def _angle_axis(ax, lo=0, hi=720, step=90):
    ax.set_xlim(lo, hi)
    ax.set_xticks(range(lo, hi + 1, step))


def _upto(n, progress):
    return max(1, int(round(n * min(max(progress, 0.0), 1.0))))


def _shade_phases(ax, d, th, top=None):
    for name, a, b in phases(d):
        ax.axvspan(a - 1, b, color=th['phase'][name], alpha=th['phase_alpha'], lw=0, zorder=0)
        if top is not None:
            ax.text((a - 1 + b) / 2, top, name, ha='center', va='top', fontsize=12,
                    color=th['muted'])


def _legend(ax, th, handles=None, **kw):
    opts = dict(frameon=False, fontsize=13, loc='best')
    opts.update(kw)
    leg = ax.legend(handles, [h.get_label() for h in handles], **opts) if handles \
        else ax.legend(**opts)
    for t in leg.get_texts():
        t.set_color(th['text'])


def _bar(d):
    return np.asarray(d['pressure']) / 1e5


def piston_velocity(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    n = _upto(len(ang), progress)
    _style(ax, th, 'Piston velocity', 'crank angle (degrees)', 'velocity (m/s)')
    ax.plot(ang[:n], d['vel_piston0'][:n], color=th['line'], lw=2.6, label='no offset')
    ax.plot(ang[:n], d['vel_piston1'][:n], color=th['line2'], lw=2.6, label='18 mm offset')
    _angle_axis(ax)
    lim = 1.25 * max(np.abs(d['vel_piston0']).max(), np.abs(d['vel_piston1']).max())
    ax.set_ylim(-lim, lim)
    _legend(ax, th, loc='upper right', ncol=2)


def piston_acceleration(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    n = _upto(len(ang), progress)
    _style(ax, th, 'Piston acceleration', 'crank angle (degrees)', 'acceleration (km/s²)')
    ax.axhline(0, color=th['axis'], lw=1)
    ax.plot(ang[:n], d['acc_piston0'][:n] / 1e3, color=th['line'], lw=2.6, label='no offset')
    ax.plot(ang[:n], d['acc_piston1'][:n] / 1e3, color=th['line2'], lw=2.6, label='18 mm offset')
    _angle_axis(ax)
    lo = min(d['acc_piston0'].min(), d['acc_piston1'].min()) / 1e3
    hi = max(d['acc_piston0'].max(), d['acc_piston1'].max()) / 1e3
    ax.set_ylim(1.15 * lo, 1.5 * hi)
    _legend(ax, th, loc='upper right', ncol=2)


def piston_displacement(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    n = _upto(len(ang), progress)
    _style(ax, th, 'Piston distance from top dead centre', 'crank angle (degrees)', 'distance (cm)')
    ax.plot(ang[:n], 100 * d['c_length0'][:n], color=th['line'], lw=2.6, label='no offset')
    ax.plot(ang[:n], 100 * d['c_length1'][:n], color=th['line2'], lw=2.6, label='18 mm offset')
    _angle_axis(ax)
    ax.set_ylim(-0.4, 1.22 * 100 * max(d['c_length0'].max(), d['c_length1'].max()))
    _legend(ax, th, loc='upper right', ncol=2)


def cylinder_volume(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    n = _upto(len(ang), progress)
    _style(ax, th, 'Cylinder volume', 'crank angle (degrees)', 'volume (cc)')
    ax.plot(ang[:n], 1e6 * d['total_volume0'][:n], color=th['line'], lw=2.6,
            label=f"no offset (compression ratio {float(d['rc0']):.1f})")
    ax.plot(ang[:n], 1e6 * d['total_volume1'][:n], color=th['line2'], lw=2.6,
            label=f"18 mm offset ({float(d['rc1']):.1f})")
    _angle_axis(ax)
    hi = 1e6 * max(d['total_volume0'].max(), d['total_volume1'].max())
    ax.set_ylim(0, 1.25 * hi)
    _legend(ax, th, loc='upper right', ncol=2)


def pressure_crank_angle(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    p = _bar(d)
    n = _upto(len(ang), progress)
    _style(ax, th, 'Cylinder pressure over one four-stroke cycle', 'crank angle (degrees)',
           'cylinder pressure (bar)')
    top = 1.1 * p.max()
    _shade_phases(ax, d, th, top=top * 0.985)
    ax.plot(ang[:n], p[:n], color=th['line'], lw=2.8)
    _angle_axis(ax)
    ax.set_ylim(0, top)


_PV_LIGHT = {'intake': '#2f6fd0', 'compression': '#3c9a5f', 'combustion': '#e55e26',
             'expansion': '#c9821a', 'exhaust': '#7a7a7a'}


def pv_diagram(ax, d, th, progress=1.0):
    v = 1e6 * np.asarray(d['vol'])
    p = _bar(d)
    n = _upto(len(v), progress)
    _style(ax, th, 'Pressure–volume diagram, one cylinder', 'cylinder volume (cc)',
           'cylinder pressure (bar)')
    for name, a, b in phases(d):
        lo, hi = a - 1, min(b, n)
        if hi <= lo:
            continue
        seg = slice(lo, min(hi + 1, len(v)))
        color = th['phase'][name] if th['phase_alpha'] < 1 else _PV_LIGHT[name]
        ax.plot(v[seg], p[seg], color=color, lw=2.8, label=name, solid_capstyle='round')
    ax.set_xlim(0, 1.05 * v.max())
    ax.set_ylim(0, 1.1 * p.max())
    if progress >= 1:
        ax.text(0.97, 0.6, f"net work {float(d['engine_work']):.0f} J per cycle",
                transform=ax.transAxes, ha='right', fontsize=14, color=th['muted'])
        _legend(ax, th, loc='upper right')


def temperature(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    n = _upto(len(ang), progress)
    T = np.asarray(d['temp'])
    _style(ax, th, 'Gas temperature in the cylinder', 'crank angle (degrees)', 'temperature (K)')
    _shade_phases(ax, d, th, top=1.12 * T.max() * 0.985)
    ax.plot(ang[:n], T[:n], color=th['line'], lw=2.8)
    ax.axhline(float(d['T_wall']), color=th['muted'], lw=1.2, ls=':')
    if progress >= 1:
        ax.text(8, float(d['T_wall']) + 40, f"wall {float(d['T_wall']):.0f} K", fontsize=12,
                color=th['muted'])
    _angle_axis(ax)
    ax.set_ylim(0, 1.12 * T.max())


def heat_release(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    lo, hi = 300, 540
    sel = (ang >= lo) & (ang <= hi)
    a = ang[sel]
    n = _upto(len(a), progress)
    _style(ax, th, 'Heat released by the burn, and lost to the walls', 'crank angle (degrees)',
           'joules per crank degree')
    qc, qw = np.asarray(d['dQ_comb'])[sel], np.asarray(d['dQ_wall'])[sel]
    ax.axhline(0, color=th['axis'], lw=1)
    ax.fill_between(a[:n], 0, qc[:n], color=th['line'], alpha=0.25, lw=0)
    ax.plot(a[:n], qc[:n], color=th['line'], lw=2.6,
            label=f"heat release (Wiebe, a = {float(d['a']):.0f}, m = {float(d['m']):.0f})")
    ax.plot(a[:n], qw[:n], color=th['line2'], lw=2.6, label='heat to the walls (Woschni)')
    ax.axvline(float(d['combustion_start']), color=th['muted'], lw=1.2, ls=':')
    if progress >= 1:
        ax.text(float(d['combustion_start']) - 3, 0.5 * qc.max(), 'spark', fontsize=12,
                color=th['muted'], ha='right')
    _angle_axis(ax, lo, hi, 30)
    ax.set_ylim(min(0, 1.2 * qw.min()), 1.3 * qc.max())
    _legend(ax, th, loc='upper right')


def gas_exchange(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    n = _upto(len(ang), progress)
    _style(ax, th, 'Gas exchange through the valves', 'crank angle (degrees)', 'mass flow (g/s)')
    mi, me = 1e3 * np.asarray(d['mdot_in']), -1e3 * np.asarray(d['mdot_ex'])
    ax.axhline(0, color=th['axis'], lw=1)
    ax2 = _twin(ax, th, 'valve lift (mm), shaded')
    ax2.fill_between(ang[:n], 0, 1e3 * np.asarray(d['lift_in'])[:n], color=th['line2'],
                     alpha=0.14, lw=0)
    ax2.fill_between(ang[:n], 0, 1e3 * np.asarray(d['lift_ex'])[:n], color=th['line'],
                     alpha=0.14, lw=0)
    ax2.set_ylim(0, 1e3 * float(d['lift_max']) * 2.6)
    ax.set_zorder(ax2.get_zorder() + 1)
    ax.patch.set_visible(False)
    l1, = ax.plot(ang[:n], mi[:n], color=th['line2'], lw=2.6, label='in through the intake valve')
    l2, = ax.plot(ang[:n], me[:n], color=th['line'], lw=2.6, label='out through the exhaust valve')
    _angle_axis(ax)
    ax.set_ylim(min(mi.min(), me.min()) * 1.4 - 5, max(mi.max(), me.max()) * 1.4)
    _legend(ax, th, [l1, l2], loc='upper center', ncol=2)


def torque_single(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    n = _upto(len(ang), progress)
    _style(ax, th, 'Crank torque of one cylinder', 'crank angle (degrees)', 'torque (N·m)')
    ax.axhline(0, color=th['axis'], lw=1)
    ax.plot(ang[:n], d['T0_gas'][:n], color=th['line'], lw=1.8, ls='--', label='gas')
    ax.plot(ang[:n], d['T0_inertia'][:n], color=th['line2'], lw=1.8, ls='--',
            label='inertia and gravity')
    ax.plot(ang[:n], d['T0'][:n], color=th['text'], lw=2.6, label='total')
    _angle_axis(ax)
    ax.set_ylim(1.5 * min(d['T0'].min(), d['T0_inertia'].min()), 1.15 * d['T0_gas'].max())
    _legend(ax, th, loc='upper right')


def torque_four(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    n = _upto(len(ang), progress)
    _style(ax, th, 'Crank torque of four cylinders, 180° apart', 'crank angle (degrees)',
           'torque (N·m)')
    ax.axhline(0, color=th['axis'], lw=1)
    for k, key in zip((1, 3, 4, 2), ('T0', 'T1', 'T2', 'T3')):
        ax.plot(ang[:n], d[key][:n], color=CYLINDER_COLOR[k], lw=1.6, alpha=th['faint'],
                label=f'cylinder {k}')
    ax.plot(ang[:n], d['T_total'][:n], color=th['text'], lw=2.8, label='engine (sum)')
    if progress >= 1:
        ax.axhline(d['T_avg_total'], color=th['line'], lw=2, ls='--',
                   label=f"mean {float(d['T_avg_total']):.0f} N·m")
    _angle_axis(ax)
    lo = min(d['T_total'].min(), d['T0'].min())
    ax.set_ylim(1.3 * lo, 1.35 * d['T_total'].max())
    _legend(ax, th, loc='upper left', ncol=3, fontsize=12)


def flywheel_energy(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    n = _upto(len(ang), progress)
    E = np.asarray(d['E_fly'])
    _style(ax, th, 'Shaft energy above the mean load', 'crank angle (degrees)', 'energy (J)')
    ax.fill_between(ang[:n], E.min(), E[:n], color=th['line'], alpha=0.12, lw=0)
    ax.plot(ang[:n], E[:n], color=th['line'], lw=2.6)
    pad = 0.18 * (E.max() - E.min())
    if progress >= 1:
        i = int(np.argmax(E))
        x = ang[i] + 30
        ax.annotate('', xy=(x, E.max()), xytext=(x, E.min()),
                    arrowprops=dict(arrowstyle='<->', color=th['line2'], lw=2))
        ax.text(x + 8, (E.max() + E.min()) / 2,
                f"swing {float(d['energy']):.0f} J\nflywheel radius "
                f"{100 * float(d['flywheel_radius']):.1f} cm",
                fontsize=13, color=th['line2'], va='center')
    _angle_axis(ax)
    ax.set_ylim(E.min() - pad, E.max() + pad)


def speed_fluctuation(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    n = _upto(len(ang), progress)
    rpm = np.asarray(d['w_fluc']) * 30 / np.pi
    _style(ax, th, 'Crank speed with the flywheel', 'crank angle (degrees)', 'crank speed (rpm)')
    ax.plot(ang[:n], rpm[:n], color=th['line'], lw=2.6, label='crank speed')
    ax.axhline(float(d['rpm']), color=th['line2'], lw=2, ls='--',
               label=f"{float(d['rpm']):.0f} rpm")
    _angle_axis(ax)
    lo, hi = rpm.min(), rpm.max()
    pad = 0.5 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)
    _legend(ax, th, loc='upper right', ncol=2)


def torque_power(ax, d, th, progress=1.0):
    n_all = np.asarray(d['rpm_sweep'], dtype=float)
    tq, pw = np.asarray(d['sweep_torque']), np.asarray(d['sweep_power']) / 1e3
    k = max(1, int(np.ceil(len(n_all) * min(max(progress, 0.0), 1.0))))
    _style(ax, th, 'Full-load torque and power', 'engine speed (rpm)', 'brake torque (N·m)')
    l1, = ax.plot(n_all[:k], tq[:k], 'o-', color=th['line'], lw=2.8, ms=7, label='brake torque')
    ax.set_xlim(n_all[0] - 250, n_all[-1] + 250)
    ax.set_ylim(0, 1.3 * tq.max())
    ax2 = _twin(ax, th, 'brake power (kW)')
    l2, = ax2.plot(n_all[:k], pw[:k], 's-', color=th['line2'], lw=2.8, ms=7, label='brake power')
    ax2.set_ylim(0, 1.3 * pw.max())
    if progress >= 1:
        i, j = int(np.argmax(tq)), int(np.argmax(pw))
        ax.annotate(f'{tq[i]:.0f} N·m at {n_all[i]:.0f} rpm', (n_all[i], tq[i]),
                    xytext=(0, 14), textcoords='offset points', ha='center', fontsize=13,
                    color=th['line'])
        ax2.annotate(f'{pw[j]:.1f} kW at {n_all[j]:.0f} rpm', (n_all[j], pw[j]),
                     xytext=(0, 14), textcoords='offset points', ha='center', fontsize=13,
                     color=th['line2'])
    _legend(ax, th, [l1, l2], loc='lower left')


def gear_ratios(ax, d, th, progress=1.0):
    iA = np.atleast_1d(d['iA']).astype(float)
    k = np.arange(1, len(iA) + 1)
    _style(ax, th, 'Eight-speed gearbox ratios', '', 'ratio')
    grow = np.clip(progress * (len(iA) + 2) - np.arange(len(iA)), 0, 1)
    direct = int(d['direct_drive'])
    colors = [th['line2'] if g == direct else th['line'] for g in k]
    ax.bar(k, iA * grow, color=colors, width=0.62, zorder=2)
    v = np.atleast_1d(d['v_gear']).astype(float)
    for g, r, s, sp in zip(k, iA, grow, v):
        if s >= 1:
            ax.text(g, r + 0.12, f'{r:.2f}', ha='center', va='bottom', fontsize=14,
                    color=th['text'])
    ax.set_xticks(k)
    ax.set_xticklabels([f'{g}{"st" if g == 1 else "nd" if g == 2 else "rd" if g == 3 else "th"}'
                        f'\n{sp:.0f} km/h' for g, sp in zip(k, v)])
    ax.set_ylim(0, 1.18 * iA.max())
    ax.grid(axis='x', visible=False)
    if progress >= 1:
        ax.text(direct, iA[direct - 1] + 0.62, 'direct drive', ha='center', fontsize=12,
                color=th['line2'])
        ax.text(0.99, 0.95, f"gear, and road speed at {float(d['rpm_max']):.0f} rpm",
                transform=ax.transAxes, ha='right', fontsize=12, color=th['muted'])


# name, title in the run video, the script that makes it, draw function
FIGURES = [
    ('piston-velocity', 'Piston velocity', 'part1.m', piston_velocity),
    ('piston-acceleration', 'Piston acceleration', 'part1.m', piston_acceleration),
    ('piston-displacement', 'Piston displacement', 'part1.m', piston_displacement),
    ('cylinder-volume', 'Cylinder volume', 'part1.m', cylinder_volume),
    ('pressure-crank-angle', 'Pressure against crank angle', 'part2.m', pressure_crank_angle),
    ('pv-diagram', 'P–V diagram', 'part2.m', pv_diagram),
    ('temperature', 'Gas temperature', 'part2.m', temperature),
    ('heat-release', 'Heat release', 'part2.m', heat_release),
    ('gas-exchange', 'Gas exchange', 'part2.m', gas_exchange),
    ('torque-single-cylinder', 'Torque, one cylinder', 'part2.m', torque_single),
    ('torque-four-cylinder', 'Torque, four cylinders', 'part2.m', torque_four),
    ('flywheel-energy', 'Flywheel', 'part2.m', flywheel_energy),
    ('speed-fluctuation', 'Crank speed', 'part2.m', speed_fluctuation),
    ('torque-power-curve', 'Torque and power', 'part2.m', torque_power),
    ('gear-ratios', 'Gear ratios', 'part2.m', gear_ratios),
]
