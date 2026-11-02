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


def _angle_axis(ax):
    ax.set_xlim(0, 720)
    ax.set_xticks(range(0, 721, 90))


def _upto(n, progress):
    return max(1, int(round(n * min(max(progress, 0.0), 1.0))))


def _shade_phases(ax, d, th, labels=True, top=None):
    for name, a, b in phases(d):
        ax.axvspan(a - 1, b, color=th['phase'][name], alpha=th['phase_alpha'], lw=0, zorder=0)
        if labels and top is not None:
            ax.text((a - 1 + b) / 2, top, name, ha='center', va='top', fontsize=12,
                    color=th['muted'])


def _legend(ax, th):
    leg = ax.legend(frameon=False, fontsize=13, loc='best')
    for t in leg.get_texts():
        t.set_color(th['text'])


def piston_velocity(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    n = _upto(len(ang), progress)
    _style(ax, th, 'Piston velocity', 'crank angle (degrees)', 'velocity (m/s)')
    ax.plot(ang[:n], d['vel_piston0'][:n], color=th['line'], lw=2.6, label='no offset')
    ax.plot(ang[:n], d['vel_piston1'][:n], color=th['line2'], lw=2.6, label='18 mm offset')
    _angle_axis(ax)
    lim = 1.12 * max(np.abs(d['vel_piston0']).max(), np.abs(d['vel_piston1']).max())
    ax.set_ylim(-lim, lim)
    _legend(ax, th)


def piston_displacement(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    n = _upto(len(ang), progress)
    _style(ax, th, 'Piston distance from top dead centre', 'crank angle (degrees)', 'distance (cm)')
    ax.plot(ang[:n], 100 * d['c_length0'][:n], color=th['line'], lw=2.6, label='no offset')
    ax.plot(ang[:n], 100 * d['c_length1'][:n], color=th['line2'], lw=2.6, label='18 mm offset')
    _angle_axis(ax)
    ax.set_ylim(-0.4, 1.12 * 100 * max(d['c_length0'].max(), d['c_length1'].max()))
    _legend(ax, th)


def cylinder_volume(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    n = _upto(len(ang), progress)
    _style(ax, th, 'Cylinder volume', 'crank angle (degrees)', 'volume (cc)')
    ax.plot(ang[:n], 1e6 * d['total_volume0'][:n], color=th['line'], lw=2.6,
            label='no offset (compression ratio 9.8)')
    ax.plot(ang[:n], 1e6 * d['total_volume1'][:n], color=th['line2'], lw=2.6,
            label='18 mm offset (10.1)')
    _angle_axis(ax)
    hi = 1e6 * max(d['total_volume0'].max(), d['total_volume1'].max())
    ax.set_ylim(0, 1.15 * hi)
    _legend(ax, th)


def pressure_crank_angle(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    p = d['pressure'] / d['p_atm']
    n = _upto(len(ang), progress)
    _style(ax, th, 'Cylinder pressure over one four-stroke cycle', 'crank angle (degrees)',
           'cylinder pressure / atmospheric')
    top = 1.1 * p.max()
    _shade_phases(ax, d, th, top=top * 0.985)
    ax.plot(ang[:n], p[:n], color=th['line'], lw=2.8)
    _angle_axis(ax)
    ax.set_ylim(0, top)


def pv_diagram(ax, d, th, progress=1.0):
    v = 1e6 * d['vol']
    p = d['pressure'] / d['p_atm']
    n = _upto(len(v), progress)
    _style(ax, th, 'Pressure–volume diagram, one cylinder', 'cylinder volume (cc)',
           'cylinder pressure / atmospheric')
    for name, a, b in phases(d):
        lo, hi = a - 1, min(b, n)
        if hi <= lo:
            continue
        seg = slice(lo, min(hi + 1, len(v)))
        color = th['phase'][name] if th['phase_alpha'] < 1 else _PV_LIGHT[name]
        ax.plot(v[seg], p[seg], color=color, lw=2.8, label=name, solid_capstyle='round')
    ax.set_xlim(0.9 * v.min(), 1.05 * v.max())
    ax.set_ylim(0, 1.1 * p.max())
    if progress >= 1:
        _legend(ax, th)


_PV_LIGHT = {'intake': '#2f6fd0', 'compression': '#3c9a5f', 'combustion': '#e55e26',
             'expansion': '#c9821a', 'exhaust': '#7a7a7a'}


def torque_single(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    n = _upto(len(ang), progress)
    _style(ax, th, 'Crank torque of one cylinder', 'crank angle (degrees)', 'torque (N·m)')
    ax.axhline(0, color=th['axis'], lw=1)
    ax.plot(ang[:n], d['T0'][:n], color=th['line'], lw=2.6, label='cylinder 1')
    if progress >= 1:
        ax.axhline(d['T_avg'], color=th['line2'], lw=2, ls='--',
                   label=f"mean {float(d['T_avg']):.1f} N·m")
    _angle_axis(ax)
    ax.set_ylim(1.25 * d['T0'].min(), 1.12 * d['T0'].max())
    _legend(ax, th)


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
    ax.set_ylim(1.3 * lo, 1.18 * d['T_total'].max())
    leg = ax.legend(frameon=False, fontsize=12, loc='upper left', ncol=3)
    for t in leg.get_texts():
        t.set_color(th['text'])


def torque_flywheel(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    n = _upto(len(ang), progress)
    T0, Tm = d['T0'], float(d['T_avg'])
    _style(ax, th, 'Torque before and after adding the flywheel', 'crank angle (degrees)',
           'torque (N·m)')
    ax.axhline(0, color=th['axis'], lw=1)
    rng = np.asarray(d['integ_range'], dtype=int)
    shown = rng[rng <= n]
    if len(shown) > 1:
        idx = shown - 1
        ax.fill_between(ang[idx], T0[idx], Tm, where=T0[idx] >= Tm, color=th['line'],
                        alpha=0.22, lw=0, label=f"surplus energy {float(d['energy']):.0f} J")
    ax.plot(ang[:n], T0[:n], color=th['line'], lw=2.6, label='without flywheel')
    ax.plot(ang[:n], np.full(n, Tm), color=th['line2'], lw=2.6,
            label=f'with flywheel ({Tm:.1f} N·m)')
    _angle_axis(ax)
    ax.set_ylim(1.25 * T0.min(), 1.12 * T0.max())
    _legend(ax, th)


def speed_fluctuation(ax, d, th, progress=1.0):
    ang = d['crank_angle']
    n = _upto(len(ang), progress)
    rpm = lambda w: np.asarray(w) * 30 / np.pi
    _style(ax, th, 'Crank speed with the flywheel (4% fluctuation)', 'crank angle (degrees)',
           'crank speed (rpm)')
    ax.plot(ang[:n], rpm(d['w_fluc'])[:n], color=th['line'], lw=2.6, label='crank speed')
    ax.axhline(rpm(d['w2']), color=th['line2'], lw=2, ls='--', label='2000 rpm')
    for key in ('w_max', 'w_min'):
        ax.axhline(rpm(d[key]), color=th['muted'], lw=1.2, ls=':')
    _angle_axis(ax)
    lo, hi = rpm(d['w_min']), rpm(d['w_max'])
    pad = 0.35 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)
    _legend(ax, th)


def gear_ratios(ax, d, th, progress=1.0):
    iA = np.atleast_1d(d['iA']).astype(float)
    k = np.arange(1, len(iA) + 1)
    _style(ax, th, 'Eight-speed gearbox ratios', 'gear', 'ratio')
    grow = np.clip(progress * (len(iA) + 2) - np.arange(len(iA)), 0, 1)
    direct = int(d['direct_drive'])
    colors = [th['line2'] if g == direct else th['line'] for g in k]
    ax.bar(k, iA * grow, color=colors, width=0.62, zorder=2)
    for g, r, s in zip(k, iA, grow):
        if s >= 1:
            ax.text(g, r + 0.12, f'{r:.2f}', ha='center', va='bottom', fontsize=14,
                    color=th['text'])
    ax.set_xticks(k)
    ax.set_xticklabels([f'{g}{"st" if g == 1 else "nd" if g == 2 else "rd" if g == 3 else "th"}'
                        for g in k])
    ax.set_ylim(0, 1.15 * iA.max())
    ax.grid(axis='x', visible=False)
    if progress >= 1:
        ax.text(direct, iA[direct - 1] + 0.62, 'direct drive', ha='center', fontsize=12,
                color=th['line2'])


# name, title in the run video, the script that makes it, draw function
FIGURES = [
    ('piston-velocity', 'Piston velocity', 'part1.m', piston_velocity),
    ('piston-displacement', 'Piston displacement', 'part1.m', piston_displacement),
    ('cylinder-volume', 'Cylinder volume', 'part1.m', cylinder_volume),
    ('pressure-crank-angle', 'Pressure against crank angle', 'part2.m', pressure_crank_angle),
    ('pv-diagram', 'P–V diagram', 'part2.m', pv_diagram),
    ('torque-single-cylinder', 'Torque, one cylinder', 'part2.m', torque_single),
    ('torque-four-cylinder', 'Torque, four cylinders', 'part2.m', torque_four),
    ('torque-flywheel', 'Flywheel', 'part2.m', torque_flywheel),
    ('speed-fluctuation', 'Speed fluctuation', 'part2.m', speed_fluctuation),
    ('gear-ratios', 'Gear ratios, as printed', 'part2.m', gear_ratios),
]
