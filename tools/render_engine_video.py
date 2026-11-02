"""Renders the engine running: four slider-cranks to scale beside live pressure, P–V and torque.

    python3 tools/render_engine_video.py [path/to/part2.mat]

Everything that moves comes from part2.mat: piston travel from s_length0, the gas column
from total_volume0, its colour from temp, the flame from the Wiebe burn of part2.m, the
valves from the phase timing in constants.m, and the charts from pressure, vol and T0..T3.
Writes docs/media/engine-simulation.mp4 (1080p60, three seamless cycles), a one-cycle GIF
for the README and a poster frame.
"""

import subprocess
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.patches import Circle, Polygon, Rectangle, Wedge

from simdata import (CYLINDER_COLOR, DARK, FIRING_ORDER, MEDIA, OFFSET, at, burned_fraction,
                     cylinder, load, phase_at, summary)

W, H, FPS = 1920, 1080, 60
STEP = 1.5            # crank degrees per frame
CYCLES = 3
START = 240.0         # engine crank angle of the first frame
TH = DARK

GAS = LinearSegmentedColormap.from_list('gas', [
    (0.00, '#1e3a5f'), (0.12, '#2f4d7d'), (0.28, '#7a3a63'),
    (0.50, '#d1452b'), (0.75, '#ff8c2a'), (1.00, '#ffe7a3')])
T_LO, T_HI = 300.0, 2850.0
PHASE_LABEL = {'intake': 'intake', 'compression': 'compression', 'combustion': 'combustion',
               'expansion': 'expansion', 'exhaust': 'exhaust'}
STEEL, STEEL_DARK, STEEL_EDGE = '#a7b0bd', '#6b7480', '#3d444f'
LIFT = 11.0          # mm, drawn valve lift
PORT = '#141920'


def valve_lift(d, phi, kind):
    """Lift in mm: a smooth bump over the window the cycle keeps that valve open."""
    x = phi % 720
    a, b = ((float(d['intake_start']) - 1, float(d['intake_end'])) if kind == 'in'
            else (float(d['exhaust_start']), float(d['exhaust_end'])))
    if not a <= x <= b:
        return 0.0
    return LIFT * np.sin(np.pi * (x - a) / (b - a)) ** 2


class Engine:
    def __init__(self, ax, d):
        self.d = d
        self.r = float(d['crank_length']) * 1000
        self.l = float(d['coupler_length']) * 1000
        self.bore = float(d['bore']) * 1000
        self.pitch = 2.2 * self.r                       # spacing animation.m uses
        self.area = np.pi / 4 * float(d['bore']) ** 2
        self.s_max = float(np.max(d['s_length0'])) * 1000
        self.s_min = float(np.min(d['s_length0'])) * 1000
        clearance = float(d['vdisp0']) / (float(d['rc0']) - 1) / self.area * 1000
        self.pin_to_crown = 34.0
        self.deck = self.s_max + self.pin_to_crown + clearance
        self.head_h = 40.0
        self.ax = ax
        self.cyl = {k: self._build(k) for k in (1, 2, 3, 4)}
        self._frame(ax)

    def _frame(self, ax):
        x0 = -self.bore / 2 - 70
        x1 = 3 * self.pitch + self.bore / 2 + 70
        ax.set_xlim(x0, x1)
        ax.set_ylim(-self.r - 120, self.deck + self.head_h + 42)
        ax.set_aspect('equal')
        ax.axis('off')

    def _build(self, k):
        ax, B = self.ax, self.bore
        xc = (k - 1) * self.pitch
        c = {}
        z = 1
        if k == 1:
            # block, head and crankshaft drawn once, behind the four cylinders
            left = -B / 2 - 16
            right = 3 * self.pitch + B / 2 + 16
            bottom = self.s_min - 26 - 8
            ax.add_patch(Rectangle((left, bottom), right - left, self.deck - bottom,
                                   fc='#20262e', ec='#363d47', lw=1.5, zorder=z))
            ax.add_patch(Rectangle((left - 6, self.deck), right - left + 12, self.head_h,
                                   fc='#2a313b', ec='#3f4752', lw=1.5, zorder=z + 6))
            ax.plot([left - 40, right + 40], [0, 0], color='#323944', lw=12,
                    solid_capstyle='round', zorder=z)
        # bore (empty below the piston) and the gas column above it
        ax.add_patch(Rectangle((xc - B / 2, self.s_min - 26 - 8), B,
                               self.deck - self.s_min + 34, fc=TH['bg'], ec='none', zorder=z + 1))
        c['gas'] = ax.add_patch(Rectangle((xc - B / 2, 0), B, 1, ec='none', zorder=z + 2))
        c['flame'] = ax.add_patch(Circle((xc, self.deck - 4), 1, fc='#ffd27a', ec='none',
                                         alpha=0, zorder=z + 3))
        c['flame'].set_clip_path(c['gas'])
        # crank web, counterweight, rod, piston, pins
        c['weight'] = ax.add_patch(Wedge((xc, 0), 40, 0, 1, fc='#4b535f', ec=STEEL_EDGE, lw=1,
                                         zorder=z + 2))
        c['web'] = ax.add_patch(Polygon(np.zeros((4, 2)), closed=True, fc=STEEL_DARK,
                                        ec=STEEL_EDGE, lw=1, zorder=z + 3))
        ax.add_patch(Circle((xc, 0), 13, fc='#2b3139', ec=STEEL_EDGE, lw=1.5, zorder=z + 4))
        c['rod'] = ax.add_patch(Polygon(np.zeros((4, 2)), closed=True, fc='#c3cad4',
                                        ec='#7d8693', lw=1, zorder=z + 4))
        c['bigend'] = ax.add_patch(Circle((0, 0), 11, fc='#c3cad4', ec='#7d8693', lw=1,
                                          zorder=z + 5))
        c['crankpin'] = ax.add_patch(Circle((0, 0), 6, fc='#59616d', ec='none', zorder=z + 6))
        c['piston'] = ax.add_patch(Rectangle((xc - B / 2 + 1.2, 0), B - 2.4, 60, fc=STEEL,
                                             ec='#d5dbe3', lw=1, zorder=z + 5))
        c['rings'] = [ax.plot([], [], color='#5f6774', lw=1.6, zorder=z + 6)[0] for _ in range(3)]
        c['pin'] = ax.add_patch(Circle((xc, 0), 7, fc='#6b7480', ec='#3d444f', lw=1,
                                       zorder=z + 6))
        # valves: a port through the head, a stem, and a face that drops into the chamber
        c['valves'] = {}
        top = self.deck + self.head_h
        for kind, dx in (('in', -19), ('ex', 19)):
            xv, out = xc + dx, (-1 if kind == 'in' else 1)
            port = ax.add_patch(Polygon([(xv - 13, self.deck), (xv + 13, self.deck),
                                         (xv + out * 15 + 10, top), (xv + out * 15 - 10, top)],
                                        closed=True, fc=PORT, ec='#3f4752', lw=1, zorder=z + 6.5))
            face = ax.add_patch(Polygon(np.zeros((4, 2)), closed=True, fc='#b9c0ca', ec='#e3e7ec',
                                        lw=0.8, zorder=z + 8))
            stem, = ax.plot([], [], color='#b9c0ca', lw=3.2, zorder=z + 7, solid_capstyle='butt')
            ends = [(xv + out * 13, top - 3), (xv + out * 3, self.deck + 7)]
            if kind == 'ex':
                ends.reverse()
            flow = ax.annotate('', xy=ends[1], xytext=ends[0], zorder=z + 9,
                               arrowprops=dict(arrowstyle='-|>', lw=2.4, mutation_scale=15,
                                               color='#ffffff'))
            tint = TH['line2'] if kind == 'in' else TH['line']
            c['valves'][kind] = (xv, face, stem, flow, port, tint)
        # spark plug
        ax.add_patch(Rectangle((xc - 4, self.deck + self.head_h), 8, 14, fc='#e9e4d8', ec='none',
                               zorder=z + 7))
        ax.add_patch(Rectangle((xc - 2.5, self.deck + self.head_h + 14), 5, 7, fc='#9aa3ae',
                               ec='none', zorder=z + 7))
        ax.plot([xc, xc], [self.deck - 3, self.deck + self.head_h], color='#9aa3ae', lw=2.2,
                zorder=z + 8)
        c['spark'], = ax.plot([xc], [self.deck - 5], marker='*', ms=30, color='#fff4b0',
                              mec='#ffd166', alpha=0, zorder=z + 10)
        c['glow'] = ax.add_patch(Circle((xc, self.deck - 5), 16, fc='#ffd166', ec='none',
                                        alpha=0, zorder=z + 9))
        c['glow'].set_clip_path(c['gas'])
        # labels under the crank
        y = -self.r - 58
        ax.text(xc, y, f'CYL {k}', ha='center', va='center', fontsize=17, fontweight='bold',
                color=CYLINDER_COLOR[k], zorder=z)
        c['phase'] = ax.text(xc, y - 24, '', ha='center', va='center', fontsize=14, zorder=z)
        c['p'] = ax.text(xc, y - 46, '', ha='center', va='center', fontsize=14,
                         color=TH['muted'], family='IBM Plex Mono', zorder=z)
        c['xc'] = xc
        return c

    def update(self, theta):
        d = self.d
        for k, c in self.cyl.items():
            phi = theta - OFFSET[k]
            xc = c['xc']
            s = float(at(d['s_length0'], phi)) * 1000
            ax_, ay_ = xc + self.r * np.sin(np.radians(phi)), self.r * np.cos(np.radians(phi))
            crown = s + self.pin_to_crown
            # gas and flame
            temp = float(at(d['temp'], phi))
            c['gas'].set_y(crown)
            c['gas'].set_height(self.deck - crown)
            c['gas'].set_facecolor(GAS(np.clip((temp - T_LO) / (T_HI - T_LO), 0, 1)))
            xb = burned_fraction(d, phi)
            burning = 0 < xb < 1
            c['flame'].set_radius(6 + 78 * np.sqrt(xb))
            c['flame'].set_alpha(0.75 * (1 - xb ** 2) if burning else 0)
            x = phi % 720
            ign = float(d['combustion_start'])
            spark = np.exp(-((x - ign) / 3.0) ** 2) if abs(x - ign) < 9 else 0
            c['spark'].set_alpha(spark)
            c['glow'].set_alpha(0.55 * spark)
            # moving parts
            alpha = np.degrees(np.arctan2(ay_, ax_ - xc))
            c['weight'].set_theta1(alpha + 180 - 62)
            c['weight'].set_theta2(alpha + 180 + 62)
            c['web'].set_xy(_bar((xc, 0), (ax_, ay_), 26, 18))
            c['rod'].set_xy(_bar((ax_, ay_), (xc, s), 16, 10))
            c['bigend'].set_center((ax_, ay_))
            c['crankpin'].set_center((ax_, ay_))
            c['piston'].set_y(s - 26)
            for j, ring in enumerate(c['rings']):
                yy = crown - 5 - 4.5 * j
                ring.set_data([xc - self.bore / 2 + 1.5, xc + self.bore / 2 - 1.5], [yy, yy])
            c['pin'].set_center((xc, s))
            for kind, (xv, face, stem, flow, port, tint) in c['valves'].items():
                lift = valve_lift(d, phi, kind)
                yb = self.deck - lift
                face.set_xy([(xv - 14, yb), (xv + 14, yb), (xv + 4, yb + 6), (xv - 4, yb + 6)])
                stem.set_data([xv, xv], [yb + 5, self.deck + self.head_h - 4])
                port.set_facecolor(_mix(PORT, tint, 0.75 * lift / LIFT))
                flow.set_visible(lift > 0.8)
                flow.arrow_patch.set_alpha(min(1, lift / 6))
            name = phase_at(d, phi)
            c['phase'].set_text(PHASE_LABEL[name])
            c['phase'].set_color(TH['phase'][name])
            c['p'].set_text(f"{float(at(d['pressure'], phi)) / float(d['p_atm']):5.1f} atm")


def _mix(a, b, t):
    """Colour a moved a fraction t of the way to colour b."""
    ca, cb = np.array(to_rgb(a)), np.array(to_rgb(b))
    return tuple(ca + (cb - ca) * float(np.clip(t, 0, 1)))


def _bar(p0, p1, w0, w1):
    """A tapered bar from p0 (width w0) to p1 (width w1) as a quad."""
    (x0, y0), (x1, y1) = p0, p1
    dx, dy = x1 - x0, y1 - y0
    n = np.hypot(dx, dy) or 1
    nx, ny = -dy / n, dx / n
    return [(x0 + nx * w0 / 2, y0 + ny * w0 / 2), (x1 + nx * w1 / 2, y1 + ny * w1 / 2),
            (x1 - nx * w1 / 2, y1 - ny * w1 / 2), (x0 - nx * w0 / 2, y0 - ny * w0 / 2)]


def _chart(fig, rect, title):
    ax = fig.add_axes(rect)
    ax.set_facecolor(TH['bg'])
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(TH['axis'])
    ax.tick_params(colors=TH['muted'], labelsize=12)
    ax.grid(True, color=TH['grid'], lw=1)
    ax.set_axisbelow(True)
    ax.set_title(title, loc='left', fontsize=16, color=TH['text'], pad=10)
    return ax


class Charts:
    def __init__(self, fig, d):
        self.d = d
        p_atm = float(d['p_atm'])
        grid = np.arange(0, 721, 1.0)
        self.p_atm = p_atm
        x0, w = 0.595, 0.375

        # pressure in each cylinder against engine crank angle
        ax = _chart(fig, [x0, 0.645, w, 0.185], 'Cylinder pressure  (atm)')
        self.p_curves = {}
        for k in (1, 2, 3, 4):
            ax.plot(grid, cylinder(d['pressure'], k, grid) / p_atm, color=CYLINDER_COLOR[k],
                    lw=1.8, alpha=0.9, label=f'cyl {k}')
        ax.set_xlim(0, 720)
        ax.set_xticks(range(0, 721, 180))
        ax.set_ylim(0, 1.12 * float(np.max(d['pressure'])) / p_atm)
        leg = ax.legend(frameon=False, fontsize=11, ncol=4, loc='upper right',
                        bbox_to_anchor=(1.0, 1.2), handlelength=1.2, columnspacing=1.0)
        for t in leg.get_texts():
            t.set_color(TH['muted'])
        self.p_cursor = ax.axvline(0, color='#ffffff', lw=1.2, alpha=0.5)
        self.p_dots = ax.scatter(np.zeros(4), np.zeros(4), s=55, zorder=5,
                                 c=[CYLINDER_COLOR[k] for k in (1, 2, 3, 4)],
                                 edgecolors='#ffffff', linewidths=1.2)

        # P–V loop of cylinder 1, the last half-turn traced bright
        ax = _chart(fig, [x0, 0.385, w, 0.185], 'Cylinder 1  ·  pressure against volume')
        v = 1e6 * np.asarray(d['vol'])
        p = np.asarray(d['pressure']) / p_atm
        ax.plot(np.append(v, v[0]), np.append(p, p[0]), color=TH['muted'], lw=1.2, alpha=0.45)
        self.trail, = ax.plot([], [], color=CYLINDER_COLOR[1], lw=3, solid_capstyle='round')
        self.pv_dot, = ax.plot([], [], 'o', ms=9, color=CYLINDER_COLOR[1], mec='#ffffff', mew=1.2)
        ax.set_xlim(0.85 * v.min(), 1.05 * v.max())
        ax.set_ylim(0, 1.12 * p.max())
        ax.set_xlabel('volume (cc)', fontsize=12, color=TH['muted'], labelpad=2)
        self.work = ax.text(0.98, 0.9, f"{float(d['engine_work']):.0f} J of work per cycle",
                            transform=ax.transAxes, ha='right', va='top', fontsize=13,
                            color=TH['muted'])

        # crank torque: each cylinder faint, the engine bold
        ax = _chart(fig, [x0, 0.125, w, 0.185], 'Crank torque  (N·m)')
        for k in (1, 2, 3, 4):
            ax.plot(grid, cylinder(d['T0'], k, grid), color=CYLINDER_COLOR[k], lw=1.2, alpha=0.45)
        self.total = cylinder(d['T0'], 1, grid) + cylinder(d['T0'], 2, grid) + \
            cylinder(d['T0'], 3, grid) + cylinder(d['T0'], 4, grid)
        ax.plot(grid, self.total, color=TH['text'], lw=2.4)
        mean = float(d['T_avg_total'])
        ax.axhline(mean, color=TH['line'], lw=1.6, ls='--')
        ax.text(0.995, 1.04, f'- - mean {mean:.0f} N·m', transform=ax.transAxes, ha='right',
                va='bottom', fontsize=12, color=TH['line'])
        ax.axhline(0, color=TH['axis'], lw=1)
        ax.set_xlim(0, 720)
        ax.set_xticks(range(0, 721, 180))
        ax.set_ylim(1.35 * self.total.min(), 1.2 * self.total.max())
        ax.set_xlabel('engine crank angle (degrees)', fontsize=12, color=TH['muted'], labelpad=2)
        self.t_cursor = ax.axvline(0, color='#ffffff', lw=1.2, alpha=0.5)
        self.t_dot, = ax.plot([], [], 'o', ms=9, color=TH['text'], mec=TH['bg'], mew=1.5)

    def update(self, theta):
        d, x = self.d, theta % 720
        self.p_cursor.set_xdata([x, x])
        self.p_dots.set_offsets([(x, float(cylinder(d['pressure'], k, theta)) / self.p_atm)
                                 for k in (1, 2, 3, 4)])
        back = np.arange(theta - 220, theta + 0.01, 1.0)
        self.trail.set_data(1e6 * at(d['vol'], back), at(d['pressure'], back) / self.p_atm)
        self.pv_dot.set_data([1e6 * float(at(d['vol'], theta))],
                             [float(at(d['pressure'], theta)) / self.p_atm])
        self.t_cursor.set_xdata([x, x])
        self.t_dot.set_data([x], [float(sum(cylinder(d['T0'], k, theta) for k in (1, 2, 3, 4)))])


class Header:
    def __init__(self, fig, d):
        s = summary(d)
        slow = (720 / STEP / FPS) / (120 / s['rpm'])
        fig.text(0.03, 0.945, 'Four-cylinder spark-ignition engine', fontsize=30,
                 fontweight='bold', color=TH['text'], va='center')
        fig.text(0.03, 0.902, f"MATLAB simulation  ·  {s['rpm']:.0f} rpm  ·  firing order "
                 f"{'-'.join(map(str, FIRING_ORDER))}  ·  shown {slow:.0f}× slower",
                 fontsize=16, color=TH['muted'], va='center')
        fig.text(0.97, 0.958, 'CRANK ANGLE', fontsize=12, color=TH['muted'], ha='right',
                 va='center')
        self.angle = fig.text(0.97, 0.912, '', fontsize=34, color=TH['text'], ha='right',
                              va='center', family='IBM Plex Mono')
        stats = [(f"{s['power_kw']:.1f} kW", 'indicated power'),
                 (f"{s['mean_torque']:.0f} N·m", 'mean torque'),
                 (f"{s['peak_pressure_atm']:.1f} atm", 'peak cylinder pressure'),
                 (f"{100 * s['efficiency']:.0f}%", 'fuel conversion efficiency'),
                 (f"{s['flywheel_cm']:.1f} cm", 'flywheel radius, 4% fluctuation')]
        for j, (big, small) in enumerate(stats):
            x = 0.03 + j * 0.19
            fig.text(x, 0.052, big, fontsize=24, color=TH['text'], fontweight='bold',
                     va='center')
            fig.text(x, 0.02, small, fontsize=13, color=TH['muted'], va='center')

    def update(self, theta):
        self.angle.set_text(f'{theta % 720:5.1f}°')


def render(d, path, frames, size=(W, H), fps=FPS, every=1, crf=18):
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100, facecolor=TH['bg'])
    engine = Engine(fig.add_axes([0.0, 0.1, 0.57, 0.77]), d)
    charts = Charts(fig, d)
    header = Header(fig, d)
    out = size
    cmd = ['ffmpeg', '-v', 'error', '-y', '-f', 'rawvideo', '-pix_fmt', 'rgba', '-s', f'{W}x{H}',
           '-r', str(fps), '-i', '-']
    if str(path).endswith('.gif'):
        cmd += ['-vf', f'scale={out[0]}:{out[1]}:flags=lanczos,split[a][b];'
                       '[a]palettegen=max_colors=192:stats_mode=diff[p];'
                       '[b][p]paletteuse=dither=bayer:bayer_scale=4:diff_mode=rectangle',
                '-loop', '0', str(path)]
    elif str(path).endswith(('.jpg', '.png', '.webp')):
        cmd += ['-frames:v', '1', '-vf', f'scale={out[0]}:{out[1]}:flags=lanczos', str(path)]
    else:
        cmd += ['-vf', f'scale={out[0]}:{out[1]}:flags=lanczos', '-c:v', 'libx264', '-preset',
                'slow', '-crf', str(crf), '-pix_fmt', 'yuv420p', '-movflags', '+faststart',
                str(path)]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    for i in range(0, frames, every):
        theta = START + i * STEP
        engine.update(theta)
        charts.update(theta)
        header.update(theta)
        fig.canvas.draw()
        proc.stdin.write(fig.canvas.buffer_rgba().tobytes())
        if i % (every * 120) == 0:
            print(f'  {path.name}: frame {i}/{frames}', flush=True)
    proc.stdin.close()
    proc.wait()
    plt.close(fig)


def main():
    plt.rcParams['font.family'] = 'IBM Plex Sans'
    d = load(sys.argv[1] if len(sys.argv) > 1 else None)
    per_cycle = int(round(720 / STEP))
    render(d, MEDIA / 'engine-simulation.mp4', CYCLES * per_cycle)
    render(d, MEDIA / 'engine-simulation.gif', per_cycle, size=(960, 540), fps=FPS // 3, every=3)
    global START
    START = 375.0           # cylinder 1 just after ignition
    render(d, MEDIA / 'engine-simulation-poster.jpg', 1)


if __name__ == '__main__':
    main()
