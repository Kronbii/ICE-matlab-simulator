"""Replays a run of the simulator: its console output on the left, its figures on the right.

    python3 tools/render_run_video.py [path/to/part2.mat] [path/to/console.log]

The console is the log of a real run (docs/octave-run.log, from tools/octave/run_all.m),
replayed line by line; each figure is drawn from the part2.mat that same run saved, in the
order the scripts open them. Writes docs/media/simulation-run.mp4 (1080p30) and a poster.
"""

import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

from plots import FIGURES
from simdata import DARK, MEDIA, ROOT, load, summary

W, H, FPS = 1920, 1080, 30
TH = DARK
INTRO, DRAW, HOLD, OUTRO = 2.6, 1.3, 1.8, 6.0
CONSOLE_LINES = 36
CONSOLE_WIDTH = 55   # characters that fit the Command Window
SCRIPTS = ('constants', 'part1', 'part2')


def parse(log):
    """The log split into (script, output lines) in the order they ran."""
    runs, cur = [], None
    for line in Path(log).read_text().splitlines():
        if line.startswith('>> '):
            cur = (line[3:].strip(), [])
            if cur[0] in SCRIPTS:           # the animation has its own video
                runs.append(cur)
        elif cur is not None:
            cur[1].extend(wrap(line.rstrip()))
    while runs and runs[-1][1] and not runs[-1][1][-1]:
        runs[-1][1].pop()
    return runs


def wrap(line):
    """A line too wide for the Command Window, broken at a space as a narrow terminal would."""
    out = []
    while len(line) > CONSOLE_WIDTH:
        cut = line.rfind(' ', 0, CONSOLE_WIDTH)
        if cut <= len(line) - len(line.lstrip()):
            break
        out.append(line[:cut].rstrip())
        line = '  ' + line[cut:].lstrip()
    out.append(line)
    return out


def schedule(runs):
    """Times for every console line, every figure and the closing card, in the order a
    script prints its lines, opens its figures and says it is done."""
    lines, typing, figs = [], [], []
    by_script = {}
    for i, (_, _, script, _) in enumerate(FIGURES):
        by_script.setdefault(script[:-2], []).append(i)
    t = INTRO + 0.4
    for script, out in runs:
        dur = 0.35 + 0.07 * len(script)
        typing.append((t, t + dur, script))
        t += dur
        lines.append((t, f'>> {script}'))
        t += 0.35
        body, tail = (out[:-1], out[-1:]) if out and out[-1] == 'done' else (out, [])
        for line in body:
            lines.append((t, line))
            sweep_row = line[:6].strip().isdigit()
            t += 0.32 if sweep_row else 0.05 if not line else 0.13
        if body:
            t += 0.5
        for i in by_script.get(script, []):
            figs.append((t, i))
            t += DRAW + HOLD
        for line in tail:
            lines.append((t, line))
            t += 0.4
    return lines, typing, figs, t, t + OUTRO


def ease(x):
    x = min(max(x, 0.0), 1.0)
    return 1 - (1 - x) ** 3


def line_color(text):
    if text.startswith('>>'):
        return TH['text']
    if text.startswith(('Kinematics check', 'Cycle converged', 'Torque check')) or text == 'done':
        return TH['line3']
    if text.startswith(('Brake torque', 'Indicated power', 'Peak torque', 'Gear ratios',
                        'Compression ratio', 'Net work')):
        return TH['line']
    if text[:6].strip().isdigit() or text.lstrip().startswith('rpm'):
        return TH['text']
    return '#c9d1d9'


def panel(fig, rect, title):
    x, y, w, h = rect
    fig.add_artist(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=0.008',
                                  transform=fig.transFigure, fc=TH['panel'], ec='#30363d',
                                  lw=1.2, zorder=0))
    fig.text(x + 0.012, y + h - 0.03, title, fontsize=15, color=TH['muted'], va='center')


def render(d, runs, path, only=None):
    lines, typing, figs, t_card, t_end = schedule(runs)
    s = summary(d)
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100, facecolor=TH['bg'])

    # top bar
    fig.text(0.017, 0.958, 'ICE-matlab-simulator', fontsize=22, fontweight='bold',
             color=TH['text'], va='center')
    fig.text(0.017, 0.925, 'The MATLAB scripts run in GNU Octave: every line they print and '
             'every figure they draw', fontsize=14, color=TH['muted'],
             va='center')
    chips = {}
    for j, name in enumerate(SCRIPTS):
        x = 0.70 + j * 0.095
        box = FancyBboxPatch((x, 0.935), 0.085, 0.036, boxstyle='round,pad=0,rounding_size=0.01',
                             transform=fig.transFigure, fc=TH['bg'], ec='#30363d', lw=1.4)
        fig.add_artist(box)
        label = fig.text(x + 0.0425, 0.953, f'{name}.m', fontsize=13, ha='center', va='center',
                         color=TH['muted'], family='IBM Plex Mono')
        chips[name] = (box, label)

    # command window
    panel(fig, (0.017, 0.035, 0.33, 0.855), 'Command Window')
    rows = [fig.text(0.028, 0.835 - k * 0.0222, '', fontsize=13.2, family='IBM Plex Mono',
                     va='center', color=TH['text']) for k in range(CONSOLE_LINES)]

    # figure window
    panel(fig, (0.362, 0.035, 0.621, 0.855), '')
    fig_title = fig.text(0.374, 0.86, '', fontsize=15, color=TH['muted'], va='center')
    nf = len(FIGURES)
    dots = [fig.text(0.968 - (nf - 1 - k) * 0.0135, 0.86, '●', fontsize=11, color='#30363d',
                     ha='center', va='center', family='DejaVu Sans') for k in range(nf)]
    ax = fig.add_axes([0.41, 0.1, 0.55, 0.69])
    own_axes = list(fig.axes)
    veil = fig.add_artist(Rectangle((0.364, 0.037), 0.617, 0.81, transform=fig.transFigure,
                                    fc=TH['panel'], ec='none', alpha=0, zorder=5))
    waiting = fig.text(0.672, 0.44, '', fontsize=17, color=TH['muted'], ha='center',
                       va='center', zorder=6)

    # closing card, laid over the figure window
    card = []
    big = [(f"{s['torque_peak']:.0f} N·m", f"peak brake torque, {s['rpm_torque_peak']:.0f} rpm"),
           (f"{s['power_peak_kw']:.1f} kW", f"peak brake power, {s['rpm_power_peak']:.0f} rpm"),
           (f"{100 * s['efficiency']:.0f}%", f"brake efficiency at {s['rpm']:.0f} rpm"),
           (f"{s['peak_pressure_bar']:.0f} bar", 'peak cylinder pressure'),
           (f"{100 * s['eta_v']:.0f}%", 'volumetric efficiency'),
           (f"{s['flywheel_cm']:.1f} cm", 'flywheel radius, steel disc 50 mm')]
    card.append(fig.text(0.39, 0.78, 'What the run found', fontsize=30, fontweight='bold',
                         color=TH['text'], va='center', zorder=7))
    for j, (v, label) in enumerate(big):
        x, y = 0.39 + (j % 3) * 0.195, 0.64 - (j // 3) * 0.2
        card.append(fig.text(x, y, v, fontsize=36, fontweight='bold', color=TH['line'],
                             va='center', zorder=7))
        card.append(fig.text(x, y - 0.06, label, fontsize=15, color=TH['muted'], va='center',
                             zorder=7))
    gears = '   '.join(f'{g:.2f}' for g in s['gears'])
    card.append(fig.text(0.39, 0.2, 'Gear ratios, 1st to 8th', fontsize=15, color=TH['muted'],
                         va='center', zorder=7))
    card.append(fig.text(0.39, 0.15, gears, fontsize=24, color=TH['text'], va='center',
                         family='IBM Plex Mono', zorder=7))
    card.append(fig.text(0.39, 0.085, f"2.0 L, compression ratio 9.8  ·  {s['bsfc']:.0f} g/kWh "
                         f"at {s['rpm']:.0f} rpm  ·  cycle periodic after {s['iterations']} "
                         f"cycles", fontsize=14, color=TH['muted'], va='center', zorder=7))
    for a in card:
        a.set_alpha(0)

    def frame(t):
        # console: committed lines, then the line being typed
        shown = [(txt, line_color(txt)) for at, txt in lines if at <= t]
        for t0, t1, script in typing:
            if t0 <= t < t1:
                n = int(len(script) * (t - t0) / (t1 - t0 - 0.2)) if t1 - t0 > 0.2 else 0
                caret = '▌' if int(t * 3) % 2 == 0 else ' '
                shown.append((f'>> {script[:n]}{caret}', TH['text']))
        if not any(t0 <= t < t1 for t0, t1, _ in typing) and t < t_end:
            caret = '▌' if int(t * 2) % 2 == 0 else ' '
            if shown and shown[-1][0].startswith('>> ') and t > lines[-1][0] + 0.2:
                shown.append(('>> ' + caret, TH['text']))
        view = shown[-CONSOLE_LINES:]
        for k, row in enumerate(rows):
            txt, col = view[k] if k < len(view) else ('', TH['text'])
            row.set_text(txt)
            row.set_color(col)

        # which script is running
        running = None
        for t0, _, script in typing:
            if t0 <= t:
                running = script
        for name, (box, label) in chips.items():
            done = running and SCRIPTS.index(name) < SCRIPTS.index(running) or t >= t_card
            active = name == running and t < t_card
            box.set_edgecolor(TH['line'] if active else TH['line3'] if done else '#30363d')
            label.set_color(TH['text'] if active or done else TH['muted'])

        # figure window
        cur = None
        for k, (t0, i) in enumerate(figs):
            if t0 <= t:
                cur = (k, t0, i)
        for extra in [a for a in fig.axes if a not in own_axes]:
            extra.remove()       # a figure's twin axis from the previous frame
        ax.clear()
        ax.set_facecolor(TH['panel'])
        if cur is None:
            ax.axis('off')
            fig_title.set_text('')
            waiting.set_text('Figures open here as the scripts draw them' if t > INTRO else '')
            veil.set_alpha(0)
        else:
            k, t0, i = cur
            name, title, script, draw = FIGURES[i]
            draw(ax, d, dict(TH, bg=TH['panel']), ease((t - t0) / DRAW))
            fig_title.set_text(f'Figure {k + 1}  ·  {title}  ·  {script}')
            for j, c in enumerate(dots):
                c.set_color(TH['line'] if j == k else TH['muted'] if j < k else '#30363d')
            fade = max(0.0, 1 - (t - t0) / 0.3)
            computing = (running == 'part2' and script == 'part1.m' and t > t0 + DRAW + HOLD)
            veil.set_alpha(max(fade, 0.72 if computing else 0))
            waiting.set_text('part2.m is running the cycle and the speed sweep…' if computing else '')
        # closing card
        c = ease((t - t_card) / 0.6)
        for a in card:
            a.set_alpha(c)
        if c > 0:
            veil.set_alpha(max(veil.get_alpha(), c))
            fig_title.set_text('Results' if c > 0.5 else fig_title.get_text())
        # intro card
        intro.set_alpha(1 - ease((t - INTRO + 0.5) / 0.5))
        for a in intro_text:
            a.set_alpha(intro.get_alpha())

    intro = fig.add_artist(Rectangle((0, 0), 1, 1, transform=fig.transFigure, fc=TH['bg'],
                                     ec='none', zorder=20))
    intro_text = [
        fig.text(0.5, 0.56, 'Running the simulator', fontsize=46, fontweight='bold',
                 color=TH['text'], ha='center', va='center', zorder=21),
        fig.text(0.5, 0.47, 'constants.m  →  part1.m  →  part2.m', fontsize=24,
                 color=TH['line'], ha='center', va='center', family='IBM Plex Mono', zorder=21),
        fig.text(0.5, 0.40, '2.0 L four-cylinder spark-ignition engine  ·  GNU Octave',
                 fontsize=18, color=TH['muted'], ha='center', va='center', zorder=21),
    ]

    if only is not None:
        frame(only)
        fig.savefig(path, facecolor=TH['bg'])
        plt.close(fig)
        return
    cmd = ['ffmpeg', '-v', 'error', '-y', '-f', 'rawvideo', '-pix_fmt', 'rgba', '-s', f'{W}x{H}',
           '-r', str(FPS), '-i', '-', '-c:v', 'libx264', '-preset', 'slow', '-crf', '20',
           '-pix_fmt', 'yuv420p', '-movflags', '+faststart', str(path)]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    n = int(t_end * FPS)
    for f in range(n):
        frame(f / FPS)
        fig.canvas.draw()
        proc.stdin.write(fig.canvas.buffer_rgba().tobytes())
        if f % 150 == 0:
            print(f'  {path.name}: frame {f}/{n}', flush=True)
    proc.stdin.close()
    proc.wait()
    plt.close(fig)
    return t_card


def main():
    plt.rcParams['font.family'] = 'IBM Plex Sans'
    d = load(sys.argv[1] if len(sys.argv) > 1 else None)
    runs = parse(sys.argv[2] if len(sys.argv) > 2 else ROOT / 'docs' / 'octave-run.log')
    t_card = render(d, runs, MEDIA / 'simulation-run.mp4')
    render(d, runs, MEDIA / 'simulation-run-poster.jpg', only=t_card + 1.5)


if __name__ == '__main__':
    main()
