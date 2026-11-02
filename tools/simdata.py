"""Loads the simulator's saved workspace (part2.mat) and the shared look of every render.

The renderers in this folder draw only what the MATLAB scripts computed: they read
part2.mat and never recompute the physics. Run the scripts first (MATLAB or GNU Octave),
then the renderers.
"""

from pathlib import Path

import numpy as np
import scipy.io as sio

ROOT = Path(__file__).resolve().parent.parent
MEDIA = ROOT / 'docs' / 'media'

# Engine crank angle at which each cylinder starts its own cycle. Cylinder 1 follows the
# simulated cycle directly (it fires at 360 degrees); the others are the same cycle shifted
# so that they fire every 180 degrees in the order 1-3-4-2, as in part2.m and animation.m.
OFFSET = {1: 0, 2: 540, 3: 180, 4: 360}
FIRING_ORDER = (1, 3, 4, 2)

LIGHT = dict(
    bg='#ffffff', text='#2b2b2b', muted='#6b6b6b', grid='#e8e8e8', axis='#9a9a9a',
    line='#e55e26', line2='#2f6fd0', line3='#3c9a5f', line4='#8a5cc2',
    phase={'intake': '#eef3fb', 'compression': '#f5f5f0', 'combustion': '#fbe3d8',
           'expansion': '#f5f5f0', 'exhaust': '#eef3fb'},
    phase_alpha=1.0,
    faint=0.9,
)

DARK = dict(
    bg='#0d1117', panel='#161b22', text='#e6edf3', muted='#8b949e', grid='#262c36', axis='#4a5260',
    line='#ff7b3a', line2='#58a6ff', line3='#3fb950', line4='#d2a8ff',
    phase={'intake': '#58a6ff', 'compression': '#3fb950', 'combustion': '#ff5a36',
           'expansion': '#f2cc60', 'exhaust': '#8b949e'},
    phase_alpha=0.10,
    faint=0.6,
)

CYLINDER_COLOR = {1: '#ff7b3a', 2: '#58a6ff', 3: '#3fb950', 4: '#d2a8ff'}


def load(path=None):
    """The workspace part2.m saves, as a dict of numpy values."""
    path = Path(path) if path else ROOT / 'part2.mat'
    raw = sio.loadmat(path, squeeze_me=True)
    return {k: v for k, v in raw.items() if not k.startswith('__')}


def phases(d):
    """The five phases of the simulated cycle as (name, first degree, last degree)."""
    return [
        ('intake', int(d['intake_start']), int(d['intake_end'])),
        ('compression', int(d['compression_start']), int(d['compression_end'])),
        ('combustion', int(d['combustion_start']), int(d['combustion_end'])),
        ('expansion', int(d['expansion_start']), int(d['expansion_end'])),
        ('exhaust', int(d['exhaust_start']), int(d['exhaust_end'])),
    ]


def phase_at(d, phi):
    """Name of the phase a cylinder is in at its own crank angle phi (degrees)."""
    x = phi % 720
    if x < 1:
        x += 720
    for name, a, b in phases(d):
        if a - 0.5 <= x < b + 0.5:
            return name
    return 'exhaust'


def at(arr, phi):
    """Value of a per-degree array (index 0 is 1 degree, index 719 is 720) at any angle."""
    arr = np.asarray(arr, dtype=float)
    n = len(arr)
    f = (np.asarray(phi, dtype=float) - 1.0) % n
    i0 = np.floor(f).astype(int)
    i1 = (i0 + 1) % n
    t = f - i0
    return arr[i0] * (1 - t) + arr[i1] * t


def cylinder(arr, k, theta):
    """A cylinder-1 quantity seen in cylinder k at engine crank angle theta."""
    return at(arr, np.asarray(theta) - OFFSET[k])


def burned_fraction(d, phi):
    """Wiebe burned fraction at a cylinder's own crank angle, as part2.m defines it."""
    a0, a1 = float(d['combustion_start']), float(d['combustion_end'])
    x = (phi % 720 - a0) / (a1 - a0)
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    return 1 - np.exp(-float(d['a']) * x ** (float(d['m']) + 1))


def summary(d):
    """The headline figures the scripts print, rounded for display."""
    p_atm = float(d['p_atm'])
    return {
        'power_kw': float(d['Engine_power']),
        'mean_torque': float(d['T_avg_total']),
        'efficiency': float(d['fuel_conv_efficiency']),
        'work': float(d['engine_work']),
        'peak_pressure_atm': float(np.max(d['pressure']) / p_atm),
        'flywheel_cm': float(d['flywheel_radius']) * 100,
        'flywheel_energy': float(d['energy']),
        'gears': [float(g) for g in np.atleast_1d(d['iA'])],
        'rpm': float(d['rpm']),
        'iterations': int(d['counter']),
    }
