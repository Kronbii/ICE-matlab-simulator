"""Writes every figure the scripts produce to docs/media as WebP.

    python3 tools/render_figures.py [path/to/part2.mat]
"""

import io
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

from plots import FIGURES
from simdata import LIGHT, MEDIA, load


def main():
    d = load(sys.argv[1] if len(sys.argv) > 1 else None)
    MEDIA.mkdir(parents=True, exist_ok=True)
    plt.rcParams['font.family'] = 'DejaVu Sans'
    for name, _, _, draw in FIGURES:
        height = 8.8 if name == 'pv-diagram' else 8.0
        fig = plt.figure(figsize=(12.8, height), dpi=100, facecolor=LIGHT['bg'])
        ax = fig.add_axes([0.085, 0.1, 0.875, 0.8])
        draw(ax, d, LIGHT)
        buf = io.BytesIO()
        fig.savefig(buf, format='png', facecolor=LIGHT['bg'])
        plt.close(fig)
        out = MEDIA / f'{name}.webp'
        Image.open(buf).convert('RGB').save(out, 'WEBP', quality=90, method=6)
        print(out.relative_to(MEDIA.parent.parent))


if __name__ == '__main__':
    main()
