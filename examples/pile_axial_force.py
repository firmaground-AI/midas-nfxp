"""Axial force along each pile (vertical beam elements) for one step, written to CSV
and optionally plotted. Needs the .mec deck next to the .nfxp for coordinates.

    python examples/pile_axial_force.py model.nfxp                 # last step
    python examples/pile_axial_force.py model.nfxp --step "Surface load" --plot

Piles are found by grouping beam-element centroids with the same (x, y). Units are
those of the model (N, m for MEC decks); --scale 1e-3 prints kN.
"""
import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from midas_nfxp import NfxpFile  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('nfxp')
ap.add_argument('--step', default='-1', help='step index or stage label (default: last step)')
ap.add_argument('--family', default='14', help='beam family code (14 = CBAR beams)')
ap.add_argument('--scale', type=float, default=1e-3, help='multiply forces by this (default 1e-3: N -> kN)')
ap.add_argument('--out', help='CSV path (default: <nfxp>_pile_axial_force.csv)')
ap.add_argument('--plot', action='store_true', help='plot N against z (needs matplotlib)')
a = ap.parse_args()
step = int(a.step) if a.step.lstrip('-').isdigit() else a.step

with NfxpFile(a.nfxp) as r:
    st = r.step(step)
    tables = r.table(st, location='element')
    fam = [k for k in tables if str(k).split('_')[0] == a.family]
    if not fam:
        sys.exit('no element family %s in %r; families are %s' % (a.family, st, list(tables)))
    piles = tables[fam[0]]

x, y, z = piles['x'], piles['y'], piles['z']
N = piles['FORCE:AXIAL FORCE'] * a.scale
heads = sorted(set(zip(np.round(x, 3), np.round(y, 3))))
out = a.out or os.path.splitext(a.nfxp)[0] + '_pile_axial_force.csv'

print('%s: %d beam elements in %d piles' % (st, len(N), len(heads)))
with open(out, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['pile', 'x', 'y', 'element', 'z', 'axial_force'])
    profiles = []
    for k, (px, py) in enumerate(heads, 1):
        on = np.nonzero((np.round(x, 3) == px) & (np.round(y, 3) == py))[0]
        on = on[np.argsort(-z[on])]                     # head to toe
        profiles.append((k, px, py, z[on], N[on]))
        print('  pile %d at (%g, %g): N = %.1f at head, %.1f at toe' % (k, px, py, N[on][0], N[on][-1]))
        for i in on:
            w.writerow([k, px, py, int(piles['id'][i]), '%.4f' % z[i], '%.4f' % N[i]])
print('written', out)

if a.plot:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 7))
    for k, px, py, zz, nn in profiles:
        ax.plot(nn, zz, marker='o', ms=3, label='pile %d (%g, %g)' % (k, px, py))
    ax.set_xlabel('axial force (compression negative)')
    ax.set_ylabel('z')
    ax.set_title(str(st))
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    png = os.path.splitext(out)[0] + '.png'
    fig.savefig(png, dpi=150)
    print('written', png)
