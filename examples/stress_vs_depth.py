"""Vertical stress of the soil elements against depth for one step: a quick check of
the in-situ stage (the slope should be close to the unit weight). Needs the .mec deck.

    python examples/stress_vs_depth.py model.nfxp [step] [--plot]

Uses EFFECTIVE NORMAL-ZZ in 3D, EFFECTIVE NORMAL-YY in 2D and NORMAL-YY in
axisymmetric models (whichever exists), at element centroids.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from midas_nfxp import NfxpFile  # noqa: E402

args = [x for x in sys.argv[1:] if x != '--plot']
step = args[1] if len(args) > 1 else '0'
step = int(step) if step.lstrip('-').isdigit() else step

with NfxpFile(args[0]) as r:
    vertical = 'z' if r.model_type == '3D' else 'y'
    names = (['EFFECTIVE NORMAL-ZZ', 'TOTAL NORMAL-ZZ'] if vertical == 'z'
             else ['EFFECTIVE NORMAL-YY', 'TOTAL NORMAL-YY', 'NORMAL-YY'])
    st = r.step(step)
    for name in names:
        try:
            ids, sv = r.element(name, st, category='STRESS')
            break
        except KeyError:
            continue
    else:
        sys.exit('no vertical stress component in %r' % st)
    xyz = r._centroids(ids)

level = xyz[:, 2 if vertical == 'z' else 1]
ok = np.isfinite(level) & np.isfinite(sv)
slope, icpt = np.polyfit(level[ok], sv[ok], 1)
r2 = np.corrcoef(level[ok], sv[ok])[0, 1] ** 2
print('%s, %s: %d elements' % (st, name, ok.sum()))
print('  d(sigma)/d%s = %.1f per unit length (R^2 = %.3f)' % (vertical, slope, r2))
print('  For level ground this is the unit weight (N/m3 in MEC units); with effective stress below the')
print('  water table it is the buoyant unit weight. Sloping ground or layering lowers R^2.')

if '--plot' in sys.argv:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5, 6))
    ax.scatter(sv[ok] / 1e3, level[ok], s=4)
    ax.set_xlabel('%s / 1000' % name)
    ax.set_ylabel(vertical)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    png = os.path.splitext(args[0])[0] + '_stress_vs_depth.png'
    fig.savefig(png, dpi=150)
    print('written', png)
