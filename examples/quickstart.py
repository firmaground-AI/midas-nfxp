"""Open a result file and print what is in it.

    python examples/quickstart.py path/to/model.nfxp
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from midas_nfxp import NfxpFile  # noqa: E402

path = sys.argv[1]
with NfxpFile(path) as r:
    print(r)
    print('mec deck:', r.mec_path or 'not found (no coordinates)')

    print('\nSaved steps')
    for s in r.steps:
        fos = '  FoS %.3f' % s.fos if s.fos is not None else ''
        print('  %3d  stage %2d %-20s %s%s' % (s.index, s.stage, s.stage_label, s.label, fos))

    print('\nResults in the last step')
    for loc, cat, comp in r.components(-1):
        if loc != 'corner':
            print('  %-8s %-7s %s' % (loc, cat, comp))

    print('\nMaximum displacement at the end of each stage')
    for st, label in r.stages.items():
        s = r.stage_steps(st)[-1]
        ids, u = r.nodal('TOTAL TRANSLATION', s)
        k = int(np.argmax(u))
        print('  %-20s %9.2f mm at node %d' % (label, u[k] * 1e3, ids[k]))
