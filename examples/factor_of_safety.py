"""Strength-reduction (SRM) and slope-analysis (SAM) factors of safety per stage.

    python examples/factor_of_safety.py model.nfxp
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from midas_nfxp import NfxpFile  # noqa: E402

with NfxpFile(sys.argv[1]) as r:
    print('%-24s %8s %8s' % ('stage', 'SRM FoS', 'SAM FoS'))
    for st, label in r.stages.items():
        steps = r.stage_steps(st)
        srm = [s.fos for s in steps if s.fos is not None]
        sam = [s.sam['fos'].min() for s in steps if s.sam is not None]
        if srm or sam:
            print('%-24s %8s %8s' % (label, '%.3f' % srm[-1] if srm else '-', '%.3f' % sam[-1] if sam else '-'))
