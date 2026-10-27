import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

from nfxp_writer import (CORNER, DISP, ELEMENT, FORCE, NODE, OTHER, STRAIN, STRESS,  # noqa: E402
                         Writer, write_mec)

# 2D plane-strain toy model: two quads and one beam along the top edge
#
#   4 ---- 5 ---- 6        beam 20 = 4-5
#   | 10   | 11   |
#   1 ---- 2 ---- 3
NODES = {1: (0., 0., 0.), 2: (1., 0., 0.), 3: (2., 0., 0.),
         4: (0., 1., 0.), 5: (1., 1., 0.), 6: (2., 1., 0.)}
ELEMENTS = {10: ('CQUAD4', [1, 2, 5, 4]), 11: ('CQUAD4', [2, 3, 6, 5]), 20: ('CBAR', [4, 5])}
NODE_IDS = np.array(sorted(NODES))

# (stage, stage label, increment label, load factor, deleted elements)
STEPS = [
    (1, 'Initial', 'INCR=1 (LOAD=0.500)', 0.5, []),
    (1, 'Initial', 'INCR=2 (LOAD=1.000)', 1.0, []),
    (2, 'Excavate', 'INCR=1 (LOAD=1.000)', 1.0, [11]),
    (3, 'Excavate-SRM', 'INCR=7 (FOS=1.4321)', 1.4321, [11]),
]


def uy(k):
    """Vertical displacement written for global step k (0-based) at every node."""
    return -0.001 * (k + 1) * NODE_IDS


def syy(k, ids):
    return -1000.0 * (k + 1) - np.asarray(ids, float)


def write_model(path, srm_record=True, eof=True):
    w = Writer()
    stage = None
    for k, (st, slabel, ilabel, factor, deleted) in enumerate(STEPS):
        if st != stage:
            w.stage(st, slabel, last=(st == STEPS[-1][0]))
            stage = st
        incr = int(ilabel.split('=')[1].split()[0])
        w.increment(k + 1, factor, ilabel, incr - 1)
        w.ids(1, NODE_IDS)
        w.result(1, NODE, DISP, 'TOTAL TRANSLATION', np.abs(uy(k)))
        w.result(2, NODE, DISP, 'T2 TRANSLATION', uy(k))
        w.result(4, NODE, DISP, 'T3 TRANSLATION', np.zeros(len(NODE_IDS)))
        w.result(3, NODE, FORCE, 'T2 SPC_FORCE', np.where(NODE_IDS <= 3, 10.0 * (k + 1), 0.0))
        quads = [e for e in (10, 11) if e not in deleted]
        w.ids(3, quads + [20])
        if deleted:
            w.deleted(deleted)
        # soil family 43 (plane strain) in two chunks with different columns
        w.ids(2, quads)
        w.result(430001, ELEMENT, STRESS, 'EFFECTIVE NORMAL-YY', syy(k, quads))
        w.result(430002, CORNER, STRESS, 'ND1- EFFECTIVE NORMAL-YY', syy(k, quads) - 1)
        if 11 in quads:
            w.ids(2, [11])
            w.result(430003, ELEMENT, STRAIN, 'NORMAL-XX', [1e-4 * (k + 1)])
        # beam family 14
        w.ecsys(1)
        w.ids(2, [20])
        w.result(140001, ELEMENT, FORCE, 'AXIAL FORCE', [-5000.0 * (k + 1)])
        w.result(140002, ELEMENT, OTHER, 'NORMAL-XX', [7.0])
        if st == 3:
            if srm_record:
                w.srm(1.4321)
            w.sam([1, 2, 3], [1.9, 1.5, 1.7])
    if eof:
        w.eof()
    return w.save(path)


@pytest.fixture
def model(tmp_path):
    """Path to a synthetic .nfxp with a matching .mec next to it."""
    p = write_model(tmp_path / 'toy.nfxp')
    write_mec(tmp_path / 'toy.mec', NODES, ELEMENTS)
    return str(p)
