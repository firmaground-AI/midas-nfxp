"""Optional checks on real GTS NX result files, which are not part of the repository.

Point NFXP_TEST_FILES at one or more .nfxp files (separated by os.pathsep, i.e. ';' on
Windows and ':' elsewhere) to run them:

    set NFXP_TEST_FILES=C:\\models\\a.nfxp;C:\\models\\b.nfxp
    python -m pytest tests/test_real_files.py -v
"""
import os

import numpy as np
import pytest

from midas_nfxp import NfxpFile

FILES = [p for p in os.environ.get('NFXP_TEST_FILES', '').split(os.pathsep) if p.strip()]
pytestmark = pytest.mark.skipif(not FILES, reason='set NFXP_TEST_FILES to run on real models')


@pytest.mark.parametrize('path', FILES, ids=[os.path.basename(p) for p in FILES])
def test_real_file(path):
    with NfxpFile(path) as r:
        assert r.steps, 'no steps indexed'
        assert r.model_type in ('2D', 'axisymmetric', '3D', '2D or axisymmetric')
        last = r.step(-1)
        assert last.blocks, 'last step has no results'
        f = r.nodal('TOTAL TRANSLATION', last)
        assert len(f.ids) == len(np.unique(f.ids))
        assert np.all(np.isfinite(f.values)) and np.all(f.values >= 0)
        # every result block points inside the file
        size = os.path.getsize(path)
        for s in r.steps:
            for b in s.blocks:
                assert b.offset + 4 * b.n <= size
        if r.mec_path:
            node = r.table(last, 'node')['node']
            assert np.isfinite(node['x']).all(), 'result nodes missing from the .mec deck'
