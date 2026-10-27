import csv
import os
import subprocess
import sys

import numpy as np
import pytest

from conftest import ELEMENTS, NODE_IDS, NODES, STEPS, syy, uy, write_model
from midas_nfxp import NfxpFile, NfxpFormatError, read_mec
from nfxp_writer import DISP, NODE, Writer, write_mec

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ---------------------------------------------------------------- index
def test_steps_and_stages(model):
    with NfxpFile(model) as r:
        assert r.complete
        assert r.run_date == '2026/10/27 12:00'
        assert len(r.steps) == len(STEPS)
        assert list(r.stages.items()) == [(1, 'Initial'), (2, 'Excavate'), (3, 'Excavate-SRM')]
        assert [s.increment for s in r.steps] == [1, 2, 1, 7]
        assert r.steps[0].factor == pytest.approx(0.5)
        assert [len(r.stage_steps('Initial')), len(r.stage_steps(3))] == [2, 1]


def test_model_type_from_mec(model):
    with NfxpFile(model) as r:
        assert r.mec_path.endswith('toy.mec')
        assert r.model_type == '2D'
    with NfxpFile(model, find_mec=False) as r:
        assert r.mec_path is None
        assert r.model_type == '2D or axisymmetric'


def test_step_selection(model):
    with NfxpFile(model) as r:
        assert r.step(-1).index == 3
        assert r.step('Initial').index == 1                # last increment of the stage
        assert r.step('initial').index == 1                # case-insensitive
        assert r.step(('Initial', 1)).index == 0           # (stage, increment)
        assert r.step((1, -1)).index == 1
        assert r.step(r.steps[2]) is r.steps[2]
        with pytest.raises(KeyError):
            r.step('No such stage')


# ---------------------------------------------------------------- results
def test_nodal(model):
    with NfxpFile(model) as r:
        for k in range(len(STEPS)):
            ids, v = r.nodal('T2 TRANSLATION', step=k)
            np.testing.assert_array_equal(ids, NODE_IDS)
            np.testing.assert_allclose(v, uy(k), rtol=1e-6)
        ids, v = r.nodal('t2 spc_force', step=0)            # names are case-insensitive
        assert v[ids <= 3].sum() == pytest.approx(30.0)


def test_element_centroid_and_corner(model):
    with NfxpFile(model) as r:
        ids, v = r.element('EFFECTIVE NORMAL-YY', step=0)
        np.testing.assert_array_equal(ids, [10, 11])
        np.testing.assert_allclose(v, syy(0, [10, 11]))
        ids, v = r.element('EFFECTIVE NORMAL-YY', step=0, corner=1)
        np.testing.assert_allclose(v, syy(0, [10, 11]) - 1)


def test_element_type_filter(model):
    with NfxpFile(model) as r:
        ids, n = r.element('AXIAL FORCE', step=-1, element_type='CBAR')
        assert list(ids) == [20] and n[0] == pytest.approx(-20000.0)
        ids, n = r.element('AXIAL FORCE', step=-1, element_type='CQUAD4')
        assert len(ids) == 0


def test_ambiguous_component_needs_category(model):
    with NfxpFile(model) as r:
        with pytest.raises(KeyError, match='category'):
            r.element('NORMAL-XX', step=0)
        assert r.element('NORMAL-XX', step=0, category='STRAIN').ids.tolist() == [11]
        assert r.element('NORMAL-XX', step=0, category='OTHER').values[0] == pytest.approx(7.0)


def test_unknown_component_suggests_similar(model):
    with NfxpFile(model) as r:
        with pytest.raises(KeyError, match='similar'):
            r.element('NORMAL-YY', step=0)


def test_components_catalogue(model):
    with NfxpFile(model) as r:
        comps = r.components(0)
        assert ('node', 'DISP', 'T2 TRANSLATION') in comps
        assert ('corner', 'STRESS', 'EFFECTIVE NORMAL-YY') in comps
        assert all(c[0] == 'node' for c in r.components(0, 'node'))


def test_history(model):
    with NfxpFile(model) as r:
        steps, v = r.history('T2 TRANSLATION', 5)
        np.testing.assert_array_equal(steps, [0, 1, 2, 3])
        np.testing.assert_allclose(v, [-0.005, -0.010, -0.015, -0.020], rtol=1e-6)
        # element 11 is excavated in stage 2, so it drops out of the history
        steps, v = r.history('EFFECTIVE NORMAL-YY', 11, location='element')
        np.testing.assert_array_equal(steps, [0, 1])


def test_active_and_deleted(model):
    with NfxpFile(model) as r:
        assert sorted(r.active_elements(0)) == [10, 11, 20]
        assert sorted(r.active_elements(-1)) == [10, 20]
        assert r.deleted_elements(0).tolist() == []
        assert r.deleted_elements(-1).tolist() == [11]


def test_factor_of_safety(model, tmp_path):
    with NfxpFile(model) as r:
        s = r.step(-1)
        assert s.srm_fos == pytest.approx(1.4321)
        assert s.fos == pytest.approx(1.4321)
        assert r.step(0).fos is None
        np.testing.assert_array_equal(s.sam['polyline_ids'], [1, 2, 3])
        assert s.sam['fos'].min() == pytest.approx(1.5)
    # without the SRM record the FoS comes from the increment label
    p = write_model(tmp_path / 'label_only.nfxp', srm_record=False)
    with NfxpFile(str(p)) as r:
        assert r.step(-1).srm_fos is None
        assert r.step(-1).fos == pytest.approx(1.4321)


# ---------------------------------------------------------------- tables / CSV
def test_table_nodes_have_coordinates(model):
    with NfxpFile(model) as r:
        t = r.table(step=0, location='node')
        assert list(t) == ['node']
        n = t['node']
        np.testing.assert_array_equal(n['id'], NODE_IDS)
        np.testing.assert_allclose(n['y'], [NODES[i][1] for i in NODE_IDS])
        assert 'DISP:T2 TRANSLATION' in n


def test_table_merges_family_chunks(model):
    with NfxpFile(model) as r:
        t = r.table(step=0, location='element')
        assert sorted(t) == ['14_CBAR', '43_CQUAD4']
        soil = t['43_CQUAD4']
        np.testing.assert_array_equal(soil['id'], [10, 11])
        np.testing.assert_allclose(soil['x'], [0.5, 1.5])          # centroids
        np.testing.assert_allclose(soil['y'], [0.5, 0.5])
        np.testing.assert_allclose(soil['STRESS:EFFECTIVE NORMAL-YY'], syy(0, [10, 11]))
        strain = soil['STRAIN:NORMAL-XX']                           # only in the second chunk
        assert np.isnan(strain[0]) and strain[1] == pytest.approx(1e-4)
        beam = t['14_CBAR']
        np.testing.assert_allclose([beam['x'][0], beam['y'][0]], [0.5, 1.0])
        assert not any(k.startswith('STRESS:ND1') for k in soil)
        corners = r.table(step=0, location='element', corners=True)['43_CQUAD4']
        assert 'STRESS:ND1- EFFECTIVE NORMAL-YY' in corners


def test_table_without_mec_uses_family_codes(model):
    with NfxpFile(model, find_mec=False) as r:
        t = r.table(step=0, location='element')
        assert sorted(t) == [14, 43]
        assert 'x' not in t[43]


def test_export_csv(model, tmp_path):
    out = tmp_path / 'out'
    with NfxpFile(model) as r:
        rows = r.export_csv(str(out), steps=[0, -1])
    names = sorted(os.listdir(str(out)))
    assert 'summary.csv' in names
    assert '000_S01_Initial_I001_nodes.csv' in names
    assert '003_S03_Excavate_SRM_I007_elements_43_CQUAD4.csv' in names
    assert '003_S03_Excavate_SRM_I007_sam_fos.csv' in names
    assert [row['fos'] for row in rows] == ['', '1.4321']
    with open(str(out / 'summary.csv')) as f:
        summary = list(csv.DictReader(f))
    assert summary[1]['at_node'] == '6' and summary[1]['sam_min_fos'] == '1.5000'


# ---------------------------------------------------------------- robustness
def test_not_a_mec_file(tmp_path):
    p = tmp_path / 'x.nfxp'
    p.write_bytes(b'\x00' * 200)
    with pytest.raises(NfxpFormatError, match='not a MecSolver'):
        NfxpFile(str(p))


def test_empty_file(tmp_path):
    p = tmp_path / 'x.nfxp'
    p.write_bytes(b'')
    with pytest.raises(NfxpFormatError, match='empty'):
        NfxpFile(str(p))


def test_unsupported_records_raise(tmp_path):
    p = Writer().stage(1, 'Eigen', last=True).raw(b'\x63\x00\x00\x00' * 8).eof().save(tmp_path / 'e.nfxp')
    with pytest.raises(NfxpFormatError, match='not supported'):
        NfxpFile(str(p))


def test_interrupted_run_warns(tmp_path):
    p = write_model(tmp_path / 'cut.nfxp', eof=False)
    with pytest.warns(UserWarning, match='end-of-file'):
        r = NfxpFile(str(p))
    assert not r.complete and len(r.steps) == len(STEPS)
    r.close()


def test_unknown_record_after_results_warns(tmp_path):
    p = tmp_path / 'tail.nfxp'
    write_model(p, eof=False)
    with open(str(p), 'ab') as f:
        f.write(b'\x63\x00\x00\x00' * 4)
    with pytest.warns(UserWarning, match='unknown record type 99'):
        r = NfxpFile(str(p))
    assert len(r.steps) == len(STEPS)
    r.close()


def test_linear_static_single_step(tmp_path):
    w = Writer().stage(1, 'Linear', last=True).ids(1, [1, 2]).result(2, NODE, DISP, 'T2 TRANSLATION', [0.0, -0.1])
    p = w.eof().save(tmp_path / 'lin.nfxp')
    with NfxpFile(str(p)) as r:
        assert len(r.steps) == 1
        s = r.steps[0]
        assert (s.stage_label, s.increment, s.factor, s.label) == ('Linear', 1, 1.0, '')
        assert r.nodal('T2 TRANSLATION').values[1] == pytest.approx(-0.1)


def test_utf8_labels(tmp_path):
    w = Writer().stage(1, 'Котлован', last=True).increment(1, 1.0, 'INCR=1 (LOAD=1.000)', 0)
    p = w.ids(1, [1]).result(2, NODE, DISP, 'T2 TRANSLATION', [0.0]).eof().save(tmp_path / 'ru.nfxp')
    with NfxpFile(str(p)) as r:
        assert r.steps[0].stage_label == 'Котлован'
        assert r.step('Котлован').index == 0


# ---------------------------------------------------------------- .mec reader
def test_read_mec_types(tmp_path):
    m = read_mec(str(write_mec(tmp_path / 'a.mec', NODES, ELEMENTS)))
    assert m['model_type'] == '2D'
    assert m['element_type'][20] == 'CBAR' and m['connectivity'][10] == [1, 2, 5, 4]
    axi = {1: ('CQUADX', [1, 2, 5, 4])}
    assert read_mec(str(write_mec(tmp_path / 'b.mec', NODES, axi)))['model_type'] == 'axisymmetric'
    nodes3d = dict(NODES)
    nodes3d.update({7: (0., 0., 1.), 8: (1., 0., 1.), 9: (1., 1., 1.), 10: (0., 1., 1.)})
    hexa = {1: ('CHEXA', [1, 2, 5, 4, 7, 8, 9, 10])}             # needs a continuation line
    m = read_mec(str(write_mec(tmp_path / 'c.mec', nodes3d, hexa)))
    assert m['model_type'] == '3D' and m['connectivity'][1] == [1, 2, 5, 4, 7, 8, 9, 10]


# ---------------------------------------------------------------- CLI
def _cli(*args):
    return subprocess.run([sys.executable, os.path.join(ROOT, 'midas_nfxp.py')] + list(args),
                          capture_output=True, text=True, encoding='utf-8')


def test_cli_info(model):
    p = _cli('info', model)
    assert p.returncode == 0, p.stderr
    assert '2D, 3 stages, 4 steps' in p.stdout
    assert 'FoS 1.4321' in p.stdout


def test_cli_history(model):
    p = _cli('history', model, 'T2 TRANSLATION', '5')
    assert p.returncode == 0, p.stderr
    lines = p.stdout.strip().splitlines()
    assert lines[0] == 'step,stage,increment,label,value' and len(lines) == 5


def test_cli_export(model, tmp_path):
    out = str(tmp_path / 'cli_out')
    p = _cli('export', model, '--out', out)
    assert p.returncode == 0, p.stderr
    assert len([f for f in os.listdir(out) if f.endswith('_nodes.csv')]) == 3     # last step of each stage
