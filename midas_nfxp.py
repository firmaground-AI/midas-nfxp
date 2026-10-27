"""Python API for MIDAS MEC solver result files (*.nfxp) written by GTS NX.

Works for 2D (plane strain), axisymmetric and 3D models, nonlinear / construction
stage / SRM / seepage runs with every saved increment, and linear static runs.
Files are memory-mapped and only indexed on open, so multi-GB 3D results open
in seconds; values are read when asked for.

    from midas_nfxp import NfxpFile
    with NfxpFile('model_CS.nfxp') as r:            # finds model_CS.mec next to it
        print(r.model_type)                         # '2D' | 'axisymmetric' | '3D'
        for s in r.steps: print(s)
        ids, uy = r.nodal('T2 TRANSLATION', step=-1)
        ids, sxx = r.element('TOTAL NORMAL-XX', step='Stage-2')
        steps, u = r.history('TOTAL TRANSLATION', 1702)
        r.export_csv('out')                         # one CSV per step / element family

Command line:  midas-nfxp info|export|history ...  (or python midas_nfxp.py; see --help)

The format is undocumented. The layout below was reverse-engineered (MEC V1.9.8)
and checked against the solvers' .out/.log summaries; an unknown record stops
indexing with a warning instead of guessing. Little-endian, int32 / float32:

  header          0x00..0x51 ("MecSolver", run date)
  10017 STAGE     108 B: type, 0, next-stage offset (-1 = last stage), ?, 12,
                  stage no., analysis code, stage label[80]
  13 INCREMENT    100 B: type, global step no., load factor / FoS (f32), analysis
                  code, increment label[80], increment no. in stage.
                  Every saved increment has one; linear static has none.
  1 / 2 / 3 IDS   kind, n, n ids. 1 = node ids, 2 = element ids of the element
                  family whose results follow, 3 = active elements of the stage
  4 DELETED       n, n element ids removed so far
  10003 ECSYS     n, n x 6 f32 (element axes)
  14 RESULT       18-int header, label[80], n f32 values. Header[1] = result id,
                  [2] = location (0 node, 1 element centroid, 2 element node
                  'NDk- ' in label), [3] = category (0 other, 1 displacement,
                  3 force, 4 stress, 5 strain). Values follow the latest IDS of
                  kind 1 (nodes) or 2 (elements).
  10025 SAM FOS   ?, post id, ?, ?, n, n polyline ids, n f64 FoS
  10032 SRM FOS   pad, f64 critical FoS
  EOF             10017, 0, -1, -1, -1, 0 (24 B)

Not covered: Civil NX / eigenvalue / spectrum / moving-load result files, which
use other record types (opening them raises NfxpFormatError).
Axisymmetric reactions are for the full 360 deg (no 2*pi factor needed).
"""
import argparse
import csv
import mmap
import os
import re
import struct
import warnings
from collections import OrderedDict, namedtuple

import numpy as np

__version__ = '0.1.0'
__all__ = ['NfxpFile', 'NfxpFormatError', 'Step', 'Block', 'Field', 'read_mec']

CATEGORY ={0: 'OTHER', 1: 'DISP', 3: 'FORCE', 4: 'STRESS', 5: 'STRAIN'}
LOCATION = {0: 'node', 1: 'element', 2: 'corner'}
SOLID = {'CTETRA', 'CPENTA', 'CHEXA', 'CPYRAM'}
AXISYMMETRIC = {'CQUADX', 'CTRIAX'}

Field = namedtuple('Field', 'ids values')


class NfxpFormatError(Exception):
    pass


def _cstr(b):
    b = bytes(b).split(b'\x00')[0]
    try:
        return b.decode('utf-8').strip()        # GTS NX writes labels as UTF-8
    except UnicodeDecodeError:
        return b.decode('latin1').strip()


class Block(object):
    """One result array: a component at nodes, element centroids or one element corner."""
    __slots__ = ('result_id', 'location', 'category', 'label', 'component', 'corner',
                 'ids_offset', 'n', 'offset')

    def __init__(self, h, label, ids_offset, n, offset):
        self.result_id, self.n, self.offset, self.ids_offset = h[1], n, offset, ids_offset
        self.location = LOCATION.get(h[2], h[2])
        self.category = CATEGORY.get(h[3], h[3])
        self.label = label
        m = re.match(r'ND(\d+)-\s*(.*)', label)
        self.corner, self.component = (int(m.group(1)), m.group(2)) if m and h[2] == 2 else (None, label)

    @property
    def family(self):
        """Element-family code from the result id (e.g. 31 solid, 41 shell, 14 beam)."""
        return self.result_id % 100000000 // 10000


class Step(object):
    """One saved solution: a stage increment (or the single linear-static solution)."""

    def __init__(self, index, stage, stage_label, analysis, step_no, increment, factor, label):
        self.index, self.stage, self.stage_label, self.analysis = index, stage, stage_label, analysis
        self.step_no, self.increment, self.factor, self.label = step_no, increment, factor, label
        self.blocks = []
        self.active = self.deleted = None       # (offset, n) of id lists
        self.srm_fos = None
        self.sam = None

    @property
    def fos(self):
        """Strength-reduction factor of safety: the SRM record when written, else the
        FoS that SRM increments carry as their load factor ('INCR=n (FOS=x)')."""
        if self.srm_fos is not None:
            return self.srm_fos
        m = re.search(r'FOS=\s*([-+0-9.Ee]+)', self.label)
        return float(m.group(1)) if m else None

    def __repr__(self):
        return '<Step %d: stage %d %r, %s>' % (self.index, self.stage, self.stage_label,
                                               self.label or 'factor %g' % self.factor)


class NfxpFile(object):
    def __init__(self, path, mec=None, find_mec=True):
        self.path = path
        if os.path.getsize(path) < 0x52 + 24:
            raise NfxpFormatError('empty or incomplete result file: %s' % path)
        self._fh = open(path, 'rb')
        self._mm = mmap.mmap(self._fh.fileno(), 0, access=mmap.ACCESS_READ)
        if self._mm[0x10:0x19] != b'MecSolver':
            self.close()
            raise NfxpFormatError('not a MecSolver result file: %s' % path)
        self.run_date = ' '.join(bytes(self._mm[0x2c:0x3c]).replace(b'\x00', b' ').decode('latin1').split())
        if mec is None and find_mec:
            cand = os.path.splitext(path)[0] + '.mec'
            mec = cand if os.path.exists(cand) else None
        self.mec_path = mec
        self._mesh = None
        self.steps = []
        self.complete = False       # True when the end-of-file record was reached
        self._index()

    # ------------------------------------------------------------------ indexing
    def _index(self):
        d, size = self._mm, len(self._mm)
        u = struct.unpack_from
        o = 0x52
        ids = {}                    # kind -> (offset of first id, n)
        stage = None
        step = None
        problem = None
        while o + 4 <= size:
            t = u('<i', d, o)[0]
            if t == 10017:
                if o + 108 > size or u('<i', d, o + 16)[0] == -1:
                    self.complete = o + 20 <= size and u('<3i', d, o + 8) == (-1, -1, -1)
                    break
                stage = (u('<i', d, o + 20)[0], _cstr(d[o + 28:o + 108]), u('<i', d, o + 24)[0])
                o += 108
                if o + 4 > size or u('<i', d, o)[0] != 13:     # no increments: linear static
                    step = Step(len(self.steps), stage[0], stage[1], stage[2], len(self.steps), 1, 1.0, '')
                    self.steps.append(step)
            elif t == 13:
                if stage is None or o + 100 > size:
                    problem = 'increment record outside a stage at 0x%x' % o
                    break
                no, = u('<i', d, o + 4)
                f, = u('<f', d, o + 8)
                label = _cstr(d[o + 16:o + 96])
                m = re.match(r'INCR=\s*(\d+)', label)       # 1-based, as in GTS NX
                inc = int(m.group(1)) if m else u('<i', d, o + 96)[0] + 1
                step = Step(len(self.steps), stage[0], stage[1], stage[2], no, inc, f, label)
                self.steps.append(step)
                o += 100
            elif t in (1, 2, 3) and o + 8 <= size and 0 <= u('<i', d, o + 4)[0] < 200000000:
                n = u('<i', d, o + 4)[0]
                if o + 8 + 4 * n > size:
                    break
                ids[t] = (o + 8, n)
                if t == 3 and step is not None:
                    step.active = ids[3]
                o += 8 + 4 * n
            elif t == 4:
                n = u('<i', d, o + 4)[0]
                if step is not None:
                    step.deleted = (o + 8, n)
                o += 8 + 4 * n
            elif t == 10003:
                o += 8 + 24 * u('<i', d, o + 4)[0]
            elif t == 14:
                if o + 152 > size:
                    break
                h = u('<18i', d, o)
                key = 1 if h[2] == 0 else 2
                if key not in ids or step is None:
                    problem = 'result block before its id list at 0x%x' % o
                    break
                ioff, n = ids[key]
                if o + 152 + 4 * n > size:
                    break
                step.blocks.append(Block(h, _cstr(d[o + 72:o + 152]), ioff, n, o + 152))
                o += 152 + 4 * n
            elif t == 10025:
                n = u('<i', d, o + 20)[0]
                step.sam = dict(post_id=u('<i', d, o + 8)[0],
                                polyline_ids=np.frombuffer(d, '<i4', n, o + 24).copy(),
                                fos=np.frombuffer(d, '<f8', n, o + 24 + 4 * n).copy())
                o += 24 + 12 * n
            elif t == 10032:
                step.srm_fos = u('<d', d, o + 8)[0]
                o += 16
            elif t == 0:
                o += 4
            else:
                problem = 'unknown record type %d at 0x%x' % (t, o)
                break
        if not self.steps or (problem and not any(s.blocks for s in self.steps)):
            self.close()
            raise NfxpFormatError('%s: %s' % (self.path, problem or 'no results found')
                                  + ' (Civil NX, eigenvalue, spectrum and moving-load files are not supported)')
        if problem:
            warnings.warn('%s: %s; indexed %d steps before it' % (os.path.basename(self.path), problem, len(self.steps)))
        elif not self.complete:
            warnings.warn('%s: file ends without the end-of-file record (run interrupted or still '
                          'writing); the last step may be incomplete' % os.path.basename(self.path))

    def close(self):
        if getattr(self, '_mm', None) is not None:
            self._mm.close()
            self._fh.close()
            self._mm = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def __repr__(self):
        return '<NfxpFile %s: %s, %d stages, %d steps>' % (
            os.path.basename(self.path), self.model_type, len(self.stages), len(self.steps))

    # ------------------------------------------------------------------ raw access
    def _ids(self, ref):
        off, n = ref
        return np.frombuffer(self._mm, '<i4', n, off).copy()

    def _values(self, b):
        return np.frombuffer(self._mm, '<f4', b.n, b.offset).astype(np.float64)

    # ------------------------------------------------------------------ steps
    @property
    def stages(self):
        """OrderedDict stage number -> stage label."""
        return OrderedDict((s.stage, s.stage_label) for s in self.steps)

    def step(self, key=-1):
        """A Step from an index (negative counts from the end), a stage label or
        number (its last saved increment), or a (stage, increment) tuple."""
        if isinstance(key, Step):
            return key
        if isinstance(key, (int, np.integer)):
            return self.steps[key]
        if isinstance(key, tuple):
            st, inc = key
            hits = [s for s in self.steps if st in (s.stage, s.stage_label)]
            if not hits:
                raise KeyError('no stage %r' % (st,))
            return hits[inc] if inc < 0 else next(s for s in hits if s.increment == inc)
        hits = [s for s in self.steps if s.stage_label == key or s.stage_label.lower() == str(key).lower()]
        if not hits:
            raise KeyError('no stage %r; stages are %s' % (key, list(self.stages.values())))
        return hits[-1]

    def stage_steps(self, stage):
        return [s for s in self.steps if stage in (s.stage, s.stage_label)]

    # ------------------------------------------------------------------ catalogue
    def components(self, step=-1, location=None):
        """[(location, category, component)] available in a step."""
        seen = OrderedDict()
        for b in self.step(step).blocks:
            if location in (None, b.location):
                seen[(b.location, b.category, b.component)] = None
        return list(seen)

    def _select(self, step, component, location, category, corner):
        st = self.step(step)
        comp = component.upper()
        hits = [b for b in st.blocks if b.location == location and b.component.upper() == comp
                and (corner is None or b.corner == corner)]
        cats = sorted({b.category for b in hits}, key=str)
        if category is not None:
            hits = [b for b in hits if b.category == category.upper()]
        elif len(cats) > 1:
            raise KeyError('%r exists as %s; pass category=' % (component, cats))
        if not hits:
            near = sorted({b.component for b in st.blocks if b.location == location and comp in b.component.upper()})
            raise KeyError('no %s result %r in %r%s' % (location, component, st,
                                                        '; similar: %s' % near if near else ''))
        return hits

    # ------------------------------------------------------------------ results
    def nodal(self, component, step=-1, category=None):
        """Field(ids, values) of a nodal result, e.g. 'TOTAL TRANSLATION', 'T2 TRANSLATION',
        'T2 SPC_FORCE'."""
        b = self._select(step, component, 'node', category, None)[0]
        return Field(self._ids((b.ids_offset, b.n)), self._values(b))

    def element(self, component, step=-1, category=None, corner=None, element_type=None):
        """Field(ids, values) of an element result at the centroid (corner=None) or at
        element node k (corner=k). Families that share the component name (e.g. pile and
        pile tip) are concatenated; element_type ('CHEXA', ...) filters with the .mec."""
        loc = 'element' if corner is None else 'corner'
        blocks = self._select(step, component, loc, category, corner)
        ids = np.concatenate([self._ids((b.ids_offset, b.n)) for b in blocks])
        vals = np.concatenate([self._values(b) for b in blocks])
        if element_type:
            types = self.mesh['element_type']
            want = {element_type.upper()} if isinstance(element_type, str) else {t.upper() for t in element_type}
            keep = np.array([types.get(int(e)) in want for e in ids], bool)
            ids, vals = ids[keep], vals[keep]
        return Field(ids, vals)

    def history(self, component, id, location='node', category=None, corner=None, steps=None):
        """(step indices, values) of one node / element over the saved steps.
        Steps where the entity is inactive or the result absent are skipped."""
        out_s, out_v = [], []
        for s in (self.steps if steps is None else [self.step(k) for k in steps]):
            try:
                f = (self.nodal(component, s, category) if location == 'node'
                     else self.element(component, s, category, corner))
            except KeyError:
                continue
            k = np.nonzero(f.ids == id)[0]
            if len(k):
                out_s.append(s.index)
                out_v.append(f.values[k[0]])
        return np.array(out_s, int), np.array(out_v)

    def active_elements(self, step=-1):
        st = self.step(step)
        return self._ids(st.active) if st.active else None

    def deleted_elements(self, step=-1):
        st = self.step(step)
        return self._ids(st.deleted) if st.deleted else np.empty(0, np.int32)

    def table(self, step=-1, location='node', corners=False):
        """Columns of every result at one location, grouped by element family for elements:
        {family: OrderedDict(column -> array)}; nodes give a single group keyed 'node'.
        Adds x, y, z (node coordinates or element centroids) when the .mec is available."""
        st = self.step(step)
        groups = OrderedDict()
        for b in st.blocks:
            if location == 'node' and b.location != 'node':
                continue
            if location == 'element' and not (b.location == 'element' or (corners and b.location == 'corner')):
                continue
            key = 'node' if b.location == 'node' else (b.ids_offset, b.n)
            g = groups.setdefault(key, OrderedDict(id=self._ids((b.ids_offset, b.n)), _family=b.family))
            name = '%s:%s' % (b.category, b.label)
            if name in g:                       # same label twice in one id group
                name += '#%d' % b.result_id
            g[name] = self._values(b)
        # One table per family: a family is written in chunks (one id list per element shape,
        # and per result category), so take the union of ids and columns; a value a chunk
        # does not carry is NaN.
        byfam = OrderedDict()
        for key, g in groups.items():
            byfam.setdefault('node' if key == 'node' else g['_family'], []).append(g)
        merged = OrderedDict()
        for fam, gs in byfam.items():
            if len(gs) == 1:
                merged[fam] = gs[0]
                continue
            ids = np.unique(np.concatenate([g['id'] for g in gs]))
            cols = OrderedDict(id=ids)
            for g in gs:
                rows = np.searchsorted(ids, g['id'])
                for k, v in g.items():
                    if k in ('id', '_family'):
                        continue
                    if k not in cols:
                        cols[k] = np.full(len(ids), np.nan)
                    cols[k][rows] = v
            merged[fam] = cols
        if self.mec_path:
            mesh = self.mesh
            named = OrderedDict()
            for fam, g in merged.items():
                xyz = (self._node_xyz(g['id']) if fam == 'node' else self._centroids(g['id']))
                cols = OrderedDict(id=g['id'])
                if fam != 'node':
                    cols['type'] = np.array([mesh['element_type'].get(int(e), '') for e in g['id']])
                    fam = '%s_%s' % (fam, '+'.join(sorted(set(cols['type'].tolist()))))
                cols['x'], cols['y'], cols['z'] = xyz.T
                cols.update((k, v) for k, v in g.items() if k not in ('id', '_family'))
                named[fam] = cols
            merged = named
        for g in merged.values():
            g.pop('_family', None)
        return merged

    def dataframe(self, step=-1, location='node', corners=False):
        """table() as pandas DataFrames (needs pandas)."""
        import pandas as pd
        t = self.table(step, location, corners)
        dfs = {k: pd.DataFrame(v) for k, v in t.items()}
        return dfs['node'] if location == 'node' else dfs

    def export_csv(self, out_dir, steps=None, corners=False):
        """Per step: <step>_nodes.csv, <step>_elements_<family>.csv, <step>_sam_fos.csv,
        plus summary.csv. steps: list of step keys (default all)."""
        os.makedirs(out_dir, exist_ok=True)
        rows = []
        for s in (self.steps if steps is None else [self.step(k) for k in steps]):
            tag = '%03d_S%02d_%s_I%03d' % (s.index, s.stage, re.sub(r'[^A-Za-z0-9]+', '_', s.stage_label).strip('_'),
                                           s.increment)
            for fam, cols in self.table(s, 'node').items():
                _write_cols(os.path.join(out_dir, tag + '_nodes.csv'), cols)
            for fam, cols in self.table(s, 'element', corners).items():
                _write_cols(os.path.join(out_dir, '%s_elements_%s.csv' % (tag, fam)), cols)
            if s.sam is not None:
                _write_cols(os.path.join(out_dir, tag + '_sam_fos.csv'),
                            OrderedDict(polyline=s.sam['polyline_ids'], fos=s.sam['fos']))
            row = OrderedDict(step=s.index, stage=s.stage, stage_label=s.stage_label, increment=s.increment,
                              label=s.label, factor='%g' % s.factor,
                              fos='' if s.fos is None else '%.4f' % s.fos,
                              sam_min_fos='' if s.sam is None else '%.4f' % s.sam['fos'].min())
            try:
                f = self.nodal('TOTAL TRANSLATION', s)
                k = int(np.argmax(f.values))
                row.update(max_total_translation='%.6e' % f.values[k], at_node=int(f.ids[k]))
            except KeyError:
                row.update(max_total_translation='', at_node='')
            rows.append(row)
        _write_rows(os.path.join(out_dir, 'summary.csv'), rows)
        return rows

    # ------------------------------------------------------------------ mesh / model type
    @property
    def mesh(self):
        if self._mesh is None:
            if not self.mec_path:
                raise RuntimeError('no .mec deck; pass mec= to NfxpFile for coordinates / element types')
            self._mesh = read_mec(self.mec_path)
        return self._mesh

    @property
    def model_type(self):
        """'3D', 'axisymmetric' or '2D' from the .mec element cards; without a .mec, '3D'
        if there are out-of-plane displacements, else '2D or axisymmetric'."""
        if self.mec_path:
            return self.mesh['model_type']
        try:
            f = self.nodal('T3 TRANSLATION', -1)
            return '3D' if np.any(f.values != 0) else '2D or axisymmetric'
        except KeyError:
            return 'unknown'

    def _node_xyz(self, ids):
        m = self.mesh
        k = np.searchsorted(m['node_ids'], ids)
        k = np.clip(k, 0, len(m['node_ids']) - 1)
        out = m['xyz'][k].copy()
        out[m['node_ids'][k] != ids] = np.nan
        return out

    def _centroids(self, ids):
        m = self.mesh
        if 'centroid' not in m:
            m['centroid'] = {e: self._node_xyz(np.array(c, np.int64)).mean(axis=0)
                             for e, c in m['connectivity'].items() if c}
        nan = np.full(3, np.nan)
        return np.array([m['centroid'].get(int(e), nan) for e in ids]).reshape(-1, 3)


def _write_cols(path, cols):
    names = list(cols)
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(names)
        arrays = [cols[k] for k in names]
        for i in range(len(arrays[0])):
            w.writerow([_fmt(a[i]) for a in arrays])


def _write_rows(path, rows):
    if not rows:
        return
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def _fmt(v):
    if isinstance(v, (float, np.floating)):
        return '' if np.isnan(v) else '%.6g' % v
    return v


# ---------------------------------------------------------------------- .mec deck
def _float(s):
    s = s.strip()
    try:
        return float(s)
    except ValueError:                      # Nastran short exponent: 1.5-3
        return float(re.sub(r'([0-9.])([+-]\d+)$', r'\1e\2', s))


def read_mec(path):
    """Nodes, element types / connectivity and model type from a MEC deck.
    Returns dict(node_ids (sorted), xyz, element_type {id: card}, connectivity {id: [nodes]},
    model_type)."""
    nids, xyz, etype, conn = [], [], {}, {}
    card = None

    def flush(c):
        name, f = c
        if name == 'GRID':
            nids.append(int(f[0]))
            xyz.append((_float(f[2]), _float(f[3]), _float(f[4])))
        else:
            eid = int(f[0])
            nodes = []
            for x in f[2:]:
                if not x or not re.fullmatch(r'[+-]?\d+', x):
                    break
                nodes.append(int(x))
            etype[eid], conn[eid] = name, nodes

    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if line[:1] == '$':
                continue
            fields = [x.strip() for x in line.rstrip('\r\n').split(',')]
            name = fields[0]
            if card is not None and name == '':
                card[1].extend(fields[1:9])
                continue
            if card is not None:
                flush(card)
                card = None
            if name == 'GRID' or (name[:1] == 'C' and len(fields) > 3 and fields[1].isdigit()
                                  and not name.startswith(('CORD', 'CONM', 'CMB', 'CBUSH1D'))):
                card = (name, fields[1:9])
    if card is not None:
        flush(card)
    order = np.argsort(nids)
    types = set(etype.values())
    if types & SOLID:
        mt = '3D'
    elif types & AXISYMMETRIC:
        mt = 'axisymmetric'
    else:
        z = np.array(xyz)[:, 2] if xyz else np.zeros(1)
        mt = '2D' if np.all(z == 0) else '3D'
    return dict(node_ids=np.array(nids, np.int64)[order], xyz=np.array(xyz, float).reshape(-1, 3)[order],
                element_type=etype, connectivity=conn, model_type=mt)


# ---------------------------------------------------------------------- legacy helpers
def read_nfxp(path):
    """Legacy: list of dicts (one per saved step) with all values loaded."""
    out = []
    with NfxpFile(path, find_mec=False) as r:
        for s in r.steps:
            nb = [b for b in s.blocks if b.location == 'node']
            eb = [b for b in s.blocks if b.location != 'node']
            out.append(OrderedDict(
                index=s.stage, label=s.stage_label, increment=s.step_no, factor=s.factor,
                increment_label=s.label, deleted_elements=r.deleted_elements(s), srm_fos=s.srm_fos,
                sam=s.sam, node_ids=r._ids((nb[0].ids_offset, nb[0].n)) if nb else None,
                element_ids=r._ids((eb[-1].ids_offset, eb[-1].n)) if eb else None,
                results=[dict(id=b.result_id, location=b.location, category=b.category, label=b.label,
                              element_ids=None if b.location == 'node' else r._ids((b.ids_offset, b.n)),
                              values=r._values(b)) for b in s.blocks]))
    return out


def read_mec_mesh(path):
    """Legacy: ({node: (x, y, z)}, {element: (card, [nodes])})."""
    m = read_mec(path)
    nodes = {int(i): tuple(p) for i, p in zip(m['node_ids'], m['xyz'])}
    return nodes, {e: (m['element_type'][e], c) for e, c in m['connectivity'].items()}


def export_csv(stages, out_dir, mesh=None):
    raise NotImplementedError('use NfxpFile(path).export_csv(out_dir)')


# ---------------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser(description='Read GTS NX / MEC .nfxp result files')
    ap.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    sub =ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('info', help='model type, stages, steps and result components')
    p.add_argument('nfxp')
    p.add_argument('--mec')
    p.add_argument('--step', default='-1', help='step whose components to list (default last)')
    p = sub.add_parser('export', help='CSV per step')
    p.add_argument('nfxp')
    p.add_argument('--mec')
    p.add_argument('--out')
    p.add_argument('--steps', default='last-of-stage', help="'all', 'last', 'last-of-stage' (default) or indices 0,5,9")
    p.add_argument('--corners', action='store_true', help='also export element-node (NDk-) values')
    p = sub.add_parser('history', help='one node/element result over all steps')
    p.add_argument('nfxp')
    p.add_argument('component')
    p.add_argument('id', type=int)
    p.add_argument('--element', action='store_true', help='id is an element (centroid value)')
    p.add_argument('--category')
    a = ap.parse_args()

    r = NfxpFile(a.nfxp, mec=getattr(a, 'mec', None))
    if a.cmd == 'info':
        print(r)
        print('run date:', r.run_date, '| mec:', r.mec_path or '-', '| complete:', r.complete)
        for s in r.steps:
            extra = (' FoS %.4f' % s.fos if s.fos is not None else '') + \
                    (' SAM min FoS %.4f' % s.sam['fos'].min() if s.sam is not None else '')
            print('  %4d  stage %2d %-28s incr %3d  %-24s blocks %d%s' % (
                s.index, s.stage, s.stage_label[:28], s.increment, s.label, len(s.blocks), extra))
        k = int(a.step) if re.fullmatch(r'-?\d+', a.step) else a.step
        print('components in', r.step(k))
        for loc, cat, comp in r.components(k):
            print('  %-8s %-7s %s' % (loc, cat, comp))
    elif a.cmd == 'export':
        if a.steps == 'all':
            steps = None
        elif a.steps == 'last':
            steps = [-1]
        elif a.steps == 'last-of-stage':
            steps = [r.stage_steps(st)[-1].index for st in r.stages]
        else:
            steps = [int(x) for x in a.steps.split(',')]
        out = a.out or os.path.splitext(a.nfxp)[0] + '_results'
        for row in r.export_csv(out, steps, a.corners):
            print(', '.join('%s=%s' % kv for kv in row.items()))
        print('written to', out)
    else:
        s, v = r.history(a.component, a.id, 'element' if a.element else 'node', a.category)
        print('step,stage,increment,label,value')
        for i, x in zip(s, v):
            st = r.steps[i]
            print('%d,%s,%d,%s,%.6g' % (i, st.stage_label, st.increment, st.label, x))
    r.close()


if __name__ == '__main__':
    main()
