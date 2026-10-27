"""Writes small synthetic .nfxp / .mec files in the layout midas_nfxp reads, so the
tests need no real (and usually confidential) model files."""
import struct

import numpy as np

STAGE, INCREMENT, RESULT, ECSYS, SAM, SRM = 10017, 13, 14, 10003, 10025, 10032
NODE, ELEMENT, CORNER = 0, 1, 2
OTHER, DISP, FORCE, STRESS, STRAIN = 0, 1, 3, 4, 5


def _label(text):
    b = text.encode('utf-8')[:79]
    return b + b'\x00' * (80 - len(b))


class Writer(object):
    def __init__(self, date='2026/10/27 12:00'):
        head = bytearray(0x52)
        head[0x10:0x19] = b'MecSolver'
        d = date.encode('latin1')[:16]
        head[0x2c:0x2c + len(d)] = d
        self.buf = bytearray(head)

    def stage(self, number, label, analysis=0, last=False):
        self.buf += struct.pack('<7i', STAGE, 0, -1 if last else 0, 0, 12, number, analysis) + _label(label)
        return self

    def increment(self, step_no, factor, label, index, analysis=0):
        self.buf += struct.pack('<iifi', INCREMENT, step_no, factor, analysis) + _label(label) \
            + struct.pack('<i', index)
        return self

    def ids(self, kind, ids):
        ids = np.asarray(ids, '<i4')
        self.buf += struct.pack('<2i', kind, len(ids)) + ids.tobytes()
        return self

    def deleted(self, ids):
        return self.ids(4, ids)

    def ecsys(self, n):
        self.buf += struct.pack('<2i', ECSYS, n) + np.zeros(6 * n, '<f4').tobytes()
        return self

    def result(self, result_id, location, category, label, values):
        h = [0] * 18
        h[0], h[1], h[2], h[3] = RESULT, result_id, location, category
        self.buf += struct.pack('<18i', *h) + _label(label) + np.asarray(values, '<f4').tobytes()
        return self

    def srm(self, fos):
        self.buf += struct.pack('<2i', SRM, 0) + struct.pack('<d', fos)
        return self

    def sam(self, polylines, fos):
        n = len(polylines)
        self.buf += struct.pack('<6i', SAM, 0, 7, 0, 0, n) + np.asarray(polylines, '<i4').tobytes() \
            + np.asarray(fos, '<f8').tobytes()
        return self

    def raw(self, b):
        self.buf += b
        return self

    def eof(self):
        self.buf += struct.pack('<6i', STAGE, 0, -1, -1, -1, 0)
        return self

    def save(self, path):
        with open(str(path), 'wb') as f:
            f.write(bytes(self.buf))
        return path


def write_mec(path, nodes, elements):
    """nodes {id: (x, y, z)}, elements {id: (card, [node ids])}."""
    lines = ['$$ synthetic test deck']
    for i, (x, y, z) in sorted(nodes.items()):
        lines.append('GRID,%d,,%r,%r,%r' % (i, x, y, z))
    for e, (card, conn) in sorted(elements.items()):
        f = [card, str(e), '1'] + [str(c) for c in conn]
        lines.append(','.join(f[:9]))
        rest = f[9:]
        while rest:
            lines.append('       ,' + ','.join(rest[:8]))
            rest = rest[8:]
    with open(str(path), 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    return path
