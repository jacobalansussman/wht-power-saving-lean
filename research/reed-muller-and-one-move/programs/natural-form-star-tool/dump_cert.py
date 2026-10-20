"""dump_cert.py <case> <F> [dmax]: run rcdp on one form and write its certificate (raw rcdp text) to stdout."""
import sys, os, struct, subprocess
import numpy as np
from rc_run import build, HERE
from xbase import Field, Forms, System
case = sys.argv[1]; F = int(sys.argv[2]); dmax = int(sys.argv[3]) if len(sys.argv) > 3 else 20
m, targets, counts, _ = build(case)
Fd = Field(m); fo = Forms(m); sy = System(fo, targets, counts)
nd = len(fo.dirs)
present = sum(1 << i for i, z in enumerate(fo.dirs) if sy.counts.get(z, 0) > 0)
masks = (fo.mask & present).astype(np.int64)
inp = os.path.join(HERE, "tmp_dc_%d.bin" % os.getpid())
with open(inp, "wb") as f:
    f.write(struct.pack("<5i", m, nd, fo.nf, 1, dmax))
    f.write(np.array([sy.counts.get(z, 0) for z in fo.dirs], dtype="<i4").tobytes())
    f.write(masks.astype("<u4").tobytes()); f.write(sy.s.astype("<i4").tobytes())
    f.write(np.array([F], dtype="<i4").tobytes())
out = subprocess.run([os.path.join(HERE, "rcdp"), inp], capture_output=True, text=True).stdout
os.remove(inp)
print("# case", case, "basis", Fd.basis, "poly", bin(5 and 0b100101))
print(out)
