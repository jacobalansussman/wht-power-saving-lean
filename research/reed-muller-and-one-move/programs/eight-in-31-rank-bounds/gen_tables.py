"""gen_tables.py: tables for scanrc.c (block 5, one move per direction, canonical Singer targets T_k = mult(alpha^2k)).
Writes tables5.bin (little-endian int32): nf, then mask[nf] (bit z-1 <-> direction z), wt[nf], vid[nf] (parity class),
then nv=32 and pk[32] (bit k set <-> <G, T_k> = 1 for the forms of that class)."""
import struct, numpy as np
from xbase import Field, Forms, canon_targets
import sys
m = int(sys.argv[1]) if len(sys.argv) > 1 else 5; Fd = Field(m); fo = Forms(m); n = (1 << m) - 1
allT = canon_targets(Fd, range(n))
P = np.stack([fo.parity(fo.tvec(T)) for T in allT], axis=1).astype(np.int64)
pm = (P << np.arange(n)[None, :]).sum(axis=1)
ids = {}; vid = np.zeros(fo.nf, dtype=np.int64)
ids[0] = 0
for q in range(fo.nf):
    k = int(pm[q]); ids.setdefault(k, len(ids)); vid[q] = ids[k]
assert len(ids) == n + 1
pk = [0] * (n + 1)
for k, i in ids.items(): pk[i] = k
with open("tables%d.bin" % m if m != 5 else "tables5.bin", "wb") as f:
    f.write(struct.pack("<i", fo.nf))
    f.write(fo.mask.astype("<u4").tobytes()); f.write(fo.wt.astype("<i4").tobytes()); f.write(vid.astype("<i4").tobytes())
    f.write(struct.pack("<i", n + 1)); f.write(np.array(pk, dtype="<u4").tobytes())
print("ok", fo.nf, "basis", Fd.basis)
