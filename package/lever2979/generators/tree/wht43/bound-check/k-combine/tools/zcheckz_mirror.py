"""k-combine: exact mirror of GS.zCheckZ (Work/KMX/ZCheckZ.lean), the function of zdesignZ.py."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zplan import Run, No, guard
from zmirror import okAZK, okBZK

def zCheck(c, p, gl, rows):
    """exact mirror of GS.zCheckZ (zplan.zCheck with the mutable-source kinds)"""
    Rn = Run(c, p, 0, len(gl))
    lo = 0
    for r in gl:
        guard(lo <= r and r + 1 <= Rn.nr, "incK")
        lo = r + 1
    guard(Rn.us + 1 <= 2 ** (Rn.sw - 1 if Rn.sw >= 1 else 0), "us >= 2^(sw-1)")
    guard(Rn.us + 1 <= 2147483648, "us >= 2^31")
    Rn.gInitL(Rn.us, gl)
    Rn.items(okAZK, okBZK)
    g = set(gl)
    for r in rows:
        if r in g:
            continue
        guard(Rn.key_ok(r), "zFin: register outside the trie")
        x = Rn.tr.get(r, 0)
        if ((x % Rn.m2) >> 64) != ((x >> Rn.K) >> 64):
            raise No("zFin: register %d holds something" % r)
    return Rn.adds
