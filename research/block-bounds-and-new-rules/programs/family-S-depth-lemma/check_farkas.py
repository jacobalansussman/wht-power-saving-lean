# pf-logging, 13 Oct 2026.  CHECKER for farkas_familyS.json.
# Part 1 (the certificate, valid for every m): for each programme, verify that every multiplier of an inequality is
#   >= 0 and that  objective - sum lambda_i g_i  is identically a constant >= 0 (exact rational arithmetic).
#   Since every g_i >= 0 on the feasible set, objective >= 0 for every feasible point and every m.
# Part 2 (model check, small m): re-enumerate the search space of familyS_ip.py / familyS_ip2.py (own loops), map each
#   enumerated configuration to the certificate's variables, and verify every listed constraint g_i >= 0 and the
#   objective >= 0 there.  This ties the constraint list to the programmes' actual semantics.
# Usage: python3 -I check_farkas.py farkas_familyS.json [MMAX]
import json, sys, itertools
from fractions import Fraction as F
cert = json.load(open(sys.argv[1])); MMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 5
allok = True
for name, prog in cert["programs"].items():
    tot = {kk: F(v) for kk, v in prog["objective"].items()}
    okmult = True
    for cname, con in prog["constraints"].items():
        lam = F(prog["multipliers"].get(cname, "0"))
        if con.get("type") != "eq" and lam < 0: okmult = False
        for var, coef in con["form"].items(): tot[var] = tot.get(var, F(0)) - lam * coef
    resid = {kk: v for kk, v in tot.items() if v != 0}
    const = resid.pop("1", F(0))
    ok = okmult and not resid and const >= 0
    allok &= ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: multipliers nonnegative: {okmult}; objective - sum lambda g = "
          f"{dict((k_, str(v)) for k_, v in resid.items()) or 'no variable terms'}, constant {const}  ({prog['objective_meaning']})")

def holds(prog, val):
    """integer evaluation of every constraint form (all coefficients are integers) at the configuration val"""
    bad = []
    for cname, con in prog["constraints"].items():
        s = sum(int(c) * (1 if v == "1" else val.get(v, 0)) for v, c in con["form"].items())
        if (con.get("type") == "eq" and s != 0) or s < 0: bad.append(cname)
    s = sum(int(c) * (1 if v == "1" else val.get(v, 0)) for v, c in prog["objective"].items())
    if s < 0: bad.append("OBJECTIVE")
    return bad

def steiner_all(d1, d2, mind):
    """all (size, L, dist) the programmes' steiner() may realise: L <= min(d1,d2), d1+d2-2L >= mind"""
    return [(d1 + d2 - L, L, d1 + d2 - 2 * L) for L in range(0, min(d1, d2) + 1) if d1 + d2 - 2 * L >= mind]

for m in range(2, MMAX + 1):
    H = 2 * m; checked = 0; fails = 0; _OPT = {}
    for a, a2, c, r1, r2 in itertools.product(range(m + 1), range(m + 1), range(1, m + 1), range(m + 1), range(m + 1)):
        if a + c + r1 < m or a2 + c + r2 < m: continue
        dJ = max(0, m - (r1 + r2), m - (a + a2))
        # ---- familyS_ip.py: both junks cancelled by pure nodes (dx >= own a, dz >= dJ, dist >= m - other) ----
        def trees(aa, other):
            for dx in range(aa, H + 1):
                for dz in range(dJ, H + 1):
                    for (s, L, dist) in steiner_all(dx, dz, m - other):
                        yield s, dx, dz, L, dist
        T1 = list(trees(a, a2)); T2 = list(trees(a2, a))
        # check every pairing would be large; check per-side minimal-size representatives plus a sample of others
        for t1 in T1[:: max(1, len(T1) // 40)] + [min(T1)]:
            for t2 in T2[:: max(1, len(T2) // 40)] + [min(T2)]:
                val = dict(a=a, a2=a2, c=c, r1=r1, r2=r2, dJ=dJ, m=m, s1=t1[0], dpar1=t1[1], dz1=t1[2], L1=t1[3], dist1=t1[4],
                           s2=t2[0], dpar2=t2[1], dz2=t2[2], L2=t2[3], dist2=t2[4])
                val["n"] = val["s1"] + val["s2"] + c + r1 + r2
                b = holds(cert["programs"]["ip1"], val); checked += 1; fails += bool(b)
        # ---- familyS_ip2.py: J1 on T2 (parent depth >= a2, K uses a), J2 on T1 (parent depth >= a, K uses a2) ----
        def options(d_par, a_other):
            key = (d_par, a_other, dJ)
            if key in _OPT: return _OPT[key]
            out = []
            for dpar in range(d_par, H + 1):
                for dz in range(dJ, H + 1):                           # pure cover
                    for dist in range(0, dpar + dz + 1):
                        if dist + a_other < m: continue
                        for (s, L, dd) in steiner_all(dpar, dz, dist):
                            out.append((s, 0, dpar, dz, L, dd))
                for q in range(1, H + 1):                            # lineage cover
                    for dyf in range(0, H + 1):
                        if q + dyf < dJ: continue
                        for dist in range(0, dpar + dyf + 1):
                            if q + dist + a_other < m: continue
                            for (s, L, dd) in steiner_all(dpar, dyf, dist):
                                out.append((s, q, dpar, dyf, L, dd))
            _OPT[key] = out
            return out
        o1 = options(a2, a); o2 = options(a, a2)
        s1o = o2[:: max(1, len(o2) // 60)] + [min(o2)]; s2o = o1[:: max(1, len(o1) // 60)] + [min(o1)]
        for (s1, q2, dp1, dz1, L1, di1) in s1o:
            for (s2, q1, dp2, dz2, L2, di2) in s2o:
                lin = max(q1, q2) if (q1 and q2) else q1 + q2
                val = dict(a=a, a2=a2, c=c, r1=r1, r2=r2, dJ=dJ, m=m, s1=s1, dpar1=dp1, dz1=dz1, L1=L1, dist1=di1,
                           s2=s2, dpar2=dp2, dz2=dz2, L2=L2, dist2=di2, q1=q1, q2=q2, lin=lin)
                val["n"] = s1 + s2 + c + r1 + r2 + lin
                b = holds(cert["programs"]["ip2"], val); checked += 1; fails += bool(b)
                if b and fails < 5: print("   model-check failure", m, val, b)
    allok &= fails == 0
    print(f"[{'PASS' if fails == 0 else 'FAIL'}] model check m={m}: {checked} enumerated configurations, constraint violations {fails}")
print("RESULT:", "PASS" if allok else "FAIL")
sys.exit(0 if allok else 1)
