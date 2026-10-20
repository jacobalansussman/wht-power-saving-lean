#!/usr/bin/env python3
"""cx-mirror / twist3.py (own code, stdlib).  A wider sample of twisted mirrors.  qcheck.py shows every twisted
mirror is an involution A = A^T with <e, A e> = 0 on the even space, so N is symmetric and a drop-free x chain inside a
support is an independent set of {N = 1}.  Sample: random products of 1 to 6 Eichler maps (all such products that pass
are kept), exact maximum independent set per support.  Reports the distribution of the worst-support chain."""
import os, json, random, importlib.util, time
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('tw', os.path.join(here, 'twist_lib.py')); tw = importlib.util.module_from_spec(spec); spec.loader.exec_module(tw)
U, Sk, v, dot, apply, compose, eichler, orth, ok, I9, mis = tw.U, tw.Sk, tw.v, tw.dot, tw.apply, tw.compose, tw.eichler, tw.orth, tw.ok, tw.I9, tw.mis
rnd = random.Random(23); even = [w for w in range(1, 512) if bin(w).count('1') % 2 == 0]
seen = set(); hist = {}; best = None; t0 = time.time(); n = 0
while time.time() - t0 < 80:
    A = I9
    for _ in range(rnd.randint(1, 6)):
        w = rnd.choice(even); z = rnd.choice([z for z in even if dot(w, z) == 0 and z != w])
        A = compose(eichler(w, z), A)
    key = tuple(A)
    if key in seen or not ok(A): continue
    seen.add(key)
    assert all(dot(A[i], 1 << j) == dot(A[j], 1 << i) for i in range(9) for j in range(9)), 'a passing A that is not symmetric'
    adj = {a: {b for b in range(v) if a != b and not (1 ^ dot(U[b], apply(A, U[a])) == 0)} for a in range(v)}
    c = [len(mis(S, adj)) for S in Sk]; n += 1
    hist[min(c)] = hist.get(min(c), 0) + 1
    if best is None or min(c) > min(best[1]): best = (A, c)
print('distinct passing A sampled:', n, '| all symmetric (asserted) | worst-support chain histogram', dict(sorted(hist.items())),
      '| best', best[1] if best else None)
json.dump(dict(sampled=n, histogram=hist, best=dict(A=best[0], chains=best[1]) if best else None), open(os.path.join(here, 'twist3.json'), 'w'), indent=1)
