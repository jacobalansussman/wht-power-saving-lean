# c16_survivors.py -- shapes of a tight 5 in 14 that pass BOTH per-irreducible filters (hook balance + 105a inequality),
# exact symmetry, c(H) inside L; and for the order-24 / order-12 groups the irreducible where 105a fails and by how much.
import numpy as np, time
from lib514 import *
import c7_symscreen2 as S7
import c14_105a as C14
t0 = time.time()
Lsubs = S7.subgroups_upto_conj(LEVI)
S4_111 = gen_closure([g for g in LEVI if int(g[7]) == 7])
groups = [(H, f'order {len(H)}') for H in Lsubs if len(H) in (6, 8, 12, 24)]
for H, name in groups:
    info = C14.analyse105(H)
    n = info['n']; mov = [m for m in info['mov'] if m['exact']]; k = info['k']
    arrsz = [n // s for s in info['Tsizes']]
    def arr_rec(i, rem, acch, acc5, cnt, norb):
        if rem < 0: return
        if i == len(arrsz):
            if rem == 0: yield acch, acc5, cnt, norb
            return
        c = 0
        while c * arrsz[i] <= rem:
            yield from arr_rec(i + 1, rem - c * arrsz[i], acch + c * info['arr_h'][i], acc5 + c * info['arr5'][i], cnt + [c], norb + c)
            c += 1
    surv = set(); worst = []
    for dh, d5, acnt, aorb in arr_rec(0, 5, np.zeros(k), np.zeros(k), [], 0):
        a = '+'.join(f"{c}x{arrsz[i]}" for i, c in enumerate(acnt) if c)
        def rec(i, rem, ah, a5, cnt, norb):
            if rem < 0 or np.any(ah - dh > 1e-6): return
            if i == len(mov):
                if rem == 0 and norb >= aorb + 2 and np.all(np.abs(ah - dh) < 1e-6):
                    m = '+'.join(f"{c}x{mov[j]['size']}(T{mov[j]['T']})" for j, c in enumerate(cnt) if c)
                    gap = a5 - d5
                    if np.all(gap > -1e-6): surv.add(f'arrays {a} | moves {m}')
                    else: worst.append((float(gap.min()), int(np.argmin(gap)), f'arrays {a} | moves {m}'))
                return
            c = 0
            while c * mov[i]['size'] <= rem:
                rec(i + 1, rem - c * mov[i]['size'], ah + c * mov[i]['hook'], a5 + c * mov[i]['m105'], cnt + [c], norb + c)
                c += 1
        rec(0, 14, np.zeros(k), np.zeros(k), [], 0)
    print(f'H {name}: shapes passing hook balance AND 105a: {len(surv)}  [{time.time()-t0:.0f}s]', flush=True)
    for s_ in sorted(surv)[:12]: print('     ', s_)
    if not surv and worst:
        worst.sort(reverse=True)
        print('   best failing shapes (largest min 105a margin, irreducible index):', worst[:4])
