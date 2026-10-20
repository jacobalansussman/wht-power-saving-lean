# skeptic-1 (wht44/ideas), 16 Oct 2026. Labelled certificates on the coupled relaxation 56+105a (a2-gct's pair).
# P0-invariant subspaces Y of V56 + V105a: each multiplicity-one P0-type taken 0 or whole; the shared type 21
# (one copy in each irreducible) taken 0, whole, or a line l in C^2 of the multiplicity space (a "label"),
# real angles and complex phases. For each Y: cap of h, demands of K (full) and K2 (two-row). Theorem R then refuses
# W_K + W_R in M iff M cap < W_K dK + W_R dR. Question: do labels refuse the Theorem-T-excluded tight mixed shapes
# 4+1 in 13, 3+2 in 12, 2+3 in 11 (allowed by every single cut)? Calibration: 5 in 13 and 1+4 in 10 refused,
# nothing beats 2.8 for full arrays, as a2-gct found on 56+70.
import numpy as np
np.seterr(all='ignore')
import importlib.util, sys
spec = importlib.util.spec_from_file_location('s1', __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 's1_depth_cuts.py'))
s1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s1)


def cdim(B):
    if B.shape[1] == 0: return 0
    s = np.linalg.svd(B, compute_uv=False)
    return int((s > 1e-8 * max(1.0, s[0])).sum())


def main():
    n, h, k, k2, T = s1.load([8, 11])
    print('types (dim, mult, irreducibles):', [(t['d'], len(t['slots']), t['js']) for t in T])
    single = [i for i, t in enumerate(T) if len(t['slots']) == 1]
    shared = [i for i, t in enumerate(T) if len(t['slots']) == 2]
    assert len(shared) == 1
    sh = T[shared[0]]
    I = np.eye(n)
    labels = [('0', None), ('all', 'all')]
    for th in np.linspace(0, np.pi, 37)[1:-1]:
        for ph in (0.0, np.pi / 2, np.pi / 3):
            labels.append((f'line(th={th:.3f},ph={ph:.2f})', (np.cos(th), np.sin(th) * np.exp(1j * ph))))
    labels.append(('line(1,0)', (1.0, 0.0))); labels.append(('line(0,1)', (0.0, 1.0)))
    shapes = [(5, 0, 14), (4, 1, 13), (3, 2, 12), (2, 3, 11), (1, 4, 10), (5, 0, 13)]
    best = {s: None for s in shapes}
    best_ratio = (0, None)
    hc = h.astype(complex); kc = k.astype(complex); k2c = k2.astype(complex)
    cnt = 0
    for mask in range(2 ** len(single)):
        base = [I[:, T[single[i]]['slots'][0]] for i in range(len(single)) if mask >> i & 1]
        for lname, l in labels:
            cols = list(base)
            if l == 'all':
                cols += [I[:, sh['slots'][0]], I[:, sh['slots'][1]]]
            elif l is not None:
                cols.append(l[0] * I[:, sh['slots'][0]] + l[1] * I[:, sh['slots'][1]])
            if not cols: continue
            Y = np.hstack(cols).astype(complex); dy = cdim(Y)
            if dy in (0, n): continue
            cap = cdim(np.hstack([Y, hc @ Y])) - dy
            dK = cdim(np.hstack([Y, kc.T @ Y])) - dy
            dR = cdim(np.hstack([Y, k2c.T @ Y])) - dy
            cnt += 1
            if cap > 0 and dK / cap > best_ratio[0]: best_ratio = (dK / cap, (mask, lname, dy, cap, dR, dK))
            for (WK, WR, M) in shapes:
                sl = M * cap - (WK * dK + WR * dR)
                if best[(WK, WR, M)] is None or sl < best[(WK, WR, M)][0]:
                    best[(WK, WR, M)] = (sl, mask, lname, dy, cap, dR, dK)
    print(f'{cnt} certificates scanned (single types 0/all, shared 21 as 0, all, or {len(labels) - 2} labels)')
    print(f'best full-array ratio dK/cap = {best_ratio[0]:.4f} at {best_ratio[1]}')
    for s in shapes:
        b = best[s]
        print(f'{s[0]}+{s[1]} in {s[2]}: least slack {b[0]} ({"REFUSED" if b[0] < 0 else "allowed"}) at '
              f'mask={b[1]} label={b[2]} dim={b[3]} cap={b[4]} dR={b[5]} dK={b[6]}')


if __name__ == '__main__':
    main()
