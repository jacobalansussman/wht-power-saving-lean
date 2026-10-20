"""h-L1: exact check of the L1 -> hosting seam in our engine's semantics (frames `core`, adapters `kg`), |H| = 9.

A host role has the frame of label V (base G_h, coordinates S_h) in ITS window X_h = ofConj d_h q_h, and its idle
climb goes to V' (same base, coordinates T_h ⊇ S_h).  One or two invocations, each in its OWN window X_q = ofConj d_q q_q
with an entrance label E (base G_e, dim e), host a slot residual on the host:

  before the call  : free adapter A = kg G_q' * kg tau * kg G_h'^-1   (G_q' = X_q.liftG G_e, G_h' = X_h.liftG G_h,
                                                                       tau a coordinate permutation of alpha with
                                                                       tau(missing(X_q, E)) = S', S' ⊆ T_h \\ S_h pending)
  the call (L1)    : resW = kg G_q' * blockMat(unit vectors of missing(X_q, E)) * kg G_q'^-1   (on whatever the slot holds)
  after the call   : A^-1 = kg G_h' * kg tau^-1 * kg G_q'^-1
  remaining climb  : kg G_h' * blockMat(unit vectors of T_h \\ (S_h ∪ S'_1 ∪ S'_2)) * kg G_h'^-1

Claim (hub PROOF.md §2 in our engine): the product equals the host's full idle climb, i.e. it takes core(G_h', S_h) to
core(G_h', T_h) exactly.  Controls that must fail: S' meeting S_h (already applied), two residuals on overlapping S',
tau not carrying S' onto the missing set, adapter without kg tau.
"""
import sys, json, random, time
sys.path.insert(0, '<verify>/checks/wht43/bound-check/v-host-0/toy')
from l1_check import Space, GI, rand_orthogonal, rand_invertible, lift_table, np

def coord_perm_table(n, tau):
    """APerm of F_2^n moving coordinate k to tau[k]"""
    T = np.zeros(1 << n, dtype=np.int64)
    for x in range(1 << n):
        y = 0
        for k in range(n):
            if x >> k & 1:
                y |= 1 << tau[k]
        T[x] = y
    return T

def inv_table(T):
    I = np.empty_like(T); I[T] = np.arange(len(T)); return I

def bits(mask, n):
    return [k for k in range(n) if mask >> k & 1]

def run(h, b, rng, two=True, control=None):
    n = h + b; sp = Space(n)
    # host window and labels V ⊆ V' in the host's base
    d_h = rand_orthogonal(n, rng); q_h = sum(1 << (h + j) for j in range(b) if rng.random() < 0.5)
    G_host = rand_invertible(h, rng)
    Gh = lift_table(h, b, d_h, G_host); Ghi = inv_table(Gh)
    sV = set(rng.sample(range(h), rng.randint(0, 1 if two else 2)))   # host label before its idle climb
    sVp = set(range(h))                                                # ... and after (full: the deepest idle climb)
    S_h = sum(1 << k for k in sV) | q_h; T_h = sum(1 << k for k in sVp) | q_h
    pending = bits(T_h & ~S_h, n)
    # invocations: own windows, entrance labels of the shapes that occur (e in {5,4,2,1,0} -> widths 4,5,7,8,9)
    calls = []
    free = list(pending); rng.shuffle(free)
    for c in range(2 if two else 1):
        widths_ok = [w for w in (4, 5, 7, 8, 9) if w <= len(free)]
        if not widths_ok:
            break
        if two and c == 0 and len(free) >= 8:
            widths_ok = [w for w in widths_ok if w <= len(free) - 4]    # leave room for a second hosted residual
        w = rng.choice(widths_ok); e = h - w
        d_q = rand_orthogonal(n, rng); q_q = sum(1 << (h + j) for j in range(b) if rng.random() < 0.5)
        G_e = rand_invertible(h, rng); s_e = set(rng.sample(range(h), e))
        Gq = lift_table(h, b, d_q, G_e); Gqi = inv_table(Gq)
        miss = bits((((1 << h) - 1) | q_q) & ~(sum(1 << k for k in s_e) | q_q), n)   # = inl(univ \\ s_e)
        assert len(miss) == w
        Sp = free[:w]; free = free[w:]
        if control == 'applied_overlap' and c == 0 and sV:
            Sp = Sp[:-1] + [min(sV)]                                   # one coordinate already applied on the host
        if control == 'double_book' and c == 1:
            Sp = Sp[:-1] + [calls[0]['Sp'][0]]                         # overlaps the first call's coordinates
        # coordinate permutation tau with tau(miss[i]) = Sp[i], completed arbitrarily.  Convention check: for
        # P = kg tau = permMat tau, P^-1 * pk(M) * P = pk(tau(M)) (tau acting on coordinate sets), so tau must carry
        # the slot's missing set ONTO the hosted pending set S' (a first version had it backwards and failed).
        tau = [None] * n
        for a_, b_ in zip(miss, Sp):
            tau[a_] = b_
        rest_src = [k for k in range(n) if tau[k] is None]; rest_dst = [k for k in range(n) if k not in Sp]
        rng.shuffle(rest_dst)
        for a_, b_ in zip(rest_src, rest_dst):
            tau[a_] = b_
        if control == 'wrong_tau' and c == 0:
            tau = list(range(n)); rng.shuffle(tau)
        Tt = coord_perm_table(n, tau); Tti = inv_table(Tt)
        calls.append(dict(Gq=Gq, Gqi=Gqi, miss=miss, Sp=Sp, Tt=Tt, Tti=Tti, w=w))
    used = set(k for c in calls for k in c['Sp'])
    rem = [k for k in pending if k not in used]
    # run on the frame (data is irrelevant: everything acts on the left)
    M = sp.core_scaled(Gh, S_h)
    scale = 0
    for c in calls:
        # A = kg Gq * kg tau * kg Gh^-1   (applied right to left)
        M = sp.kg_left(Ghi, M)
        if control != 'no_tau':
            M = sp.kg_left(c['Tt'], M)
        M = sp.kg_left(c['Gq'], M)
        # the call: resW = kg Gq * blockMat(miss) * kg Gq^-1
        M = sp.kg_left(c['Gqi'], M); M = sp.block2_left([1 << k for k in c['miss']], M); M = sp.kg_left(c['Gq'], M)
        scale += c['w']
        # A^-1 = kg Gh * kg tau^-1 * kg Gq^-1
        M = sp.kg_left(c['Gqi'], M)
        if control != 'no_tau':
            M = sp.kg_left(c['Tti'], M)
        M = sp.kg_left(Gh, M)
    # remaining climb
    M = sp.kg_left(Ghi, M); M = sp.block2_left([1 << k for k in rem], M); M = sp.kg_left(Gh, M)
    scale += len(rem)
    target = sp.core_scaled(Gh, T_h).scal(1 << scale)
    return M.eq(target), dict(h=h, b=b, pending=len(pending), widths=[c['w'] for c in calls], remaining=len(rem),
                             control=control)

if __name__ == '__main__':
    rng = random.Random(20261015)
    res = []; ok = True
    for b in (1, 2):
        for t in range(6):
            r, info = run(9, b, rng, two=True); res.append((r, info)); ok &= r
            print('hosted == full climb:', r, info, flush=True)
        for t in range(3):
            r, info = run(9, b, rng, two=False); res.append((r, info)); ok &= r
            print('hosted == full climb:', r, info, flush=True)
        for ctl in ('applied_overlap', 'double_book', 'wrong_tau', 'no_tau'):
            r, info = run(9, b, rng, two=True, control=ctl)
            applicable = not (ctl == 'applied_overlap' and info['pending'] == 9) and not (ctl == 'double_book' and len(info['widths']) < 2)
            res.append((not r, info)); ok &= (not r) or not applicable
            print(f'control {ctl} (must fail): equality={r} applicable={applicable}', info, flush=True)
    print(json.dumps(dict(all_pass=bool(ok), runs=len(res))))
    print('ALL SEAM CHECKS PASS' if ok else 'SOME SEAM CHECK FAILED')
    sys.exit(0 if ok else 1)
