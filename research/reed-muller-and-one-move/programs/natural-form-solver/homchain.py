"""homchain.py (rm2-solve, 15 Oct): deletion CHAIN from a solved plain program, one move at a time.

usage: homchain.py SPEC ORDER NSUB MAXDEL CANDS SEED TLIMIT [X0.npy] [EXTRA]
  SPEC  'k:0,2,4,8,...' canonical Singer targets T_k (plain frames alpha^k f_j), or 'I:W' (slice, all targets I,
        W copies of the standard frame).  Block from env HM (default 5).
  ORDER time order of the plain moves (randN | lex | keep (grouped by array)).
  EXTRA number of redundant extra moves (a copy of a random direction already present, inserted at a random time)
        -- the machinery control: a redundant move must be removable.
Stage 0: solve the plain W in mW (+EXTRA) from random starts (or load X0, which must match SPEC and ORDER).
Stage k: delete one more move. Candidates: every move (env ANYDEL=0: only moves whose direction still occurs at
  least twice) whose deletion passes span + deletion; ranked by the exact
  tests of what remains (span, deletion, then fewest tight forms) and by the move's size in the balanced gauge
  |in_i| |out_i| (smallest first); try up to CANDS of them. For each: penalty path mu |in_i|^2 + |out_i|^2 on the
  move (mu 1e-3 .. 1e4, 10 steps x 15 LM iterations on the subset), drop it, polish (LM, 250 iterations on the
  subset), then evaluate on ALL forms; success = ALL residual < 1e-12 after 4 more LM iterations on all forms.
On success the chain continues from the new program; when every candidate fails it stops ("STUCK").
By Theorem B(a) any directions are allowed: a multi-target program with Singer targets is a no-phase unit with the
same W and M. The final shape W in M is reported with its saving ln(mW/M)/ln m.
"""
import sys, os, time, warnings, random
import numpy as np
warnings.filterwarnings("ignore")
from inst import canonical, slice_I, Instance
from glm import Model, lm, random_start

a = sys.argv[1:]
spec, order = a[0], a[1]; nsub = int(a[2]); maxdel = int(a[3]); ncand = int(a[4]); seed = int(a[5])
tlimit = float(a[6]); x0file = a[7] if len(a) > 7 and a[7] != "-" else None
extra = int(a[8]) if len(a) > 8 else 0
MM = int(os.environ.get("HM", "5"))
ANYDEL = int(os.environ.get("ANYDEL", "1"))   # 1: any move may be deleted (exact tests decide); 0: duplicates only
kind, rest = spec.split(":")
if kind == "k":
    ks = [int(v) for v in rest.split(",")]
    inst, fd = canonical(MM, ks, dirs="frames", order=order, seed=seed)
elif kind == "f":   # 'f:k1,k2,../c1,c2,..': array l uses frame number c_l (0..5) of its target (packing2 order)
    from packing2 import all_frames
    kpart, cpart = rest.split("/")
    ks = [int(v) for v in kpart.split(",")]; ch = [int(v) for v in cpart.split(",")]
    d0 = [z for k, c_ in zip(ks, ch) for z in sorted(sorted(F) for F in all_frames(k))[c_]]
    inst, fd = canonical(MM, ks, dirs=d0, order=order, seed=seed)
else:
    inst, fd = slice_I(MM, int(rest), dirs="frames", order=order, seed=seed)
tag = f"hc {spec} {order} s{seed}"
rng = np.random.default_rng(seed); pr = random.Random(seed)
T, dirs = inst.T, list(inst.dirs)
if extra:
    for e in range(extra):
        z = pr.choice(dirs); dirs.insert(pr.randrange(len(dirs) + 1), z)
W = inst.W
t0 = time.time()


def model_for(dirs):
    it = Instance(MM, T, dirs)
    md = Model(it.t, it.p, it.W, it.M)
    idx = np.sort(rng.choice(it.Nf, nsub, replace=False)) if 0 < nsub < it.Nf else None
    return it, md, idx


def log(s):
    print(f"[{tag}] {s} t {time.time() - t0:.0f}s", flush=True)


def full_check(md, x, its=4):
    c = md.cost(x)
    if c < 1e-6:
        x, c, _ = lm(md, x, iters=its, verbose=0)
    return x, c


it, md, idx = model_for(dirs)
log(f"plain W={W} M={it.M} dirs={dirs} exact={it.exact_tests()}")
if x0file and not extra:
    x = np.load(x0file); c = md.cost(x)
    log(f"loaded {x0file}: ALL {c:.3e}")
else:
    for att in range(6):
        x, c, h = lm(md, random_start(md, rng, 1, 0.3), iters=400, idx=idx, verbose=0, tlimit=1500)
        call = md.cost(x)
        log(f"plain solve attempt {att}: subset {c:.3e} ALL {call:.3e} its {len(h)}")
        if call < 1e-14: break
x, c = full_check(md, x)
log(f"STAGE 0 W={W} M={md.M} ALL {c:.3e}")
if c > 1e-12:
    log("plain program not solved; stop"); sys.exit(0)
np.save(f"hc_{spec.replace(':', '').replace(',', '-').replace('/', '_')}_{order}_s{seed}_M{md.M}.npy", x)


def balanced_size(md, x, i):
    A, B, H = md.unpack(x)
    nin = np.sqrt(np.sum(np.abs(A[:, i]) ** 2) + np.sum(np.abs(H[:i, i]) ** 2))
    nout = np.sqrt(np.sum(np.abs(B[:, i]) ** 2) + np.sum(np.abs(H[i, i + 1:]) ** 2))
    return nin * nout


def move_mask(md, i):
    M = md.M; nA = md.nA; mask = np.zeros(md.nvar)
    for c_ in range(md.W):
        mask[c_ * M + i] = 1; mask[nA + c_ * M + i] = 1
    iu0, iu1 = md.iu
    mask[2 * nA + np.nonzero((iu0 == i) | (iu1 == i))[0]] = 1
    return mask


mus = np.geomspace(1e-3, 1e4, 10)
for stage in range(1, maxdel + 1):
    if time.time() - t0 > tlimit:
        log("time limit"); break
    M = md.M
    cnt = {}
    for z in dirs: cnt[z] = cnt.get(z, 0) + 1
    cands = []
    for i in range(M):
        if cnt[dirs[i]] < 2 and not ANYDEL: continue
        rest_d = dirs[:i] + dirs[i + 1:]
        et = Instance(MM, T, rest_d).exact_tests()
        if not et["span_ok"] or et["violations"] > 0: continue
        cands.append((et["tight"], balanced_size(md, x, i) * (1 + 0.05 * pr.random()), i))
    cands.sort()
    if not cands:
        log(f"STUCK at W={W} M={M}: no candidate passes span + deletion"); break
    log(f"stage {stage}: {len(cands)} candidates; first {[(c_[2], dirs[c_[2]], c_[0], round(c_[1], 3)) for c_ in cands[:ncand]]}")
    done = False; tried = []
    for (tight, sz, i) in cands[:ncand]:
        mask = move_mask(md, i)
        xp = x.copy()
        for mu in mus:
            md.pen = mu * mask
            xp, cp, _ = lm(md, xp, iters=15, idx=idx, verbose=0, tlimit=1200)
            md.pen = None
        res = md.cost(xp, idx); pn = float((mask * np.abs(xp) ** 2).sum())
        keep = [j for j in range(M) if j != i]
        A, B, H = md.unpack(xp)
        nd = [dirs[j] for j in keep]
        it2, md2, idx2 = model_for(nd)
        x2 = md2.pack(A[:, keep], B[:, keep], H[np.ix_(keep, keep)])
        x2, c2, h2 = lm(md2, x2, iters=250, idx=idx2, verbose=0, tlimit=1800)
        x2, call = full_check(md2, x2)
        tried.append((i, dirs[i], res, pn, c2, call))
        log(f"  delete move {i} (dir {dirs[i]}, tight {tight}): penalised residual {res:.3e} size {pn:.2e}; "
            f"polished subset {c2:.3e} ALL {call:.3e}")
        if call < 1e-12:
            x, md, idx, dirs = x2, md2, idx2, nd
            s_ = np.log(MM * W / md.M) / np.log(MM)
            log(f"REACHED W={W} M={md.M} (saving {s_:.4f}) ALL {call:.3e} exact={it2.exact_tests()}")
            np.save(f"hc_{spec.replace(':', '').replace(',', '-').replace('/', '_')}_{order}_s{seed}_M{md.M}.npy", x)
            done = True; break
    if not done:
        log(f"STUCK at W={W} M={M}: all {len(tried)} candidates failed; best ALL {min(t_[5] for t_ in tried):.3e}")
        break
log(f"END W={W} M={md.M} dirs={dirs}")
