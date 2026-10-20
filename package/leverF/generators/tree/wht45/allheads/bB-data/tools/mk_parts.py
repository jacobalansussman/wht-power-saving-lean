"""p-hostck (wht45/lever, 19 Oct 2026): a host table (lg-a2a4 format v1, as written by hub_to_host.py / p-own) -> the
Lean part modules of its kernel host check (p-a3's checker Work/XHostT/KCheck.lean, one `decide +kernel` theorem per
host kind), the range glue of each part (Work/XHostT/KGlue.lean) and the whole-table module (Work/XHostT/KWhole.lean).

Built on p-a3's tools/mk_kdata.py (same Calls data and the same packed streams, so the same decoders); new here: the
count tree `cntT` (the claimed record count of every host kind, checked by the kernel in each part), the claims
`= some (cK κ)`, one range theorem `rng<j> : KOK .. a b` per part, the Whole module (balanced tree of the ranges, the
global count fact, the order fact, `whole`), and negative controls.

The generator REFUSES a table that fails any check (it re-runs every check in Python with the checker's rules before
writing anything): exit code 3 and a message naming the first failing host kind.

Negative controls (`--perturb load|win`): one entry `hrows[k][q]` of the table is changed (`load`: moved to another bin
of the same host kind, so two loads go wrong; `win`: moved to the same bin of another replica whose window misses the
call).  Without `--emit-neg` the generator must refuse it.  With `--emit-neg` it writes, under TAG, the perturbed Calls
module and, for every failing kind, NegSome<κ>.lean (claims the check passes: must FAIL to compile) and
NegNone<κ>.lean (claims the checker returns `none`: must compile).

usage: python3 -I mk_parts.py HOSTTAB.json.gz TAG KINDS_PER_PART OUTDIR [--perturb load|win] [--emit-neg]
OUTDIR receives Calls.lean, Order.lean, P<j>.lean, Whole.lean, MANIFEST.json (module Work.XHostT.Data.<TAG>.*)."""
import sys, json, gzip, os, collections, hashlib

args = [a for a in sys.argv[1:] if not a.startswith('--')]
src, tag, kpp, outdir = args[0], args[1], int(args[2]), args[3]
perturb = None
if '--perturb' in sys.argv:
    perturb = sys.argv[sys.argv.index('--perturb') + 1]
    args.remove(perturb)
    assert perturb in ('load', 'win')
emit_neg = '--emit-neg' in sys.argv
raw = gzip.open(src).read()
d = json.loads(raw)
T, R, NH = d['T'], d['R'], d['NH']
NC = 9 * T
NK = 9 * T
hr = [list(r) for r in d['hrows']]
time, disp, rem = d['time'], d['disp'], list(d['hrem'])
res = {int(p): v for p, v in d['res'].items()}
KZ = [0, 2, 0, 4, 0, 6, 3, 5, 7]
KK = [2, 1, 4, 1, 6, 1, 6, 4, 2]
RK = [16, 8, 32, 8, 48, 8, 48, 32, 16]
assert d.get('kz', KZ) == KZ and d.get('kk', KK) == KK and d.get('ranks', RK) == RK, 'kinds differ from XHost.kz/kK/rankOf'


def refuse(msg):
    print('REFUSED:', msg)
    sys.exit(3)


# encodability (the checker's field widths)
if not (NC <= 4096 and R <= 1024 and NH == 2016 * T and len(hr) == NC and len(time) == NC and len(disp) == NC):
    refuse('shape: need 9T calls <= 4096, R <= 1024, NH = 2016 T')
if not all(len(r) == R and all(0 <= e < NH for e in r) for r in hr):
    refuse('every slot must be hosted (nb = 0) at a host index below NH')
if not (all(0 <= x <= 15 for p in res for x in res[p]) and set(disp) <= set(res) and all(len(res[p]) == R for p in res)):
    refuse('residuals must be 0..15, R per provider, every call dispatched to a listed provider')
if not (len(rem) == NH and all(0 <= x <= 63 for x in rem)):
    refuse('remainders must be 0..63, one per host index')
if not (sorted(time) == list(range(NC)) and all(time[k] < time[k + 1] for k in range(NC) if k % 9 != 8)):
    refuse('global order: times must be a permutation of 0..NC-1, increasing along each replica')


def win(kap):
    rh, c = divmod(kap, 9)
    lo = 0 if KZ[c] == 0 else time[9 * rh + KZ[c] - 1] + 1
    hi = time[9 * rh + KZ[c] + KK[c] - 1]
    return lo, hi


perturbed = None
if perturb:
    done = False
    for k in range(NC):
        for q, e in enumerate(hr[k]):
            kap, b = divmod(e, 224)
            if perturb == 'load':
                e2 = 224 * kap + (b + 1) % 224
                if e2 not in hr[k] and res[disp[k]][q] > 0:
                    perturbed = dict(mode='load', k=k, q=q, old=e, new=e2)
                    hr[k][q] = e2; done = True
            else:
                c = kap % 9
                for o in range(T):
                    lo, hi = win(9 * o + c)
                    e2 = 2016 * o + e % 2016
                    if not (lo <= time[k] <= hi) and e2 not in hr[k]:
                        perturbed = dict(mode='win', k=k, q=q, old=e, new=e2)
                        hr[k][q] = e2; done = True; break
            if done: break
        if done: break
    assert done, 'no entry to perturb'
    print('perturbed one entry:', perturbed)

recs = collections.defaultdict(list)
for k in range(NC):
    for q, e in enumerate(hr[k]):
        recs[e].append((k, q))
if max(len(v) for v in recs.values()) > 15:
    refuse('a host bin with more than 15 records (4-bit counts in the checker)')


def pychk(kap):
    """Python mirror of chkPieceCount .. kap 1: the record count, or None with the reason."""
    lo, hi = win(kap)
    rk = RK[kap % 9]
    cnt = 0
    for i in range(224 * kap, 224 * kap + 224):
        load = 0
        kmin = 0
        for (k, q) in sorted(recs.get(i, [])):
            if not (kmin <= k < NC and q < R and hr[k][q] == i and lo <= time[k] <= hi):
                return None, f'bin {i}: record (call {k}, slot {q}) outside the window [{lo}, {hi}] (time {time[k]})'
            kmin = k + 1
            load += res[disp[k]][q]
        if load + rem[i] != rk:
            return None, f'bin {i}: load {load} + remainder {rem[i]} != rank {rk}'
        cnt += len(recs.get(i, []))
    return cnt, None


counts, bad = [], []
for kap in range(NK):
    c, why = pychk(kap)
    if c is None:
        bad.append((kap, why))
        c = len([1 for i in range(224 * kap, 224 * kap + 224) for _ in recs.get(i, [])])
    counts.append(c)
if bad and not (perturb and emit_neg):
    refuse(f'{len(bad)} host kind(s) fail the Python mirror; first: kind {bad[0][0]}: {bad[0][1]}')
if not bad and sum(counts) != NC * R:
    refuse('record counts do not sum to NC * R')
if perturb and not bad:
    raise SystemExit('perturbation not detected by the mirror (bug)')


def pack(fields, w):
    v = 0
    for j, f in enumerate(fields):
        assert 0 <= f < 2 ** w
        v |= f << (w * j)
    return v


def chunks(fields, w, per):
    return [pack(fields[a:a + per], w) for a in range(0, len(fields), per)]


def tree(lo, hi, leaf):
    if hi - lo == 1:
        return leaf(lo)
    m = (lo + hi) // 2
    return f"(.node {m} {tree(lo, m, leaf)} {tree(m, hi, leaf)})"


callleaf = lambda k: f"(.leaf {hex(pack(hr[k], 20))} {time[k]} {disp[k]})"
remleaf = lambda kap: f"(.leaf {hex(pack(rem[224 * kap: 224 * kap + 224], 6))} 0 0)"
cntleaf = lambda kap: f"(.leaf {counts[kap]} 0 0)"
os.makedirs(outdir, exist_ok=True)
NS = f"OAI.PowerSaving.XHostT.K.D{tag}"
mod = f"Work.XHostT.Data.{tag}"
HDR = "set_option maxHeartbeats 0\nset_option maxRecDepth 100000\n"
resT = sum(res[p][q] << (4 * (1024 * p + q)) for p in res for q in range(R))
srcsha = hashlib.sha256(raw).hexdigest()
with open(f"{outdir}/Calls.lean", 'w') as f:
    f.write(f"import Work.XHostT.KGlue\n/-! p-hostck mk_parts.py: data of host table {os.path.basename(src)} (sha256 {srcsha[:16]}..., "
            f"T = {T}, {NC} calls, R = {R}{', PERTURBED: ' + json.dumps(perturbed) if perturbed else ''}). -/\n")
    f.write(HDR + f"namespace {NS}\nopen OAI.PowerSaving.XHostT.K\n")
    f.write(f"def NC : Nat := {NC}\ndef R : Nat := {R}\ndef NK : Nat := {NK}\n")
    f.write(f"noncomputable def resT : Nat := {hex(resT)}\n")
    f.write(f"noncomputable def ct : CT := {tree(0, NC, callleaf)}\n")
    f.write(f"noncomputable def rt : CT := {tree(0, NK, remleaf)}\n")
    f.write(f"/-- the claimed record count of every host kind (checked by the kernel in the parts). -/\n")
    f.write(f"noncomputable def cntT : CT := {tree(0, NK, cntleaf)}\n")
    f.write(f"noncomputable def cK : Nat → Nat := cntOf cntT\n")
    f.write(f"end {NS}\n")


def kind_body(kap, claim):
    bf, rf = [], []
    for i in range(224 * kap, 224 * kap + 224):
        rr = sorted(recs.get(i, []))
        bf.append(len(rr))
        rf.extend(k | (q << 12) for (k, q) in rr)
    out = [f"noncomputable def binS{kap} : List Nat := [{', '.join(hex(x) for x in chunks(bf, 4, 512))}]",
           f"noncomputable def recS{kap} : List Nat := [{', '.join(hex(x) for x in chunks(rf, 22, 128))}]",
           f"theorem pc{kap} : chkPieceCount ct rt resT NC R {kap} 1 binS{kap} recS{kap} = {claim} := by decide +kernel"]
    return out, len(rf)


def rng_term(a, b):
    """a balanced KOK.app tree over the kinds a..b-1 (leaves KOK.one)."""
    if b - a == 1:
        return f"(KOK.one (a := {a}) (b := {b}) rfl binS{a} recS{a} pc{a})"
    m = (a + b) // 2
    return f"(KOK.app (a := {a}) (b := {m}) (c := {b}) {rng_term(a, m)} {rng_term(m, b)})"


if perturb and emit_neg:
    neg = []
    for kap, why in bad:
        for kind, claim, must in (('NegSome', f"some (cK {kap})", 'FAIL'), ('NegNone', 'none', 'PASS')):
            body, nr = kind_body(kap, claim)
            fn = f"{kind}{kap}.lean"
            with open(f"{outdir}/{fn}", 'w') as f:
                f.write(f"import {mod}.Calls\n/-! p-hostck mk_parts.py NEGATIVE control on {tag} (perturbed {json.dumps(perturbed)}): "
                        f"host kind {kap} ({why}); this module must {must} to compile. -/\n")
                f.write(HDR + f"namespace {NS}\nopen OAI.PowerSaving.XHostT.K\n" + "\n".join(body) + f"\nend {NS}\n")
            neg.append(dict(file=fn, kind=kap, must=must, why=why))
    json.dump(dict(src=os.path.basename(src), sha256=srcsha, perturbed=perturbed, negatives=neg),
              open(f"{outdir}/MANIFEST.json", 'w'), indent=1)
    print('NEGATIVE: wrote Calls +', len(neg), 'modules:', [(n['file'], n['must']) for n in neg])
    sys.exit(0)

with open(f"{outdir}/Order.lean", 'w') as f:
    f.write(f"import {mod}.Calls\n/-! p-hostck mk_parts.py: the global order of {tag} is valid (times distinct, below NC, increasing along each replica). -/\n")
    f.write(HDR + f"namespace {NS}\nopen OAI.PowerSaving.XHostT.K\n")
    f.write(f"theorem order_ok : chkOrder ct NC NC 0 0 = true := by decide +kernel\nend {NS}\n")
parts = []
for j, k0 in enumerate(range(0, NK, kpp)):
    k1 = min(k0 + kpp, NK)
    body, tot = [], 0
    for kap in range(k0, k1):
        b, nr = kind_body(kap, f"some (cK {kap})")
        assert nr == counts[kap]
        body += b
        tot += nr
    body.append(f"/-- part {j}: host kinds {k0} ≤ κ < {k1} pass the kernel check. -/")
    body.append(f"theorem rng{j} : KOK ct rt resT NC R cK {k0} {k1} :=\n  {rng_term(k0, k1)}")
    with open(f"{outdir}/P{j}.lean", 'w') as f:
        f.write(f"import {mod}.Calls\n/-! p-hostck mk_parts.py: part {j} of {tag}: host kinds {k0}..{k1 - 1} (one kernel theorem per kind), {tot} records. -/\n")
        f.write(HDR + f"namespace {NS}\nopen OAI.PowerSaving.XHostT.K\n" + "\n".join(body) + f"\nend {NS}\n")
    parts.append(dict(module=f"P{j}", kinds=[k0, k1], records=tot))


def whole_term(lo, hi):
    if hi - lo == 1:
        return f"rng{lo}"
    m = (lo + hi) // 2
    a, b, c = parts[lo]['kinds'][0], parts[m]['kinds'][0], parts[hi - 1]['kinds'][1]
    return f"(KOK.app (a := {a}) (b := {b}) (c := {c}) {whole_term(lo, m)} {whole_term(m, hi)})"


NP = len(parts)
with open(f"{outdir}/Whole.lean", 'w') as f:
    f.write("import Work.XHostT.KWhole\n" + f"import {mod}.Order\n" + "".join(f"import {mod}.P{j}\n" for j in range(NP)))
    f.write(f"/-! p-hostck mk_parts.py: the whole host check of {tag} ({os.path.basename(src)}, sha256 {srcsha[:16]}...): "
            f"{NP} parts, {NK} host kinds, {NC} calls x {R} slots = {NC * R} records. -/\n")
    f.write(HDR + f"namespace {NS}\nopen OAI.PowerSaving.XHostT.K\n")
    f.write(f"/-- every host kind of the table passes the kernel check (assembled from the {NP} parts). -/\n")
    f.write(f"theorem all_kinds : KOK ct rt resT NC R cK 0 NK :=\n  {whole_term(0, NP)}\n")
    f.write(f"/-- global fact: the record counts sum to `NC * R`. -/\n")
    f.write(f"theorem sum_ok : ((List.range NK).map cK).sum = NC * R := by decide +kernel\n")
    f.write(f"/-- **the whole table**: a valid cross-replica host schedule (TableOK) with a valid global order. -/\n")
    f.write(f"theorem whole : TableOK ct rt resT NC R NK ∧ chkOrder ct NC NC 0 0 = true :=\n"
            f"  ⟨tableOK_of_KOK cK all_kinds sum_ok, order_ok⟩\n")
    f.write(f"end {NS}\n#print axioms {NS}.whole\n")
# global parts: the remainder histogram and the residual sums (Work/XHostT/KHistC.lean, KHist.lean), in chunks of GK
# host kinds (core Lean modules GlobC<j>, decided in parallel) and one assembly module Glob (Mathlib)
rc = collections.Counter(rem)
V = sorted(rc)
rs = {p: sum(res[p]) for p in sorted(res)}
GK = max(dv for dv in range(1, 121) if NK % dv == 0)
NG = NK // GK
gcnt = {v: [sum(1 for i in range(224 * GK * j, 224 * GK * (j + 1)) if rem[i] == v) for j in range(NG)] for v in V}
assert all(sum(gcnt[v]) == rc[v] for v in V)
for j in range(NG):
    with open(f"{outdir}/GlobC{j}.lean", 'w') as f:
        f.write(f"import Work.XHostT.KHistC\nimport {mod}.Calls\n/-! p-hostck mk_parts.py: global chunk {j} of {tag}: host kinds {GK * j}..{GK * j + GK - 1}. -/\n")
        f.write(HDR + f"namespace {NS}\nopen OAI.PowerSaving.XHostT.K\n")
        f.write(f"theorem gin{j} : inR rt {V} {GK * j} {GK} = true := by decide +kernel\n")
        for v in V:
            f.write(f"theorem gc{v}_{j} : cntR rt {v} {GK * j} {GK} = {gcnt[v][j]} := by decide +kernel\n")
        f.write(f"end {NS}\n")
with open(f"{outdir}/Glob.lean", 'w') as f:
    f.write("import Work.XHostT.KHist\n" + "".join(f"import {mod}.GlobC{j}\n" for j in range(NG)))
    f.write(f"/-! p-hostck mk_parts.py: global facts of {tag}: remainder histogram {dict(sorted(rc.items()))}, residual sums {rs} ({NG} chunks of {GK} kinds). -/\n")
    f.write(HDR + f"namespace {NS}\nopen OAI.PowerSaving.XHostT.K\n")
    f.write(f"/-- every remainder is one of {V}. -/\ntheorem rem_in : inR rt {V} 0 NK = true :=\n"
            f"  inR_list rt {V} {GK} {NG} 0 (fun j hj => by\n    interval_cases j\n    exacts [{', '.join(f'gin{j}' for j in range(NG))}])\n")
    for v in V:
        cs = '[' + ', '.join(str(x) for x in gcnt[v]) + ']'
        f.write(f"theorem rem_c{v} : cntR rt {v} 0 NK = {rc[v]} := by\n"
                f"  have h := cntR_list rt {v} {GK} {cs} 0 {NG} rfl (fun j hj => by\n    interval_cases j\n    exacts [{', '.join(f'gc{v}_{j}' for j in range(NG))}])\n"
                f"  rw [show ({cs} : List ℕ).sum = {rc[v]} by decide] at h\n  exact h\n")
    rhs = " + ".join(f"({rc[v]} : ℕ) • f {v}" for v in V)
    f.write(f"/-- **the remainder histogram** of {tag}: any function of the remainders summed over all host indices. -/\n")
    f.write(f"theorem rem_hist {{M : Type*}} [AddCommMonoid M] (f : ℕ → M) :\n    ∑ i ∈ Finset.range (224 * NK), f (remAt rt i) = {rhs} := by\n")
    f.write(f"  rw [rem_sum rt NK {V} (by decide) rem_in f]\n")
    f.write(f"  simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, {', '.join(f'rem_c{v}' for v in V)}, add_zero, add_assoc]\n")
    for p in sorted(res):
        f.write(f"/-- the residual sum of provider {p}. -/\ntheorem res_sum{p} : ∑ q ∈ Finset.range R, resAt resT {p} q = {rs[p]} := by\n  rw [← resS_eq]; decide +kernel\n")
    f.write(f"end {NS}\n#print axioms {NS}.rem_hist\n" + "".join(f"#print axioms {NS}.res_sum{p}\n" for p in sorted(res)))
v0 = V[1] if len(V) > 1 else V[0]
with open(f"{outdir}/NegGlob.lean", 'w') as f:
    f.write(f"import Work.XHostT.KHistC\nimport {mod}.Calls\n/-! p-hostck mk_parts.py NEGATIVE control: one chunk's remainder count one too high; this module must FAIL to compile. -/\n")
    f.write(HDR + f"namespace {NS}\nopen OAI.PowerSaving.XHostT.K\ntheorem gbad : cntR rt {v0} 0 {GK} = {gcnt[v0][0] + 1} := by decide +kernel\nend {NS}\n")
json.dump(dict(src=os.path.basename(src), sha256=srcsha, T=T, NC=NC, R=R, NH=NH, NK=NK, kinds_per_part=kpp,
               rem_hist={str(v): rc[v] for v in V}, res_sums={str(p): rs[p] for p in rs}, glob_chunks=NG, glob_kinds_per_chunk=GK,
               parts=parts, total_records=sum(p['records'] for p in parts), max_records_per_bin=max(len(v) for v in recs.values()),
               providers=sorted(res), disp_counts=dict(collections.Counter(disp))),
          open(f"{outdir}/MANIFEST.json", 'w'), indent=1)
print('wrote', NP, 'parts,', NK, 'kinds, records', sum(p['records'] for p in parts), '= NC*R', NC * R)
