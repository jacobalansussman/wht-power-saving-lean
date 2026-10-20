#!/usr/bin/env python3
"""replay.py (pk2-orbits): rebuild design number TRIAL of an orb.py run (same mode, m, depth, seed, kinds; the runs
are deterministic) and decide it with a larger closure cap; optionally pickle it.
Usage: replay.py mode m depth seed trial cap [kinds] [save]"""
import sys, os, random, pickle, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "heights")); sys.path.insert(0, os.path.join(HERE, "..", "frame-game")); sys.path.insert(0, HERE)
from fpa import *
import orb

def build(mode, m, depth, seed, trial, kinds):
    W = 3 if mode == "three" else 2
    S = Sp(m); F = orb.FShim(S); rng = random.Random(seed)
    for tr in range(trial + 1):
        P = orb.make_prog(S, m, W, rng, mode)
        if not P.ok: continue
        O = orb.Orbit(S, P, W, rng, kinds)
        for step in range(depth + 1):
            res, st = orb.evaluate(S, F, P, W, O.nvar, rng, do_exact=True)
            if all(x == "Y" for x in res["exact"]) or step == depth: break
            st_num = orb.subs(st, res["vals"])
            dl = orb.deliver_idx(P, W); rows = orb.junk_rows(S, st_num, dl, W)
            orph = []
            for l in range(W):
                r = end_solve(S, st_num, {(l, S.K): Fr(1)}, cap_cols=6000, cap_rows=20000)
                if r[0] == "no" and "singles" in r[1] and r[1]["singles"]:
                    dmin = r[1]["singles"][0][0]
                    orph += [row for (dd, row) in r[1]["singles"] if dd <= dmin + 1]
            target_rows = orph if orph else rows
            if not target_rows: break
            rng.shuffle(target_rows)
            done = False
            for row in target_rows[:4]:
                if O.add_orbit(row, rows + orph): done = True; break
            if not done: break
    return S, P, O, res, st

if __name__ == "__main__":
    mode, m, depth, seed, trial, cap = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
    kinds = sys.argv[7].split(",") if len(sys.argv) > 7 and sys.argv[7] != "all" else ["stab", "cycle", "twof", "perm", "rand"]
    S, P, O, res, st = build(mode, m, depth, seed, trial, kinds)
    W = 3 if mode == "three" else 2
    print(f"replayed {mode} m={m} seed={seed} trial={trial}: a={getattr(P, 'a', (P.a1, P.a2))} orbits={O.log} "
          f"verdicts psilin={res['psilin']} chlin={res['chlin']} exact={res['exact']} top={res['top']}", flush=True)
    stn = orb.subs(st, res["vals"])
    print("   states:", len(stn), "terms:", sum(len(x) for x in stn), "levels of delivering states' walkers:",
          [sorted(S.level(g) for (l, g) in stn[i]) for i in orb.deliver_idx(P, W)], flush=True)
    for l in range(W):
        t0 = time.time()
        print(f"   output {l}: layered screen (cap {cap}): {orb.layer_screen(S, stn, l, cap_cols=cap, cap_rows=4 * cap)} "
              f"({time.time() - t0:.0f}s)", flush=True)
        t0 = time.time()
        r = end_solve(S, stn, {(l, S.K): Fr(1)}, cap_cols=cap, cap_rows=4 * cap)
        info = {k: v for k, v in r[1].items() if k not in ("solution", "singles")}
        print(f"   output {l}: {r[0]} {info} ({time.time() - t0:.0f}s)", flush=True)
        if r[0] == "yes":
            ok, _ = verify_combo(S, stn, r[1]["solution"], {(l, S.K): Fr(1)}); print("      re-verified:", ok)
    if "save" in sys.argv:
        with open(os.path.join(HERE, f"design_{mode}_m{m}_s{seed}_t{trial}.pkl"), "wb") as fh:
            pickle.dump(dict(states=stn, regs=P.regs, order=P.order, inj=P.inj, vals=res["vals"], log=O.log), fh)
