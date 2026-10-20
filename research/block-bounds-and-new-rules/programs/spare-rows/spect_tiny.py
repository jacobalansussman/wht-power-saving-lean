# spect_tiny.py (d-G3-SPECTATOR/more-spectators, wave 17).  The equal-frame search tools of wht33/tinygame (imported
# read-only: frames.py, game.py, cutgen.py, core.c) with RANK-2 sessions: a payload starting at frame s must end at
# some frame at distance exactly 2 from s (the target B_2 (x) I_(m-2): free letters move the end frame through every
# frame at distance 2, since the stabiliser of L_0 is transitive on each distance class).  Budget 2W - 1 moves.
# Per pattern class: cut-LP bound (valid for every field), then the exact level search over GF(p) if not excluded.
# Usage: python3 -I spect_tiny.py KIND M W P MINUTES     e.g. spect_tiny.py sub 4 4 3 40
import sys, time, json, os
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# (needs numpy; the original loaded a private copy of numpy and scipy here)
import frames, game

kind, m, W, p, minutes = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
t0 = time.time()
G = frames.graph(kind, m); n, adj = G['n'], G['adj']
classes = frames.pattern_classes(G, W)
print(f"{kind} m={m}: {n} frames; W={W}: {len(classes)} pattern classes ({time.time()-t0:.1f}s)", flush=True)
dist = {}
def d2(s):
    if s not in dist: dist[s] = sorted(j for j, d in frames.bfs_dist(adj, s).items() if d == 2)
    return dist[s]
rows = []
for rep, size in classes:
    starts = list(rep); tsets = [d2(s) for s in starts]
    lp = game.cut_lp(n, adj, starts, tsets)
    rows.append((lp if lp is not None else -1, starts, tsets, size))
rows.sort(key=lambda r: r[0])
budget = 2 * W - 1
summary = {}
out = open(f"res_tiny_{kind}{m}_W{W}_p{p}.jsonl", "w")
for lp, starts, tsets, size in rows:
    rec = {"graph": kind, "m": m, "W": W, "class": starts, "orbit": size, "cut_lp": lp, "budget": budget, "field": p}
    if lp > budget + 1e-6:
        rec["status"] = "excluded by cut"
    elif time.time() - t0 > minutes * 60:
        rec["status"] = "not run (time box)"
    else:
        res = game.run_core(n, adj, starts, tsets, p, budget, mem_mb=800, time_s=240, cuts="lp")
        rec["status"] = {"none": "searched: no scheme", "found": "SCHEME FOUND (GF(p))"}.get(res.get("status"), res.get("status"))
        rec["states"] = res.get("states"); rec["time"] = res.get("time")
        if res.get("status") == "found": rec["moves"] = res.get("moves")
    summary[rec["status"]] = summary.get(rec["status"], 0) + 1
    out.write(json.dumps(rec) + "\n"); out.flush()
    if rec["status"] not in ("excluded by cut",):
        print(f"  class {starts} (dims {[len(G['frames'][s]).bit_length()-1 for s in starts]}): cut LP {lp:.4f} -> {rec['status']}"
              f" states={rec.get('states')} ({time.time()-t0:.0f}s)", flush=True)
lps = [r[0] for r in rows]
print(f"SUMMARY {kind} m={m} W={W} GF({p}) budget {budget}: {summary}; cut-LP min {min(lps):.4f} max {max(lps):.4f} ({time.time()-t0:.0f}s)")
