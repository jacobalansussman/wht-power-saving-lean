"""xw3-extract: one row per FOLLOWUP.md flag, read from the RESULT lines of the fu.py runs (out_*.txt here).

Rule (mechanical, from the fu.py stages):
 - f=3,4 collapse space X (nt vectors, = path_ext's NONTRIVIAL count when the seeds agree);
 - f=5 (fresh inputs, same coefficients, complete span): SURVIVORS = vectors of X whose residual is < 1e-6 (the
   subspace spectrum, not the per-vector residuals, decides);
 - survivors are class (b) when NO collapse of the irreducible columns survives against lower + reducible columns
   (Bplus survivors 0 and Bplus N_count 0): every survivor is then, at f = 5, a combination of reducible columns
   (words with a direction applied twice, D_z D_z = X^z) that lands in the menu span through the P identity;
 - the non-surviving directions are class (a); any Bplus survivor would be a class (c) candidate.
C5 ("1;2;4,A.T2.Ai") is not a separate run: C = A.T2.Ai commutes with every X4_j (A maps row-2 translations to
row-2 translations, T2 commutes with them), so C D4 C^-1 = D4 at every f [proved; out_frame_check.txt: 0.0 at
f = 3..6], and batch C's numbers for C5 are identical to C2's.
"""
import json, glob, os, sys
here = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))  # folder with the fu.py outputs out_*.txt
FLAGS = [("B1", "batch B: 1;2;4;7, W=(1, K5, H15), B=H15, S=1", 1),
         ("B2", "batch B: 1;2;4;4, W=(1, K5, H15), B=H15, S=1", 5),
         ("B3", "batch B: 1;2;4;7, W=(1, K5, H30), B=H30, S=2", 1),
         ("C1", "batch C: 1,B;2;4,A, W=(K5, H30), B=H30, S=1", 2),
         ("C2", "batch C: 1;2;4, W=(K9, H30), B=H30, S=2", 8),
         ("C3", "batch C: 1,B;2;4,A, W=(K9, H30), B=H30, S=2", 10),
         ("C4s601", "batch C: 1;2;7, W=(K9, H30), B=H30, S=2 (f=5 seed 601)", 7)]
rows = {}
for p in sorted(glob.glob(os.path.join(here, "out_*.txt"))):
    for line in open(p):
        if line.startswith("RESULT "):
            r = json.loads(line[7:]); rows[r["tag"]] = r
fmt = lambda xs: "[" + ", ".join(f"{x:.2g}" for x in xs) + "]" if xs else "-"
print("flag | configuration | flagged | f=3,4 extracted | f=5 subspace spectrum (relative distance, ascending) | "
      "survivors f=5 | f=6 residual of survivors | irreducible vs lower+reducible: f=3,4 collapses -> f=5 survivors | class")
out = {}
for tag, conf, flagged in FLAGS:
    r = rows[tag]; A, B, Bp, Ap = r["A"], r["B"], r.get("Bplus", {}), r.get("Aplus", {})
    nt, surv = A["nt"], B["survivors"]
    cand_c = Bp.get("survivors", 0) > 0 or Bp.get("N_count", 0) > 0
    if cand_c: cls = f"(c) CANDIDATE: {Bp.get('survivors')}"
    elif surv == 0: cls = f"(a) all {nt}"
    else: cls = f"(b) {surv}, (a) {nt - surv}"
    out[tag] = dict(flagged=flagged, nt=nt, survivors=surv, cls=cls)
    print(f"{tag} | {conf} | {flagged} | {nt} | {fmt(B['spectrum'])} | {surv} | {fmt(r.get('F6', []))} | "
          f"{Ap.get('nt')} -> {Bp.get('survivors', '-')} | {cls}")
c2 = out["C2"]
print(f"C5 | batch C: 1;2;4,A.T2.Ai, W=(K9, H30), B=H30, S=2 | 8 | = C2 (frame commutes with D4) | as C2 | "
      f"{c2['survivors']} | as C2 | as C2 | {c2['cls']}")
