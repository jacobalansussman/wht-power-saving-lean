# v-e3-boost vcmp.py: certlib's A, C_k and distance-set profiles vs mine (all M <= 8, every S; intervals up to 16);
# projection closed forms vs the exact counts in e3-boost/proj/N*.txt.
import sys, re
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), ".."))
import certlib
import vsurv as S
from fractions import Fraction
mine = {}
n = bad = 0
for name, M, P, st in certlib.PROFILES:
    if name.startswith("Dist"):
        tag = name[len(f"Dist{M}S"):]
        Sx = tuple(map(int, tag)) if M <= 9 else tuple(range(int(tag.split("-")[0]), int(tag.split("-")[1]) + 1))
        Q = S.Dprof(M, Sx)
    elif re.fullmatch(r"C\d+k\d+", name):
        k = int(name.split("k")[1]); Q = S.Cprof(M, k)
    elif re.fullmatch(r"A\d+", name):
        Q = S.Aprof(M)
    else: continue
    n += 1
    if [Fraction(x) for x in P] != [Fraction(x) for x in Q]: bad += 1; print("MISMATCH", name, P[:4], Q[:4])
print(f"{n} certlib profiles (A, C_k, Dist) compared with mine: mismatches {bad}")
import os, glob
for f in sorted(glob.glob(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "..", "proj", "N*.txt"))):
    txt = open(f).read()
    print(os.path.basename(f), txt[:300].replace("\n", " | "))
    M = int(re.search(r"N(\d+)", os.path.basename(f)).group(1))
    N = certlib.Nproj(M)
    print("   closed form N_r:", {r: N[r] for r in range(1, M + 1)})
