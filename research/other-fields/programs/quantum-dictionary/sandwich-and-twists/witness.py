#!/usr/bin/env python3
"""witness.py -- rerun twist.py's m=3 machinery and print one explicit twisted witness for the NON-orthogonal pair
u_S=001, u_T=111 (timed version), with every frame given by a basis of Pauli vectors (x|z) and its excess."""
import sys, runpy
sys.argv = ["twist.py", "3"]
import builtins
src = open(__file__.replace("witness.py", "twist.py")).read()
src = src[:src.index('out("\\nTwisted pair search')]
g = {"__file__": __file__.replace("witness.py", "twist.py"), "__name__": "twistlib"}
exec(compile(src, "twist.py", "exec"), g)
np = g["np"]; P = g["P"]; D = g["D"]; frames = g["frames"]; I_mask = g["I_mask"]; port = g["port"]; leqF = g["leqF"]
iF, i0, exc, desc = g["iF"], g["i0"], g["exc"], g["G"]["desc"]; m = 3
uS, uT = 0b001, 0b111
best, found, (aS, eS, aT, eT) = g["run_pair"](uS, uT, timed=True)
l, a0, ast, k = found[0]
phi = P[k]
XS = np.nonzero(I_mask(aS, iF))[0]; YS = np.nonzero(I_mask(i0, eS))[0]
XT = I_mask(aT, iF); YT = I_mask(i0, eT)
def show(f):
    b = ["(%s|%s)" % (format(v & 7, "03b")[::-1], format(v >> 3, "03b")[::-1]) for v in desc[f]]
    return "frame %3d basis %s excess %d%s" % (f, " ".join(b), exc[f], "" if exc[f] else " (subspace frame)")
print("pair u_S=001, u_T=111 (u_S.u_T = 1: S does not know T in the old sense); copy cost l =", l)
print("phi = element %d of Sp(6,2) (twist of x_T relative to x_S inside the total R)" % k)
print("R reads X_S at a0:      ", show(a0)); print("   and X_T at phi(a0):  ", show(phi[a0]))
print("R delivers x_S at a*:   ", show(ast)); print("   and x_T at phi(a*):  ", show(phi[ast]))
mus = [mu for mu in YS if XT[phi[mu]]]
nus = [nu for nu in XS if YT[phi[nu]]]
print("Y_S cancels x_T with X_T itself when Y_S is at mu and X_T at phi(mu):")
for mu in mus: print("   mu:", show(mu), "\n   phi(mu):", show(phi[mu]))
print("Y_T cancels x_S with X_S itself when X_S is at nu and Y_T at phi(nu):")
for nu in nus: print("   nu:", show(nu), "\n   phi(nu):", show(phi[nu]))
# R's walk a0 -> a*: frames of R's two terms along every geodesic
mid = [z for z in range(len(frames)) if D[a0, z] == 1 and D[z, ast] == 1]
print("geodesic midpoints of R's walk (x_S term) and the x_T term's frame there:")
for z in mid: print("   x_S:", show(z), "| x_T:", show(phi[z]))
