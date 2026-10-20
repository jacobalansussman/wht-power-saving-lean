# count.py (as-5): step 1 of the brief, by hand [derived]. Butterfly-stage resources of a one-sample-per-clock
# radix-2^2 DIF FFT, natural-order input, W = 16-bit components, counted in W-bit real operations and in words.
# Cost model (ECP5/NX with Yosys): an adder or adder/subtractor bit = 1 LUT4 position (the add/sub select rides on
# the carry cell's C input); a 2:1 mux bit = 1 LUT4; a 3:1 mux bit ~ 2 LUT4; flip-flop bits counted apart.
import sys
W = 16
def sdf(n):
    st = []
    for s in range(1, n + 1):
        even = s % 2 == 0
        add = 4; mux = 4 + (2 + 2 if even else 0)   # out+memin muxes; BF2II: b swap (2) + imag S/D role swap (2)
        st.append(dict(add=add, mux=mux, ff=2, mem=2 ** (n - s)))
    return st
def sdfx(n):  # muxes folded into the carry cells (ECP5 CCU2C), only the BF2II b swap stays
    st = []
    for s in range(1, n + 1):
        even = s % 2 == 0
        st.append(dict(add=4, mux=2 if even else 0, ff=2, mem=2 ** (n - s)))
    return st
def sc(n):    # serial-butterfly SC as built here: 1 complex add/sub, bit exchange (2 complex muxes), 3 complex regs
    P = list(range(n)); st = []
    for s in range(1, n + 1):
        k = P.index(n - s)
        mem = (2 ** k - 1) if k > 0 else 0
        if k > 0:
            P[0], P[k] = P[k], P[0]
        even = s % 2 == 0
        st.append(dict(add=2, mux=(4 if k > 0 else 0) + (2 if even else 0) + 1, ff=8, mem=mem))  # +1: carry-in select
    return st
for n in (6, 8, 10):
    print('N = %d' % (1 << n))
    for nm, f in (('R2^2SDF (generic)', sdf), ('R2^2SDF, muxes in carry cells', sdfx), ('serial-butterfly SC', sc)):
        st = f(n)
        add = sum(x['add'] for x in st) * W; mux = sum(x['mux'] for x in st) * W
        ff = sum(x['ff'] for x in st) * W; mem = sum(x['mem'] for x in st)
        print('  %-32s adder bits %5d  mux bits %5d  LUT4 (adder+mux) %5d  FF bits %5d  memory words %5d' % (nm, add, mux, add + mux, ff, mem))
