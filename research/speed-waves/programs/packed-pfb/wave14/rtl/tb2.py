#!/usr/bin/env python3
# as-4: bit-exact simulation of one graded-group packed PFB lane (gen2.py) against a numpy integer model.
#   tb2.py FAM SCHEME W NFRAMES COEFFILE [MODULE]      (COEFFILE: EFFECTIVE coefficients c_eff = c_raw << shift, 14,336)
#   UNISIM=1: Xilinx's own DSP48E1 / DSP48E2 model (XilinxUnisimLibrary, Apache-2.0) instead of my model.
# Frames: random 8-bit a, b with the raw coefficients of a random branch; corners: -128/127 on both signals, samples
# that drive every group's sums to their extremes on the largest-L1 branches, and synthetic coefficient sets that
# sit exactly on the per-group bound (127P + 128N = 2^(s-1) - 1, 128P + 127N = 2^(s-1)).
# Expected: ya = sum_t a_t c_eff_t, yb = sum_t b_t c_eff_t. Every frame must match at one latency.
import sys, os, subprocess, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
import gen2
T, N = 14, 1024
fam, sname, W, NF, cfile = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
mod = sys.argv[6] if len(sys.argv) > 6 else f"glane_{fam}_{sname}_w{W}"
s = 19 if fam == 'e2' else 17
sc = gen2.schemes()[sname]; groups = gen2.order(sc); ng = len(groups); ow = gen2.out_width(fam, sc, W); UG = gen2.unpacked(sc)
sh = np.zeros(T, int); gid = np.zeros(T, int); pos = np.zeros(T, int); glen = np.zeros(T, int)
for g, taps, s_ in groups:
    for j, t in enumerate(taps): sh[t] = s_; gid[t] = g; pos[t] = j; glen[t] = len(taps)
Gmax = max(len(t) for _, t, _ in groups)
UNISIM = os.environ.get('UNISIM') == '1'
work = os.path.join(D, '..', 'sim', mod + ('_unisim' if UNISIM else '') + ('_yosys' if os.environ.get('MODEL') == 'yosys' else '') + os.environ.get('TAG', '')); os.makedirs(work, exist_ok=True)
rng = np.random.default_rng(int(os.environ.get('SEED', '12345')))
ce = np.load(cfile).astype(np.int64).reshape(T, N)
raw = ce >> sh[:, None]; assert np.all((raw << sh[:, None]) == ce), "coefficient file does not match the scheme's shifts"
cc = raw.T                                       # cc[k, t] raw
assert np.abs(cc).max() <= (1 << (W - 1)) - 1, "raw coefficients do not fit W bits"
lim = 1 << (s - 1)
def bound_ok(c):
    for g in range(ng):
        if g in UG: continue
        cg = c[gid == g]; P = cg[cg > 0].sum(); Nn = -cg[cg < 0].sum()
        if not (127 * P + 128 * Nn <= lim - 1 and 128 * P + 127 * Nn <= lim): return False
    return True
assert all(bound_ok(cc[k]) for k in range(N)), "coefficient file breaks the packing bound"
C, A, B = [], [], []
def add(c, a, b): C.append(np.asarray(c, np.int64)); A.append(np.asarray(a, np.int64)); B.append(np.asarray(b, np.int64))
L1g = np.stack([np.abs(cc[:, gid == g]).sum(1) for g in range(ng)], 1)
ks = sorted(set([int(x) for x in L1g.argmax(0)] + [0, N // 2, N - 1, int(np.abs(cc).sum(1).argmax())]))
for k in ks:
    c = cc[k]
    for a, b in ((-128, -128), (127, 127), (-128, 127), (127, -128), (0, -128), (-128, 0), (-1, -1)):
        add(c, np.full(T, a), np.full(T, b))
    hi = np.where(c > 0, 127, -128); lo = np.where(c > 0, -128, 127)
    for a, b in ((hi, hi), (lo, lo), (hi, lo), (lo, hi)): add(c, a, b)
pk = (1 << (W - 1)) - 1
synth = []
for _ in range(60):          # random sets pushed onto the bound of every group
    v = rng.integers(-pk, pk + 1, T)
    for g in range(ng):
        m = gid == g
        if g in UG: continue
        while True:
            P = v[m & (v > 0)].sum(); Nn = -v[m & (v < 0)].sum()
            if 127 * P + 128 * Nn <= lim - 1 and 128 * P + 127 * Nn <= lim: break
            v[m] = (v[m] * 0.97).astype(np.int64)
        # top up one positive coefficient so that 127P + 128N is as close to the bound as W allows
        P = v[m & (v > 0)].sum(); Nn = -v[m & (v < 0)].sum(); room = (lim - 1 - 127 * P - 128 * Nn) // 127
        idx = np.nonzero(m)[0][0]; v[idx] = min(pk, v[idx] + max(0, room)) if v[idx] >= 0 else v[idx]
    if bound_ok(v): synth.append(v.copy())
for g in [g_ for g_ in range(ng) if g_ not in UG]:          # one group exactly on 128P + 127N = 2^(s-1) territory with a single big coefficient
    v = np.zeros(T, np.int64); t0 = np.nonzero(gid == g)[0][0]; v[t0] = min(pk, (lim - 1) // 127); synth.append(v.copy()); synth.append(-v)
synth = [v for v in synth if bound_ok(v) and np.abs(v).max() <= pk]
for c in synth:
    hi = np.where(c > 0, 127, -128); lo = np.where(c > 0, -128, 127)
    for a, b in ((hi, hi), (lo, lo), (hi, lo), (lo, hi), (np.full(T, -128), np.full(T, -128))): add(c, a, b)
ncorner = len(C)
kk = rng.integers(0, N, NF)
for i in range(NF): add(cc[kk[i]], rng.integers(-128, 128, T), rng.integers(-128, 128, T))
C = np.array(C); A = np.array(A); B = np.array(B); nfr = len(C)
CE = C << sh[None, :]
YA = (A * CE).sum(1); YB = (B * CE).sum(1)
SG = np.stack([C[:, gid == g].sum(1) for g in range(ng)], 1)
NCYC = nfr + Gmax + 30
def tw(v, w): return int(v) & ((1 << w) - 1)
ports = [('a', T * 8), ('b', T * 8), ('c', T * W), ('negs', 16 * ng)]
IW = sum(w for _, w in ports)
off_t = pos + (Gmax - glen)                                   # arrival offset of tap t
off_g = [Gmax - len(taps) + 1 for _, taps, _ in groups]       # C value of group g: one cycle after its first DSP
lines = []
for k in range(NCYC):
    val = 0; o = 0
    for nm, w in ports:
        f = 0
        if nm == 'negs':
            for g in range(ng):
                n = k - off_g[g]
                if 0 <= n < nfr: f |= tw(-SG[n, g], 16) << (16 * g)
        else:
            ew = 8 if nm in ('a', 'b') else W
            for t in range(T):
                n = k - off_t[t]
                if 0 <= n < nfr: f |= tw({'a': A, 'b': B, 'c': C}[nm][n, t], ew) << (t * ew)
        val |= f << o; o += w
    lines.append(f"{val:0{(IW + 3) // 4}x}")
open(f"{work}/stim.hex", 'w').write('\n'.join(lines) + '\n')
conn = ', '.join(f".{nm}(st[{o_ + w - 1}:{o_}])" for (nm, w), o_ in zip(ports, np.cumsum([0] + [w for _, w in ports])[:-1]))
tb = f"""`timescale 1ns/1ps
module tb;
  reg clk = 0; reg [{IW-1}:0] mem [0:{NCYC-1}]; reg [{IW-1}:0] st = 0; integer k, f;
  wire [{ow-1}:0] ya, yb;
  {mod} dut (.clk(clk), {conn}, .ya(ya), .yb(yb));
  initial begin
    $readmemh("{work}/stim.hex", mem); f = $fopen("{work}/out.txt", "w");
    {'#300;' if UNISIM else ''}
    for (k = 0; k < {NCYC}; k = k + 1) begin
      st = mem[k]; #1 clk = 1; #1 clk = 0;
      $fdisplay(f, "%h %h", ya, yb);
    end
    $fclose(f); $finish;
  end
endmodule
"""
open(f"{work}/tb.v", 'w').write(tb)
U = os.path.join(D, '..', 'unisim', 'src'); prim = 'DSP48E1' if fam == 'e1' else 'DSP48E2'
YOSYSM = os.environ.get('MODEL') == 'yosys'   # third model: Yosys' own behavioural DSP48E1 (cells_sim.v)
src = ([f"{U}/{prim}.v", f"{U}/glbl.v"] if UNISIM else [f"{D}/yosys_dsp48e1.v"] if YOSYSM else [f"{D}/{'dsp48e1_sim.v' if fam == 'e1' else 'dsp48e2_sim.v'}"]) + [f"{D}/{mod}.v"]
subprocess.run(['iverilog', '-g2012', '-o', f"{work}/sim.vvp"] + (['-s', 'tb', '-s', 'glbl'] if UNISIM else []) + [f"{work}/tb.v"] + src, check=True)
r = subprocess.run(['vvp', '-n', f"{work}/sim.vvp"], check=True, capture_output=True, text=True)
msgs = [l for l in (r.stdout + r.stderr).splitlines() if l.strip()]
if msgs: print("simulator messages:", *msgs[:8], sep='\n  ')
def sx(h_, w):
    if 'x' in h_ or 'z' in h_: return None
    v = int(h_, 16); return v - (1 << w) if v >> (w - 1) else v
got = [tuple(sx(x, ow) for x in l.split()) for l in open(f"{work}/out.txt")]
exp = [(YA[i], YB[i]) for i in range(nfr)]
best = None
for lat in range(0, 40):
    bad = sum(1 for i in range(nfr) if i + lat >= len(got) or got[i + lat] != exp[i])
    if best is None or bad < best[1]: best = (lat, bad)
lat, bad = best
first = next((i for i in range(nfr) if got[i + lat] != exp[i]), None)
print(f"{mod} [{'Xilinx unisim ' + prim if UNISIM else 'Yosys cells_sim ' + prim if YOSYSM else 'own ' + prim + ' model'}]: {nfr} frames ({ncorner} corner, {NF} random), "
      f"latency {lat}, mismatching frames {bad}" + (f", first at frame {first}: got {got[first + lat]} expected {exp[first]}" if bad else "")
      + ("  PASS" if bad == 0 else "  FAIL"))
