#!/usr/bin/env python3
# Bit-exact simulation of one PFB lane against a numpy integer model (Icarus Verilog + rtl/dsp48e2_sim.v).
#   tb.py VARIANT W NFRAMES COEFFILE [MODULE]     VARIANT in base, packed, pertap
# Frames: random 8-bit samples on both signals with the coefficients of a random branch of COEFFILE (the
# hardware's own quantised coefficients, 14,336 taps), plus corner frames: -128 / 127 on both signals,
# the samples that drive every branch's sums to their positive and negative extremes, alternating signs,
# and synthetic coefficient sets that sit exactly on the packing bound (127P + 128N = 2^18 - 1, etc.).
# Expected: ya = sum_t a_t c_t, yb = sum_t b_t c_t (exact integers). Every frame must match at one latency.
import sys, os, subprocess, numpy as np
D = os.path.dirname(os.path.abspath(__file__))
T, N = 14, 1024
var, W, NF, cfile = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
G = int(var[7:]) if var.startswith('packedg') else T          # taps per packed cascade group
NG = (T + G - 1) // G
mod = sys.argv[5] if len(sys.argv) > 5 else f"{var}_lane_w{W}"
UNISIM = os.environ.get('UNISIM') == '1'      # 1: simulate with Xilinx's own DSP48E2 model (XilinxUnisimLibrary, Apache-2.0)
work = os.path.join(D, '..', 'sim', mod + ('_unisim' if UNISIM else '') + os.environ.get('TAG', '')); os.makedirs(work, exist_ok=True)
rng = np.random.default_rng(12345)
cc = np.load(cfile).astype(np.int64).reshape(T, N).T          # cc[k, t] = c[k + N t]
assert np.abs(cc).max() <= (1 << (W - 1)) - 1, "coefficients do not fit W bits (symmetric range)"

def bound_ok(c):
    for g0 in range(0, T, G):
        cg = c[g0:g0 + G]; P = cg[cg > 0].sum(); Nn = -cg[cg < 0].sum()
        if not (127 * P + 128 * Nn <= (1 << 18) - 1 and 128 * P + 127 * Nn <= (1 << 18)): return False
    return True

# --- frames -----------------------------------------------------------------------------------------------
C, A, B = [], [], []
def add(c, a, b): C.append(np.asarray(c, np.int64)); A.append(np.asarray(a, np.int64)); B.append(np.asarray(b, np.int64))
L1 = np.abs(cc).sum(1); kmax = int(L1.argmax())
for k in (kmax, 0, N // 2, N - 1):
    c = cc[k]
    for a, b in ((-128, -128), (127, 127), (-128, 127), (127, -128), (0, -128), (-128, 0), (-1, -1)):
        add(c, np.full(T, a), np.full(T, b))
    hi = np.where(c > 0, 127, -128); lo = np.where(c > 0, -128, 127)
    for a, b in ((hi, hi), (lo, lo), (hi, lo), (lo, hi)): add(c, a, b)
# synthetic coefficient sets exactly on the bound (packed legality) and at the W-bit extremes
pk = (1 << (W - 1)) - 1
synth = []
for _ in range(40):
    v = rng.integers(-pk, pk + 1, T)
    while True:
        Pp = v[v > 0].sum(); Nn = -v[v < 0].sum()
        if not var.startswith('packed') or bound_ok(v): break
        v = (v * 0.97).astype(np.int64)
    synth.append(v)
synth.append(np.full(T, pk)); synth.append(np.full(T, -pk))
if var.startswith('packed'):  # largest legal sets
    synth = [s for s in synth if bound_ok(s)]
    synth.append(np.array([2047] + [0] * (T - 1))); synth.append(np.array([-2047] + [0] * (T - 1)))
    synth.append(np.array([2047, -1] + [0] * (T - 2)))  # 127*2047 + 128*1 = 2^18 - 1 - 0 ... checked below
    synth = [s for s in synth if bound_ok(s) and np.abs(s).max() <= pk]
for c in synth:
    hi = np.where(c > 0, 127, -128); lo = np.where(c > 0, -128, 127)
    for a, b in ((hi, hi), (lo, lo), (hi, lo), (lo, hi), (np.full(T, -128), np.full(T, -128))): add(c, a, b)
ncorner = len(C)
ks = rng.integers(0, N, NF)
for i in range(NF):
    add(cc[ks[i]], rng.integers(-128, 128, T), rng.integers(-128, 128, T))
C = np.array(C); A = np.array(A); B = np.array(B); nfr = len(C)
if var.startswith('packed'):
    assert all(bound_ok(c) for c in C), "a frame breaks the packing bound"
YA = (A * C).sum(1); YB = (B * C).sum(1); SG = np.stack([C[:, g0:g0 + G].sum(1) for g0 in range(0, T, G)], 1)

# --- stimulus: one line per cycle, ports concatenated (first port in the low bits) --------------------------
skew = var != 'pertap'
NCYC = nfr + T + 30
def tw(v, w): return int(v) & ((1 << w) - 1)
ports = [('a', T * 8), ('b', T * 8), ('c', T * W)] if var != 'base' else [('x', T * 8), ('c', T * W)]
if var.startswith('packed'): ports.append(('negs', 16 * NG))
if var == 'pertap': ports.append(('negc', T * (W + 1)))
IW = sum(w for _, w in ports)
lines = []
for k in range(NCYC):
    val = 0; off = 0
    for nm, w in ports:
        f = 0
        if nm in ('a', 'b', 'x', 'c', 'negc'):
            ew = 8 if nm in ('a', 'b', 'x') else (W if nm == 'c' else W + 1)
            for t in range(T):
                # tap t arrives (t mod G) cycles after its group's first tap, and a short last group starts later
                # so that every group ends on the same cycle (address offsets in the real filter)
                n = k - (t % G) - (G - min(G, T - (t // G) * G)) if skew else (k - 1 if nm == 'negc' else k)
                if 0 <= n < nfr:
                    v = {'a': A, 'b': B, 'x': A, 'c': C}[nm][n, t] if nm != 'negc' else -C[n, t]
                    f |= tw(v, ew) << (t * ew)
        elif nm == 'negs':
            for g in range(NG):
                n = k - 1 - (G - min(G, T - g * G))
                if 0 <= n < nfr: f |= tw(-SG[n, g], 16) << (16 * g)
        val |= f << off; off += w
    lines.append(f"{val:0{(IW + 3) // 4}x}")
open(f"{work}/stim.hex", 'w').write('\n'.join(lines) + '\n')

outs = ['y'] if var == 'base' else ['ya', 'yb']
ow = {'base': 48, 'packed': 19, 'pertap': 23}.get(var, 19 + (NG - 1).bit_length())
conn = ', '.join(f".{nm}(st[{o + w - 1}:{o}])" for (nm, w), o in zip(ports, np.cumsum([0] + [w for _, w in ports])[:-1]))
conn += ', ' + ', '.join(f".{o}({o})" for o in outs)
tb = f"""`timescale 1ns/1ps
module tb;
  reg clk = 0; reg [{IW-1}:0] mem [0:{NCYC-1}]; reg [{IW-1}:0] st = 0; integer k, f;
  {' '.join(f'wire [{ow-1}:0] {o};' for o in outs)}
  {mod} dut (.clk(clk), {conn});
  initial begin
    $readmemh("{work}/stim.hex", mem); f = $fopen("{work}/out.txt", "w");
    {'#300;  // wait for glbl GSR (100 ns) to be released' if UNISIM else ''}
    for (k = 0; k < {NCYC}; k = k + 1) begin
      st = mem[k]; #1 clk = 1; #1 clk = 0;
      $fdisplay(f, "{' '.join(['%h'] * len(outs))}", {', '.join(outs)});
    end
    $fclose(f); $finish;
  end
endmodule
"""
open(f"{work}/tb.v", 'w').write(tb)
U = os.path.join(D, '..', 'unisim', 'src')
src = ([f"{U}/DSP48E2.v", f"{U}/glbl.v"] if UNISIM else [f"{D}/dsp48e2_sim.v"]) + [f"{D}/{mod}.v"]
subprocess.run(['iverilog', '-g2012', '-o', f"{work}/sim.vvp"] + (['-s', 'tb', '-s', 'glbl'] if UNISIM else []) + [f"{work}/tb.v"] + src, check=True)
r = subprocess.run(['vvp', '-n', f"{work}/sim.vvp"], check=True, capture_output=True, text=True)
msgs = [l for l in (r.stdout + r.stderr).splitlines() if l.strip()]
if msgs: print("simulator messages:", *msgs[:8], sep='\n  ')
raw = [l.split() for l in open(f"{work}/out.txt")]
def sx(h, w):
    if 'x' in h or 'z' in h: return None
    v = int(h, 16); return v - (1 << w) if v >> (w - 1) else v
got = [[sx(h, ow) for h in r] for r in raw]
exp = [(YA[i],) if var == 'base' else (YA[i], YB[i]) for i in range(nfr)]
best = None
for lat in range(0, 30):
    bad = sum(1 for i in range(nfr) if i + lat >= len(got) or tuple(got[i + lat]) != exp[i])
    if best is None or bad < best[1]: best = (lat, bad)
lat, bad = best
first = next((i for i in range(nfr) if tuple(got[i + lat]) != exp[i]), None)
print(f"{mod}{' [Xilinx unisim DSP48E2]' if UNISIM else ' [own DSP48E2 model]'}: {nfr} frames ({ncorner} corner, {NF} random), best latency {lat} cycles (cycle k shows frame k - lat), "
      f"mismatching frames {bad}" + (f", first at frame {first}: got {got[first + lat]} expected {exp[first]}" if bad else "")
      + ("  PASS" if bad == 0 else "  FAIL"))
