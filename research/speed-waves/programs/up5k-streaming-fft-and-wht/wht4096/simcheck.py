#!/usr/bin/env python3
"""simcheck.py -- simulate a generated wht_sdf design with Icarus Verilog (Yosys iCE40 cell models:
SB_MAC16, SB_SPRAM256KA, SB_LUT4, SB_CARRY) and compare every complete output frame with a numpy model.
Usage: simcheck.py design.v LOGN --samples 100000 [--seed 1] [--work DIR]
Frames: corner frames first (all max, all min, alternating max/min, impulse at 0 and at N-1, all -1,
all 1, a ramp), then uniform random 16-bit frames until --samples input samples are reached.
Model: per stage v >>= 1 (floor), pairs (a, b) at distance D -> (a - b, a + b) (swapped butterfly);
the result equals the reversed truncated WHT, which is also checked."""
import argparse, os, subprocess, sys
import numpy as np

CELLS = '$OSS_CAD_SUITE/share/yosys/ice40/cells_sim.v'
BIN = '$OSS_CAD_SUITE/bin'

def model(fr, logn):
    v = fr.astype(np.int64).copy()
    N = 1 << logn
    for i in range(logn):
        D = N >> (i + 1)
        v = v >> 1
        v = v.reshape(-1, 2, D)
        a, b = v[:, 0, :].copy(), v[:, 1, :].copy()
        v = np.stack([a - b, a + b], axis=1).reshape(N)
        assert v.min() >= -32768 and v.max() <= 32767
    return v

def model_std_reversed(fr, logn):
    v = fr.astype(np.int64).copy()
    N = 1 << logn
    for i in range(logn):
        D = N >> (i + 1)
        v = v >> 1
        v = v.reshape(-1, 2, D)
        a, b = v[:, 0, :].copy(), v[:, 1, :].copy()
        v = np.stack([a + b, a - b], axis=1).reshape(N)
    return v[::-1]

def frames(logn, total, seed):
    N = 1 << logn
    rng = np.random.default_rng(seed)
    F = []
    F.append(np.full(N, 32767)); F.append(np.full(N, -32768))
    F.append(np.where(np.arange(N) % 2 == 0, 32767, -32768))
    e = np.zeros(N, dtype=np.int64); e[0] = 32767; F.append(e)
    e = np.zeros(N, dtype=np.int64); e[-1] = -32768; F.append(e)
    F.append(np.full(N, -1)); F.append(np.full(N, 1))
    F.append((np.arange(N) * 7919 % 65536) - 32768)
    F.append(rng.choice([-32768, 32767, -1, 0, 1], size=N))
    while len(F) * N < total:
        F.append(rng.integers(-32768, 32768, size=N))
    return [f.astype(np.int64) for f in F]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('design'); ap.add_argument('logn', type=int)
    ap.add_argument('--samples', type=int, default=100000)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--work', default=None)
    ap.add_argument('--std', action='store_true', help='design uses the standard (unswapped) butterfly: natural output order')
    ap.add_argument('--netlist', action='store_true', help='design is a post-synthesis netlist (cells_sim already needed)')
    A = ap.parse_args()
    logn = A.logn; N = 1 << logn
    work = A.work or os.path.splitext(A.design)[0] + '_sim'
    os.makedirs(work, exist_ok=True)
    F = frames(logn, A.samples, A.seed)
    stream = np.concatenate(F)
    with open(os.path.join(work, 'in.hex'), 'w') as f:
        for s in stream:
            f.write('%04x\n' % (int(s) & 0xffff))
    ncyc = len(stream) + N + logn + 8
    tb = f'''`timescale 1ns/1ps
module tb;
  reg clk = 0, rst = 1; reg [15:0] x = 0; wire [15:0] y; wire y_first;
  reg [15:0] mem [0:{len(stream)-1}];
  integer k, fo;
  wht_sdf dut (.clk(clk), .rst(rst), .x(x), .y(y), .y_first(y_first));
  always #5 clk = ~clk;
  initial begin
    $readmemh("{os.path.abspath(os.path.join(work, 'in.hex'))}", mem);
    fo = $fopen("{os.path.abspath(os.path.join(work, 'out.txt'))}", "w");
    @(posedge clk); @(posedge clk); #1;
    @(negedge clk); rst = 0; x = mem[0];
    for (k = 1; k < {ncyc}; k = k + 1) begin
      @(negedge clk);
      $fwrite(fo, "%h %b\\n", y, y_first);
      x = (k < {len(stream)}) ? mem[k] : 16'h0;
    end
    $fclose(fo); $finish;
  end
endmodule
'''
    # NB: rst is 1 at the posedge before the first negedge where x = mem[0] is applied -> n = 0 while mem[0] is on x
    tbp = os.path.join(work, 'tb.v')
    open(tbp, 'w').write(tb)
    exe = os.path.join(work, 'sim.vvp')
    cmd = [BIN + '/iverilog', '-g2012', '-DICE40_U', '-o', exe, '-s', 'tb', tbp, A.design, CELLS]
    r = subprocess.run(['/usr/sbin/taskpolicy', '-b'] + cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout, r.stderr); sys.exit(2)
    r = subprocess.run(['/usr/sbin/taskpolicy', '-b', BIN + '/vvp', '-n', exe], capture_output=True, text=True)
    lines = open(os.path.join(work, 'out.txt')).read().split('\n')
    ys, fs = [], []
    for ln in lines:
        if not ln.strip():
            continue
        a, b = ln.split()
        ys.append(None if ('x' in a or 'z' in a) else int(a, 16)); fs.append(b == '1')
    # The first output line corresponds to the cycle after x = mem[0] was applied, i.e. n = 1 there?
    # Find y_first markers; frame f output index 0 is at the f-th marker that has a full frame after it.
    marks = [i for i, f in enumerate(fs) if f]
    nfr = len(F)
    ok = 0; bad = 0; first_bad = None; checked = 0
    # match: the marker whose following N outputs equal model(F[0]) defines frame 0
    exp = [model(fr, logn) for fr in F]
    if A.std:
        exp = [e[::-1].copy() for e in exp]
    for fr in F[:3]:
        assert np.array_equal(model(fr, logn), model_std_reversed(fr, logn))
    # alignment: use the theory (marker m, m+N, ...) and report which marker matched frame 0
    base = None
    for m in marks[:3]:
        seg = ys[m:m + N]
        if len(seg) == N and None not in seg:
            got = np.array([(v ^ 0x8000) - 0x8000 for v in seg])
            if np.array_equal(got, exp[0]):
                base = m; break
    if base is None and len(marks) > 1:
        print('no marker aligns frame 0; falling back to the theoretical alignment (second marker)')
        base = marks[1]
    if base is None:
        print('RESULT FAIL: no marker aligns frame 0; markers', marks[:4])
        # still report the first marker comparison
        if marks:
            seg = ys[marks[0]:marks[0] + N]
            got = [None if v is None else (v ^ 0x8000) - 0x8000 for v in seg]
            print('first outputs', got[:8], 'expected', exp[0][:8].tolist())
        sys.exit(1)
    for f in range(nfr):
        seg = ys[base + f * N: base + (f + 1) * N]
        if len(seg) < N:
            break
        if not fs[base + f * N]:
            bad += 1; first_bad = first_bad or (f, 'marker'); continue
        got = np.array([(-99999 if v is None else (v ^ 0x8000) - 0x8000) for v in seg])
        checked += 1
        if np.array_equal(got, exp[f]):
            ok += 1
        else:
            bad += 1
            if first_bad is None:
                j = int(np.nonzero(got != exp[f])[0][0])
                first_bad = (f, j, int(got[j]), int(exp[f][j]))
    # accuracy vs the exact WHT scaled by 1/N (reversed order), in output LSBs, over the random frames
    errs = []
    for f in range(9, min(nfr, 30)):
        h = F[f].astype(np.float64).copy(); hl = N // 2
        while hl >= 1:
            h = h.reshape(-1, 2, hl); h = np.stack([h[:, 0] + h[:, 1], h[:, 0] - h[:, 1]], 1).reshape(N); hl //= 2
        errs.append(np.abs(exp[f] - h[::-1] / N).max())
    if errs:
        print('model vs exact WHT/N (random frames): max abs error %.3f LSB' % max(errs))
    print('frames %d, samples %d, checked %d, ok %d, bad %d, first_bad %s, base line %d (marker index %d)'
          % (nfr, len(stream), checked, ok, bad, first_bad, base, marks.index(base)))
    print('RESULT', 'PASS' if bad == 0 and ok == nfr else 'FAIL')

if __name__ == '__main__':
    main()
