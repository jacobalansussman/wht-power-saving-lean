# sim.py (as-5): bit-exact check of a generated FFT (gen.py) against the numpy/python fixed-point model.
#   python sim.py rtl_sc_256.v [--frames 400] [--seed 1]
# Frames: 6 corner frames (impulse +32767, impulse -32768, full-scale complex tone, all -32768, all +32767,
# full-scale real cosine) then random full-range frames, back to back (no gaps). Prints mismatches and SQNR.
import sys, os, json, subprocess, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen import model

ap = argparse.ArgumentParser()
ap.add_argument('rtl')
ap.add_argument('--frames', type=int, default=400)
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--work', default='simwork')
ap.add_argument('--lib', action='append', default=[])
ap.add_argument('--incdir', action='append', default=[])
a = ap.parse_args()
meta = json.load(open(a.rtl.replace('_broken', '') + '.json')) if not os.path.exists(a.rtl + '.json') else json.load(open(a.rtl + '.json'))
n = meta['n']; N = 1 << n; lat = meta['lat']; P = meta['P']; name = meta['name']
os.makedirs(a.work, exist_ok=True)
rng = np.random.default_rng(a.seed)
F = a.frames
xr = np.zeros((F, N), dtype=np.int64); xi = np.zeros((F, N), dtype=np.int64)
t = np.arange(N)
xr[0, 0] = 32767
xr[1, 5] = -32768; xi[1, 5] = -32768
xr[2] = np.clip(np.round(32767 * np.cos(2 * np.pi * 3 * t / N)), -32768, 32767)
xi[2] = np.clip(np.round(32767 * np.sin(2 * np.pi * 3 * t / N)), -32768, 32767)
xr[3] = -32768; xi[3] = -32768
xr[4] = 32767; xi[4] = 32767
xr[5] = np.clip(np.round(32767 * np.cos(2 * np.pi * 17.3 * t / N)), -32768, 32767)
xr[6:] = rng.integers(-32768, 32768, size=(F - 6, N)); xi[6:] = rng.integers(-32768, 32768, size=(F - 6, N))
# stream: F frames + 1 zero frame + latency of zeros
tot = (F + 1) * N + lat + 4
sr = np.zeros(tot, dtype=np.int64); si = np.zeros(tot, dtype=np.int64)
sr[:F * N] = xr.reshape(-1); si[:F * N] = xi.reshape(-1)
with open(os.path.join(a.work, name + '_in.hex'), 'w') as f:
    for k in range(tot):
        f.write('%04x%04x\n' % (int(sr[k]) & 0xFFFF, int(si[k]) & 0xFFFF))
tb = '''`timescale 1ns/1ps
module tb;
  reg clk = 0, rst = 1; reg [31:0] mem [0:%d]; reg signed [15:0] in_re, in_im; wire signed [15:0] out_re, out_im;
  integer k, fo;
  %s dut(.clk(clk), .rst(rst), .in_re(in_re), .in_im(in_im), .out_re(out_re), .out_im(out_im));
  always #5 clk = ~clk;
  initial begin
    $readmemh("%s", mem); fo = $fopen("%s", "w"); in_re = 0; in_im = 0;
    repeat (3) @(posedge clk); #1 rst = 0;
    for (k = 0; k < %d; k = k + 1) begin in_re = mem[k][31:16]; in_im = mem[k][15:0]; @(posedge clk); #1 $fwrite(fo, "%%0d %%0d\\n", out_re, out_im); end
    $fclose(fo); $finish;
  end
endmodule
''' % (tot - 1, name, os.path.join(a.work, name + '_in.hex'), os.path.join(a.work, name + '_out.txt'), tot)
# NOTE on timing: rst is sampled high at the edge before k=0, so gcnt = 0 while sample 0 is applied; the value
# written after the edge that ends cycle k is the output register's value during cycle k+1.
tbf = os.path.join(a.work, name + '_tb.v')
open(tbf, 'w').write(tb)
exe = os.path.join(a.work, name + '.vvp')
subprocess.check_call(['iverilog', '-g2012', '-o', exe] + ['-I' + d for d in a.incdir] + [tbf, a.rtl] + a.lib)
subprocess.check_call(['vvp', '-n', exe], stdout=subprocess.DEVNULL)
def _p(v):
    try:
        return int(v)
    except ValueError:
        return 999999  # x/z (pipeline not yet filled)
out = np.array([[_p(v) for v in ln.split()] for ln in open(os.path.join(a.work, name + '_out.txt'))], dtype=np.int64)
# slot -> position map
posmap = np.zeros(N, dtype=np.int64)
for q in range(N):
    p = 0
    for tb_ in range(n):
        if (q >> tb_) & 1:
            p |= 1 << P[tb_]
    posmap[q] = p
bad = 0; checked = 0; first = None
err = []; ref = []
for f in range(F):
    mr, mi = model(xr[f], xi[f], n)
    mr = np.array(mr); mi = np.array(mi)
    # output during cycle t (t from 0 = sample 0's cycle) is out[t-1]; slot of cycle t is t - lat
    base = f * N + lat - 1
    hr = out[base:base + N, 0][np.argsort(posmap)]  # by position
    hi = out[base:base + N, 1][np.argsort(posmap)]
    m = (hr != mr) | (hi != mi)
    checked += N
    if m.any():
        bad += int(m.sum())
        if first is None:
            first = (f, int(np.argmax(m)))
    # precision vs float FFT (positions are bit-reversed bins), scaled by 1/N
    if f >= 6:
        X = np.fft.fft(xr[f] + 1j * xi[f]) / N
        br = np.array([int(format(i, '0%db' % n)[::-1], 2) for i in range(N)])
        err.append((mr + 1j * mi) - X[br]); ref.append(X[br])
e = np.concatenate(err); r = np.concatenate(ref)
sqnr = 10 * np.log10(np.mean(np.abs(r) ** 2) / np.mean(np.abs(e) ** 2))
print(json.dumps({'design': name, 'samples_checked': checked, 'mismatches': bad, 'first_mismatch_frame_pos': first,
                  'result': 'PASS' if bad == 0 else 'FAIL', 'sqnr_db_model_vs_float_random_frames': round(float(sqnr), 3)}))
