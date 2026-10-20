# sqnr64.py (as-5): SQNR of the published r22sdf FFT64 (MIT) and of our 64-point model on the same inputs
# (complex uniform noise inside the circle |z| <= 0.7 * 32767, which r22sdf needs: it has no saturation).
import sys, os, subprocess, numpy as np
sys.path.insert(0, '.')
from gen import model
N = 64; F = 400; rng = np.random.default_rng(7)
r = 0.7 * 32767 * np.sqrt(rng.random((F, N))); th = 2 * np.pi * rng.random((F, N))
xr = np.round(r * np.cos(th)).astype(np.int64); xi = np.round(r * np.sin(th)).astype(np.int64)
os.makedirs('simwork', exist_ok=True)
with open('simwork/r22_in.hex', 'w') as f:
    for a, b in zip(xr.reshape(-1), xi.reshape(-1)):
        f.write('%04x%04x\n' % (int(a) & 0xFFFF, int(b) & 0xFFFF))
tb = '''`timescale 1ns/1ps
module tb; reg clock = 0, reset = 1, di_en = 0; reg [15:0] di_re = 0, di_im = 0; wire do_en; wire [15:0] do_re, do_im;
reg [31:0] mem [0:%d]; integer k, fo;
FFT dut(.clock(clock), .reset(reset), .di_en(di_en), .di_re(di_re), .di_im(di_im), .do_en(do_en), .do_re(do_re), .do_im(do_im));
always #5 clock = ~clock;
always @(posedge clock) if (do_en) $fwrite(fo, "%%0d %%0d\\n", $signed(do_re), $signed(do_im));
initial begin $readmemh("simwork/r22_in.hex", mem); fo = $fopen("simwork/r22_out.txt", "w");
 repeat (3) @(posedge clock); #1 reset = 0;
 for (k = 0; k < %d; k = k + 1) begin di_en = 1; di_re = mem[k][31:16]; di_im = mem[k][15:0]; @(posedge clock); #1; end
 di_en = 0; repeat (400) @(posedge clock); $fclose(fo); $finish; end
endmodule
''' % (F * N - 1, F * N)
open('simwork/r22_tb.v', 'w').write(tb)
src = ['ext_r22sdf/' + s for s in ('FFT64.v', 'SdfUnit.v', 'Butterfly.v', 'DelayBuffer.v', 'Multiply.v', 'Twiddle64.v')]
subprocess.check_call(['iverilog', '-g2012', '-o', 'simwork/r22.vvp', 'simwork/r22_tb.v'] + src)
subprocess.check_call(['vvp', '-n', 'simwork/r22.vvp'], stdout=subprocess.DEVNULL)
o = np.loadtxt('simwork/r22_out.txt', dtype=np.int64)
br = np.array([int(format(i, '06b')[::-1], 2) for i in range(N)])
e1 = []; e2 = []; ref = []
for f in range(min(F, len(o) // N)):
    X = np.fft.fft(xr[f] + 1j * xi[f]) / N
    h = o[f * N:(f + 1) * N]; hz = h[:, 0] + 1j * h[:, 1]
    e1.append(hz - X[br]); ref.append(X[br])
    mr, mi = model(xr[f], xi[f], 6); e2.append(np.array(mr) + 1j * np.array(mi) - X[br])
s = lambda e: 10 * np.log10(np.mean(np.abs(np.concatenate(ref)) ** 2) / np.mean(np.abs(np.concatenate(e)) ** 2))
print('frames', len(e1), 'SQNR r22sdf FFT64 (published, MIT) = %.2f dB; ours (V6, 18-bit twiddles) = %.2f dB' % (s(e1), s(e2)))
