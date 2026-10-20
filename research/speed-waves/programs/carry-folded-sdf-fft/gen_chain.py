# gen_chain.py (as-5): the FMCW range-Doppler chain around two generated FFTs, for area shares.
#   window (Hann ROM, 2 multiplies, round) -> range FFT (fft_<arch>_256) -> corner turn (double-buffered
#   2 x 64 chirps x 256 bins x 32 bits, block RAM) -> Doppler FFT (fft_<arch>_64) -> |.|^2 (2 multiplies, add).
# One complex 16-bit ADC sample per clock, 256 samples per chirp, 64 chirps per frame, continuous.
# Output: power (top 16 of 33 bits). Only for area/clock shares; the FFTs themselves are checked by sim.py.
import sys, math
arch = sys.argv[1]; out = sys.argv[2]
NR, ND = 256, 64
L = []
w = L.append
w('module chain_%s (input clk, input rst, input signed [15:0] adc_re, input signed [15:0] adc_im, output [15:0] pwr);' % arch)
# window
w('  reg [7:0] wn; always @(posedge clk) wn <= rst ? 8\'d0 : wn + 8\'d1;')
w('  reg [15:0] wv; always @(posedge clk) case (wn)')
for k in range(NR):
    v = int(round(32767 * (0.5 - 0.5 * math.cos(2 * math.pi * k / NR))))
    w('    8\'d%d: wv <= 16\'d%d;' % (k, v))
w('  endcase')
w('  reg signed [15:0] a_re, a_im; always @(posedge clk) begin a_re <= adc_re; a_im <= adc_im; end')
w('  reg signed [32:0] p_re, p_im; always @(posedge clk) begin p_re <= a_re * $signed({1\'b0, wv}); p_im <= a_im * $signed({1\'b0, wv}); end')
w('  reg signed [15:0] x_re, x_im; always @(posedge clk) begin x_re <= (p_re + 33\'sd16384) >>> 15; x_im <= (p_im + 33\'sd16384) >>> 15; end')
w('  reg rst2, rst3, rst4; always @(posedge clk) begin rst2 <= rst; rst3 <= rst2; rst4 <= rst3; end')
w('  wire signed [15:0] r_re, r_im;')
w('  fft_%s_256 u_range (.clk(clk), .rst(rst4), .in_re(x_re), .in_im(x_im), .out_re(r_re), .out_im(r_im));' % arch)
# corner turn: write counter (starts after the range FFT latency is irrelevant for area; use a free counter)
w('  reg [14:0] wc; always @(posedge clk) wc <= rst ? 15\'d0 : wc + 15\'d1;  // {bank, chirp[5:0], slot[7:0]}')
w('  reg [31:0] ct0 [0:16383]; reg [31:0] ct1 [0:16383];')
w('  wire [13:0] wa = wc[13:0];')
w('  wire [13:0] ra = {wc[7:0], wc[13:8]};  // read bins in order, 64 chirps each: transpose of the write address')
w('  reg [31:0] q0, q1;')
w('  always @(posedge clk) begin if (!wc[14]) ct0[wa] <= {r_re, r_im}; else ct1[wa] <= {r_re, r_im}; end')
w('  always @(posedge clk) begin q0 <= ct0[ra]; q1 <= ct1[ra]; end')
w('  reg bsel; always @(posedge clk) bsel <= wc[14];')
w('  wire [31:0] q = bsel ? q0 : q1;')
w('  wire signed [15:0] d_re, d_im;')
w('  fft_%s_64 u_doppler (.clk(clk), .rst(rst4), .in_re(q[31:16]), .in_im(q[15:0]), .out_re(d_re), .out_im(d_im));' % arch)
w('  reg signed [31:0] s_re, s_im; always @(posedge clk) begin s_re <= d_re * d_re; s_im <= d_im * d_im; end')
w('  reg [32:0] pw; always @(posedge clk) pw <= s_re + s_im;')
w('  assign pwr = pw[32:17];')
w('endmodule')
open(out, 'w').write('\n'.join(L) + '\n')
