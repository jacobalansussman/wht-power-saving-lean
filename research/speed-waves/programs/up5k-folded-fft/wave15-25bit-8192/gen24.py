#!/usr/bin/env python3
"""gen24.py -- av-5 (wave 15): a FOLDED N-point complex FFT for the iCE40-UP5K with 25-bit data words
(24-bit-class dynamic range; stored as 32-bit words in 2 SPRAMs side by side), for a real-input vibration analyser
(packed pair: two real frames per complex FFT, PSD = (|Z_k|^2 + |Z_{N-k}|^2)/2, as-0 / at-cross H3).

Algorithm (as as-6 ffold.py, widened): in-place radix-2 DIF, every stage halves (V6 rounding), Q14 cosine table over
half a turn, output X/N in bit-reversed order.  Inputs are 16-bit ADC samples, stored as sample << 8 (pre = 8, as
at-cross/h3/bfp.py).  All values stay below 2^24 in magnitude (complex magnitude <= 2^23 sqrt 2), so 25 bits suffice
and nothing wraps.
  s = (a>>1) + (b>>1) + (b&1);  d = (a>>1) - (b>>1)
  Y_re = (2^13 + d_re C + d_im S) >> 14;  Y_im = (2^13 + d_im C - d_re S) >> 14      (exact, no wrap)
  C = ROM[k], S = (k[L-2] ? +1 : -1) ROM[k ^ N/4], ROM[k] = round(16384 cos(pi k / (N/2))), k = (j << stage) mod N/2

Memory: 2 SPRAMs = one 16K x 32 single-port memory; complex word u at addresses 2u (re) and 2u+1 (im).
Slot of 8 clocks per butterfly j (one memory access per clock):
  ph0 R a.re | ph1 R b.re | ph2 R a.im | ph3 W y.re(j-1) | ph4 R b.im | ph5 W s.re(j) | ph6 W s.im(j) | ph7 W y.im(j-1)
DSP version ("dsp", 6 SB_MAC16):
  MS   32-bit adder: C:D = a>>1 (input registers, captured ph1 and ph3), A:B = b>>1 (SPRAM out), CI = b[0] -> s
  MDR  32-bit subtractor, C:D captured ph1 (a.re), result at the end of ph2 -> d.re (held a whole slot)
  MDI  the same, captured ph3, result at the end of ph5 -> d.im
  ML   multiply-accumulate of the LOW 9 bits of d (unsigned, as a positive 16-bit signed A) x ROM: 2^13 + dlo*C +- dlo'*|S|
  MH   multiply-accumulate of the HIGH 16 bits d[24:9] x ROM: dhi*C +- dhi'*|S|          (ML, MH in lockstep)
  MC   32-bit adder: H + (Lo >> 9);  Y = that >> 5  ==  (2^13 + d C + d' S) >> 14 exactly (proof in the report)
  ML/MH use their A/B input registers (A_REG = B_REG = 1), no pipeline register; OLOAD loads the start value.
Usage: gen24.py gen --var dsp --log 13 -o f24.v ; gen24.py sim f24.v --log 13 --frames 4 --work DIR
"""
import sys, os, argparse, subprocess, math
import numpy as np

BIN = '$OSS_CAD_SUITE/bin'
CELLS = '$OSS_CAD_SUITE/share/yosys/ice40/cells_sim.v'


def rom_of(L, tw=14):
    H = 1 << (L - 1)
    return [min(int(round((1 << tw) * math.cos(math.pi * k / H))), 32767 if tw <= 15 else (1 << tw)) for k in range(H)]


def model(xr, xi, L, broken=None, pre=8, W=None, tw=14):
    """bit-true in-place DIF on samples << pre; returns (re, im) in memory order (address u holds X[bitrev(u)]/N << pre).
    W: optional word width with wrap-around (None = exact; the 25-bit hardware never wraps)."""
    N = 1 << L; H = N >> 1
    rom = np.array(rom_of(L, tw), dtype=np.int64)
    if broken == 'tw':
        rom = rom.copy(); rom[H // 3] ^= 1
    re = np.array(xr, dtype=np.int64) << pre; im = np.array(xi, dtype=np.int64) << pre
    j = np.arange(H, dtype=np.int64)
    for s in range(L):
        q = L - 1 - s; h = 1 << q
        a = ((j >> q) << (q + 1)) | (j & (h - 1)); b = a | h
        k = (j << s) & (H - 1)
        ar, ai, br, bi = re[a], im[a], re[b], im[b]
        cr = br & 1; ci = bi & 1
        if broken == 'carry':
            cr = cr * 0
        sr = (ar >> 1) + (br >> 1) + cr; si = (ai >> 1) + (bi >> 1) + ci
        dr = (ar >> 1) - (br >> 1); di = (ai >> 1) - (bi >> 1)
        C = rom[k]; sg = np.where((k >> (L - 2)) & 1 == 1, 1, -1); S = sg * rom[k ^ (H >> 1)]
        rc = 1 << (tw - 1)
        yr = (rc + dr * C + di * S) >> tw; yi = (rc + di * C - dr * S) >> tw
        if W:
            lim = 1 << (W - 1); wr = lambda v: ((v + lim) & (2 * lim - 1)) - lim
            sr, si, yr, yi = wr(sr), wr(si), wr(yr), wr(yi)
        re[a] = sr; im[a] = si; re[b] = yr; im[b] = yi
    return re, im


def bitrev(n, bits):
    r = 0
    for i in range(bits):
        r = (r << 1) | ((n >> i) & 1)
    return r


MACP = ("SB_MAC16 #(.C_REG(1'b%(cr)s), .A_REG(1'b%(ar)s), .B_REG(1'b%(ar)s), .D_REG(1'b%(cr)s), .TOP_8x8_MULT_REG(1'b0), .BOT_8x8_MULT_REG(1'b0),"
        " .PIPELINE_16x16_MULT_REG1(1'b0), .PIPELINE_16x16_MULT_REG2(1'b0), .TOPOUTPUT_SELECT(2'd1), .TOPADDSUB_LOWERINPUT(2'd%(low)s),"
        " .TOPADDSUB_UPPERINPUT(1'b%(up)s), .TOPADDSUB_CARRYSELECT(2'd2), .BOTOUTPUT_SELECT(2'd1), .BOTADDSUB_LOWERINPUT(2'd%(low)s),"
        " .BOTADDSUB_UPPERINPUT(1'b%(up)s), .BOTADDSUB_CARRYSELECT(2'd%(bcs)s), .MODE_8x8(1'b%(m8)s), .A_SIGNED(1'b1), .B_SIGNED(1'b1))\n"
        "   %(name)s (.CLK(clk), .CE(1'b1), .C(%(C)s), .A(%(A)s), .B(%(B)s), .D(%(D)s), .AHOLD(1'b0), .BHOLD(1'b0), .CHOLD(%(ch)s), .DHOLD(%(ch)s),"
        " .IRSTTOP(1'b0), .IRSTBOT(1'b0), .ORSTTOP(1'b0), .ORSTBOT(1'b0), .OLOADTOP(%(ld)s), .OLOADBOT(%(ld)s), .ADDSUBTOP(%(sub)s), .ADDSUBBOT(%(sub)s),"
        " .OHOLDTOP(%(oh)s), .OHOLDBOT(%(oh)s), .CI(%(ci)s), .ACCUMCI(1'b0), .SIGNEXTIN(1'b0), .O(%(O)s));\n")


def gen(var, L):
    N = 1 << L; H = N >> 1
    J = L - 1          # bits of j and k
    TWQ = int(os.environ.get('TW', '14'))
    rom = rom_of(L, TWQ)
    brk_tw = os.environ.get('BRK_TW'); brk_carry = os.environ.get('BRK_CARRY')
    o = []
    w = o.append
    w("""// av-5 (wave 15) folded %(N)d-point complex FFT, 25-bit data in 32-bit words, variant %(var)s. Generated by gen24.py.
// Protocol: while idle, we_in writes {din << 8} (sign-extended) at word address addr_in (2u = re, 2u+1 = im), and dout
// shows the 32-bit word at addr_in one clock later; a pulse on start runs the FFT in place (busy high); afterwards
// word u holds X[bitrev%(L)d(u)] / %(N)d (scaled by 2^8).
module f24 (
  input clk, input start, input we_in, input [%(L)d:0] addr_in, input [15:0] din,
  output [31:0] dout, output reg busy);
  reg [2:0] ph; reg [%(J1)d:0] j; reg [%(L1)d:0] th; reg primed, drain; reg [%(L1)d:0] bprev;
  reg [%(J1)d:0] k, step;
  wire [%(L1)d:0] a, oh; genvar x;
  generate for (x = 0; x < %(L)d; x = x + 1) begin : ga
    if (x == 0) begin
      assign a[0] = th[0] ? j[0] : 1'b0;  assign oh[0] = ~th[0];
    end else if (x == %(L1)d) begin
      assign a[%(L1)d] = th[%(L2)d] ? 1'b0 : j[%(L2)d]; assign oh[%(L1)d] = th[%(L2)d];
    end else begin
      assign a[x] = th[x] ? j[x] : (th[x-1] ? 1'b0 : j[x-1]); assign oh[x] = ~th[x] & th[x-1];
    end
  end endgenerate
  wire [%(L1)d:0] b = a | oh;
  reg [%(L)d:0] caddr;
  always @* case (ph)
    3'd0: caddr = {a, 1'b0};     3'd1: caddr = {b, 1'b0};     3'd2: caddr = {a, 1'b1};  3'd3: caddr = {bprev, 1'b0};
    3'd4: caddr = {b, 1'b1};     3'd5: caddr = {a, 1'b0};     3'd6: caddr = {a, 1'b1};  default: caddr = {bprev, 1'b1};
  endcase
  wire [13:0] maddr = busy ? caddr : addr_in;
  wire wy = ph[1] & ph[0];                                  // ph3, ph7: y of the previous butterfly
  wire ws = (ph == 3'd5) | (ph == 3'd6);                    // s of this butterfly
  wire mwe = busy ? ((wy & primed) | (ws & ~drain)) : we_in;
  wire [31:0] mo;                                           // SPRAM pair output
  assign dout = mo;
""" % dict(N=N, var=var, L=L, L1=L - 1, L2=L - 2, J1=J - 1))
    if var == 'lcs':
        w(LCS % dict(J1=J - 1, J2=J - 2))
    if var == 'beh':
        w(BEH % dict(J1=J - 1, J2=J - 2))
    if var == 'dsp3':
        w("""
  // ---- dsp3: the butterfly in LOGIC CELLS (25-bit adder and subtractor, a held in ar, d.re and d.im in registers
  // with the same timing as dsp's MDR/MDI outputs), the twiddle exactly as dsp (ML, MH, MC): 3 SB_MAC16
  reg [24:0] ar, dr25, di25, so25;
  wire [24:0] mo25 = mo[24:0];
  wire [24:0] ha = {ar[24], ar[24:1]}, hb25 = {mo25[24], mo25[24:1]};
  wire [24:0] sv = ha + hb25 + mo25[0];
  wire [24:0] dv = ha - hb25;
  always @(posedge clk) if (busy) begin
    if (ph == 3'd1 || ph == 3'd3) ar <= mo25;
    if (ph == 3'd2 || ph == 3'd5) so25 <= sv;
    if (ph == 3'd2) dr25 <= dv;
    if (ph == 3'd5) di25 <= dv;
  end
  wire [31:0] so = {{7{so25[24]}}, so25}, dre = {{7{dr25[24]}}, dr25}, dim = {{7{di25[24]}}, di25};
  wire [31:0] ml, mh, mc;
""")
    if var in ('dsp', 'dsp3'):
      if var == 'dsp':
        w("""
  // ---- butterfly: 3 SB_MAC16 as 32-bit adders. hb = b >> 1 straight from the SPRAMs; a >> 1 held in C/D input registers
  wire [31:0] hb = {mo[31], mo[31:1]};
  wire [31:0] so, dre, dim, ml, mh, mc;
""")
      if var == 'dsp':
        w('  ' + MACP % dict(cr=1, ar=0, low=0, up=1, bcs=3, m8=1, name='ms', C='hb[31:16]', A='hb[31:16]', B='hb[15:0]', D='hb[15:0]',
                             ch="~(ph == 3'd1 || ph == 3'd3)", ld="1'b0", sub="1'b0", oh="~(ph == 3'd2 || ph == 3'd5)",
                             ci=("1'b0" if brk_carry else 'mo[0]'), O='so'))
      if var == 'dsp':
        w('  ' + MACP % dict(cr=1, ar=0, low=0, up=1, bcs=0, m8=1, name='mdr', C='hb[31:16]', A='hb[31:16]', B='hb[15:0]', D='hb[15:0]',
                             ch="~(ph == 3'd1)", ld="1'b0", sub="1'b1", oh="~(ph == 3'd2)", ci="1'b0", O='dre'))
      if var == 'dsp':
        w('  ' + MACP % dict(cr=1, ar=0, low=0, up=1, bcs=0, m8=1, name='mdi', C='hb[31:16]', A='hb[31:16]', B='hb[15:0]', D='hb[15:0]',
                             ch="~(ph == 3'd3)", ld="1'b0", sub="1'b1", oh="~(ph == 3'd5)", ci="1'b0", O='dim'))
      if True:
        w("""
  // ---- twiddle: ML (low 9 bits of d) and MH (d[24:9]) in lockstep, A/B input registers on (operands presented one clock
  // before the product).  load ph6, +d.re*C ph7, +-d.im*S ph0, combine re ph1 (and load), +d.im*C ph2, -+d.re*S ph3,
  // combine im ph4.  Operand pins: ph6 d.re, ph7 d.im, ph1 d.im, ph2 d.re.  ROM address: ph5 k, ph6 k^N/4, ph0 k, ph1 k^N/4.
  reg [15:0] romq; wire [%(J1)d:0] radr = {k[%(J1)d] ^ (ph == 3'd6 || ph == 3'd1), k[%(J2)d:0]};
  wire selim = (ph == 3'd7) || (ph == 3'd1);
  wire [31:0] dsel = selim ? dim : dre;
  wire tld = (ph == 3'd6) || (ph == 3'd1);
  wire thold = (ph == 3'd4) || (ph == 3'd5);
  wire tsub = (ph == 3'd0 && ~k[%(J1)d]) || (ph == 3'd3 && k[%(J1)d]);
""" % dict(J1=J - 1, J2=J - 2))
        w('  ' + MACP % dict(cr=0, ar=1, low=2, up=0, bcs=0, m8=0, name='mlo', C="16'h0000", A="{7'b0, dsel[8:0]}", B='romq', D=("16'h%04x" % (1 << (TWQ - 1))),
                             ch="1'b0", ld='tld', sub='tsub', oh='thold', ci="1'b0", O='ml'))
        w('  ' + MACP % dict(cr=0, ar=1, low=2, up=0, bcs=0, m8=0, name='mhi', C="16'h0000", A='dsel[24:9]', B='romq', D="16'h0000",
                             ch="1'b0", ld='tld', sub='tsub', oh='thold', ci="1'b0", O='mh'))
        w('  ' + MACP % dict(cr=0, ar=0, low=0, up=1, bcs=0, m8=1, name='mcb', C='mh[31:16]', A="{{9{ml[31]}}, ml[31:25]}", B='ml[24:9]', D='mh[15:0]',
                             ch="1'b0", ld="1'b0", sub="1'b0", oh="~(ph == 3'd1 || ph == 3'd4)", ci="1'b0", O='mc'))
        w("""
  wire [31:0] yw = {{7{mc[%(T)d]}}, mc[%(T)d:%(B)d]};
""" % dict(T=TWQ + 15, B=TWQ - 9) + """  wire [31:0] dext = {{8{din[15]}}, din, 8'h00};
  wire [31:0] mdin = busy ? (wy ? yw : so) : dext;
""")
    w("""
  SB_SPRAM256KA spl (.ADDRESS(maddr), .DATAIN(mdin[15:0]), .MASKWREN(4'b1111), .WREN(mwe), .CHIPSELECT(1'b1),
                     .CLOCK(clk), .STANDBY(1'b0), .SLEEP(1'b0), .POWEROFF(1'b1), .DATAOUT(mo[15:0]));
  SB_SPRAM256KA sph (.ADDRESS(maddr), .DATAIN(mdin[31:16]), .MASKWREN(4'b1111), .WREN(mwe), .CHIPSELECT(1'b1),
                     .CLOCK(clk), .STANDBY(1'b0), .SLEEP(1'b0), .POWEROFF(1'b1), .DATAOUT(mo[31:16]));
  // twiddle ROM: Q14 cosine over half a turn (%(H)d x 16 bits; block RAM)
  reg [15:0] rom [0:%(H1)d];
  initial begin
%(ri)s  end
  always @(posedge clk) romq <= rom[radr];
  initial busy = 0;
  always @(posedge clk) begin
    if (!busy) begin
      if (start) begin busy <= 1; ph <= 0; j <= 0; th <= {1'b0, {%(L1)d{1'b1}}}; primed <= 0; drain <= 0; step <= 1; end
    end else begin
      ph <= ph + 3'd1;
      if (ph == 3'd4) k <= (j == 0) ? 0 : k + step;
      if (ph == 3'd7) begin
        bprev <= b; primed <= 1; j <= j + 1'b1;
        if (&j) begin th <= {1'b0, th[%(L1)d:1]}; step <= {step[%(J2)d:0], 1'b0}; end
        if (&j && th == 0) drain <= 1;
        if (drain) busy <= 0;
      end
    end
  end
endmodule
""" % dict(H=H, H1=H - 1, L1=L - 1, J2=J - 2,
           ri=''.join("    rom[%d] = 16'h%04x;\n" % (i, (v ^ (1 if (brk_tw and i == H // 3) else 0)) & 0xFFFF) for i, v in enumerate(rom))))
    v = ''.join(o)
    if var == 'lcs':
        reps = [("reg [2:0] ph;", "reg [5:0] ph;"),
                ("    3'd0: caddr = {a, 1'b0};     3'd1: caddr = {b, 1'b0};     3'd2: caddr = {a, 1'b1};  3'd3: caddr = {bprev, 1'b0};\n"
                 "    3'd4: caddr = {b, 1'b1};     3'd5: caddr = {a, 1'b0};     3'd6: caddr = {a, 1'b1};  default: caddr = {bprev, 1'b1};",
                 "    6'd0: caddr = {a, 1'b0};     6'd1: caddr = {b, 1'b0};     6'd2: caddr = {a, 1'b1};  6'd3: caddr = {bprev, 1'b1};\n"
                 "    6'd4: caddr = {b, 1'b1};     6'd5: caddr = {a, 1'b0};     6'd6: caddr = {a, 1'b1};  default: caddr = {bprev, 1'b0};"),
                ("wire wy = ph[1] & ph[0];", "wire wy = (ph == 6'd3) || (ph == 6'd7);"),
                ("wire ws = (ph == 3'd5) | (ph == 3'd6);", "wire ws = (ph == 6'd5) | (ph == 6'd6);"),
                ("ph <= ph + 3'd1;", "ph <= (ph == 6'd39) ? 6'd0 : ph + 6'd1;"),
                ("if (ph == 3'd4) k <= (j == 0) ? 0 : k + step;", "if (ph == 6'd1) k <= (j == 0) ? 0 : k + step;"),
                ("if (ph == 3'd7) begin", "if (ph == 6'd39) begin")]
        for a, b in reps:
            assert v.count(a) == 1, a
            v = v.replace(a, b)
    if var == 'beh':
        reps = [("3'd3: caddr = {bprev, 1'b0};", "3'd3: caddr = {bprev, 1'b1};"), ("default: caddr = {bprev, 1'b1};", "default: caddr = {bprev, 1'b0};"),
                ("if (ph == 3'd4) k <= (j == 0) ? 0 : k + step;", "if (ph == 3'd2) k <= (j == 0) ? 0 : k + step;")]
        for a, b in reps:
            assert v.count(a) == 1, a
            v = v.replace(a, b)
    return v



LCS = r"""
  // ---- lcs: LOGIC-ONLY baseline at the same widths (25-bit data): butterfly adders in logic cells and a SERIAL
  // twiddle multiply (as-6's fold_lcs, widened): one bit of each twiddle per clock, both products of a component at
  // once, H <- (H + b1_i*A1 + b2_i*A2 + c) >> 1, rounding carry at step 13, bits 14, 15 kept.  40-clock slot:
  // ph0..ph7 memory as the DSP design (but y.im is written at ph3 and y.re at ph7), ROM k at ph4 (C -> B1 end ph5),
  // k^N/4 at ph5 (S -> B2 end ph6), H cleared ph6, re steps ph7..22, Y_re kept + swap ph23 (H cleared), im steps ph24..39.
  reg [24:0] ar, dreg, P, Q, so, yr; reg [15:0] B1, B2, romq; reg sk, t14, t15; reg signed [26:0] H;
  wire [24:0] mo25 = mo[24:0];
  wire [24:0] ha = {ar[24], ar[24:1]}, hb = {mo25[24], mo25[24:1]};
  wire [24:0] sv = ha + hb + mo25[0];
  wire [24:0] dv = ha - hb;
  wire [%(J1)d:0] radr = {k[%(J1)d] ^ (ph == 6'd5), k[%(J2)d:0]};
  wire stR = ph >= 6'd7 && ph <= 6'd22, stI = ph >= 6'd24 && ph <= 6'd39;
  wire [5:0] im7 = ph - 6'd7, im24 = ph - 6'd24;
  wire [3:0] i = stR ? im7[3:0] : im24[3:0];
  wire neg1 = (i == 4'd15), neg2 = (stR ? ~sk : sk) ^ (i == 4'd15);
  wire [27:0] t1 = {28{B1[0]}} & ({{3{P[24]}}, P} ^ {28{neg1}});
  wire [27:0] t2 = {28{B2[0]}} & ({{3{Q[24]}}, Q} ^ {28{neg2}});
  wire c1 = (B1[0] & neg1) | (i == 4'd13), c2 = B2[0] & neg2;
  wire [27:0] s1 = {H[26], H} + t1 + c1;
  wire [27:0] s2 = s1 + t2 + c2;
  wire [24:0] yi = {H[22:0], t15, t14};
  wire [31:0] dext = {{8{din[15]}}, din, 8'h00};
  wire [24:0] wsel = (ph == 6'd3) ? yi : (ph == 6'd7) ? yr : so;
  wire [31:0] mdin = busy ? {{7{wsel[24]}}, wsel} : dext;
  always @(posedge clk) begin
    if (busy) begin
      if (ph == 6'd1 || ph == 6'd3) ar <= mo25;
      if (ph == 6'd2 || ph == 6'd5) so <= sv;
      if (ph == 6'd2) dreg <= dv;
      if (ph == 6'd5) begin P <= dreg; Q <= dv; B1 <= romq; sk <= k[%(J1)d]; end
      if (ph == 6'd6) B2 <= romq;
      if (ph == 6'd23) begin P <= Q; Q <= P; yr <= yi; end
      if (ph == 6'd6 || ph == 6'd23) H <= 0;
      else if (stR || stI) begin
        H <= s2[27:1]; B1 <= {B1[0], B1[15:1]}; B2 <= {B2[0], B2[15:1]};
        if (i == 4'd14) t14 <= s2[0];
        if (i == 4'd15) t15 <= s2[0];
      end
    end
  end
"""

BEH = r"""
  // ---- beh: the OBVIOUS design at the same widths and the same 8-clock slot: butterfly in logic cells, ONE behavioural
  // multiply-accumulate (25 x 16 product, 41-bit accumulator; synth_ice40 -dsp maps the multiply to SB_MAC16s), one
  // product per clock: ph6 2^13 + d.re*C, ph7 +-d.im*S (Y_re kept at ph0), ph0 2^13 + d.im*C, ph1 -+d.re*S (Y_im in acc
  // ph2..ph6).  y.im written at ph3, y.re at ph7.  ROM: k at ph5, k^N/4 at ph6, k at ph7, k^N/4 at ph0.
  reg [24:0] ar, dreg, qreg, so, yr; reg [15:0] romq;
  wire [24:0] mo25 = mo[24:0];
  wire [24:0] ha = {ar[24], ar[24:1]}, hb = {mo25[24], mo25[24:1]};
  wire [24:0] sv = ha + hb + mo25[0];
  wire [24:0] dv = ha - hb;
  wire [%(J1)d:0] radr = {k[%(J1)d] ^ (ph == 3'd6 || ph == 3'd0), k[%(J2)d:0]};
  wire [24:0] mA = (ph == 3'd6 || ph == 3'd1) ? dreg : qreg;
  wire ldm = (ph == 3'd6 || ph == 3'd0);
  wire subm = (ph == 3'd7) ? ~k[%(J1)d] : (ph == 3'd1) ? k[%(J1)d] : 1'b0;
  reg signed [40:0] acc;
  wire signed [40:0] pm = $signed(mA) * $signed(romq);
  wire signed [40:0] base = ldm ? 41'sd8192 : acc;
  wire [24:0] yi = acc[38:14];
  wire [31:0] dext = {{8{din[15]}}, din, 8'h00};
  wire [24:0] wsel = (ph == 3'd3) ? yi : (ph == 3'd7) ? yr : so;
  wire [31:0] mdin = busy ? {{7{wsel[24]}}, wsel} : dext;
  always @(posedge clk) begin
    if (busy) begin
      if (ph == 3'd1 || ph == 3'd3) ar <= mo25;
      if (ph == 3'd2 || ph == 3'd5) so <= sv;
      if (ph == 3'd2) dreg <= dv;
      if (ph == 3'd5) qreg <= dv;
      if (ph == 3'd0) yr <= acc[38:14];
      if (ph == 3'd6 || ph == 3'd7 || ph == 3'd0 || ph == 3'd1) acc <= subm ? base - pm : base + pm;
    end
  end
"""


TB = r"""
`timescale 1ns/1ps
module tb;
  reg clk = 0; always #5 clk = ~clk;
  reg start = 0, we = 0; reg [%(L)d:0] ad = 0; reg [15:0] di = 0;
  wire [31:0] o; wire busy;
  %(top)s dut (.clk(clk), .start(start), .we_in(we), .addr_in(ad), .din(di), .dout(o), .busy(busy));
  reg [15:0] vin [0:%(nw)d]; integer fr, i, out, cyc;
  initial begin
    $readmemh("in.hex", vin);
    out = $fopen("out.txt", "w");
    for (fr = 0; fr < %(nf)d; fr = fr + 1) begin
      @(negedge clk);
      for (i = 0; i < %(N2)d; i = i + 1) begin
        we = 1; ad = i; di = vin[fr*%(N2)d+i]; @(negedge clk);
      end
      we = 0; start = 1; @(negedge clk); start = 0; cyc = 1;
      while (busy) begin @(negedge clk); cyc = cyc + 1; end
      $fwrite(out, "C %%0d\n", cyc);
      ad = 0; @(negedge clk);
      for (i = 0; i < %(N2)d; i = i + 1) begin
        $fwrite(out, "%%08x\n", o);
        ad = i + 1; @(negedge clk);
      end
    end
    $fclose(out); $finish;
  end
endmodule
"""


def frames(nf, L, seed=1):
    N = 1 << L
    rng = np.random.default_rng(seed)
    t = np.arange(N)
    F = []
    F.append((np.full(N, 32767), np.zeros(N, int)))                                            # full-scale DC
    F.append((np.full(N, -32768), np.full(N, -32768)))                                         # all minimum
    F.append((np.where(t % 2 == 0, 32767, -32768), np.where(t % 2 == 0, -32768, 32767)))       # alternating extremes
    imp = np.zeros(N, int); imp[0] = 32767; F.append((imp, np.zeros(N, int)))                  # impulse
    F.append((rng.choice([-32768, 32767, -1, 0, 1], N), rng.choice([-32768, 32767, -1, 0, 1], N)))
    F.append((np.round(32767 * np.cos(2 * np.pi * 17 * t / N)).astype(int), np.round(32767 * np.sin(2 * np.pi * 17 * t / N)).astype(int)))
    F.append((np.round(32767 * np.cos(2 * np.pi * (N / 2 - 0.5) * t / N)).astype(int), np.zeros(N, int)))   # bin-edge tone
    # packed pair: two real Hann-windowed records (tone at a bin edge + noise) in re and im
    w = 0.5 - 0.5 * np.cos(2 * np.pi * t / N)
    for _ in range(2):
        f1 = rng.uniform(10, N / 2 - 10) + 0.5; f2 = rng.uniform(10, N / 2 - 10) + 0.5
        x = np.clip(np.round(w * (32767 * 10 ** (-3 / 20) * np.cos(2 * np.pi * f1 * t / N + rng.uniform(0, 6.3)) + rng.normal(0, 30, N))), -32768, 32767)
        y = np.clip(np.round(w * (32767 * 10 ** (-3 / 20) * np.cos(2 * np.pi * f2 * t / N + rng.uniform(0, 6.3)) + rng.normal(0, 30, N))), -32768, 32767)
        F.append((x.astype(int), y.astype(int)))
    while len(F) < nf:
        if len(F) % 2:
            F.append((rng.integers(-32768, 32768, N), rng.integers(-32768, 32768, N)))
        else:
            g = np.clip(np.round(rng.normal(0, 8192, (2, N))), -32768, 32767).astype(int)
            F.append((g[0], g[1]))
    return F[:nf]


def sim(design, L, nf, work, top='f24', broken=None, skip=0):
    N = 1 << L
    work = os.path.abspath(work); os.makedirs(work, exist_ok=True)
    F = frames(nf + skip, L)[skip:]
    with open(os.path.join(work, 'in.hex'), 'w') as f:
        for xr, xi in F:
            for a, b in zip(xr, xi):
                f.write('%04x\n%04x\n' % (int(a) & 0xFFFF, int(b) & 0xFFFF))
    with open(os.path.join(work, 'tb.v'), 'w') as f:
        f.write(TB % dict(top=top, nw=nf * 2 * N - 1, nf=nf, L=L, N2=2 * N))
    exe = os.path.join(work, 'a.out')
    srcs = [os.path.abspath(design), os.path.join(work, 'tb.v'), CELLS]
    r = subprocess.run([BIN + '/iverilog', '-g2012', '-o', exe, '-DNO_ICE40_DEFAULT_ASSIGNMENTS'] + srcs,
                       capture_output=True, text=True, cwd=work)
    if r.returncode:
        print(r.stderr[-3000:]); return False
    r = subprocess.run([BIN + '/vvp', '-n', exe], capture_output=True, text=True, cwd=work)
    lines = open(os.path.join(work, 'out.txt')).read().split('\n')
    bad = 0; badfr = 0; pos = 0; cyc = []
    for fi, (xr, xi) in enumerate(F):
        cyc.append(int(lines[pos].split()[1])); pos += 1
        mr, mi = model(xr, xi, L, broken=broken, tw=int(os.environ.get('TW', '14')))
        got = lines[pos:pos + 2 * N]; pos += 2 * N
        nb = 0
        for u in range(N):
            vr, vi = got[2 * u], got[2 * u + 1]
            if 'x' in vr + vi or 'z' in vr + vi:
                nb += 1; continue
            if int(vr, 16) != (int(mr[u]) & 0xFFFFFFFF) or int(vi, 16) != (int(mi[u]) & 0xFFFFFFFF):
                nb += 1
        bad += nb; badfr += (nb > 0)
        if nb and badfr <= 3:
            print('frame %d: %d mismatching samples' % (fi + skip, nb))
    print('N %d frames %d (skip %d) samples %d mismatching samples %d frames with errors %d clocks per frame %s' % (
        N, nf, skip, nf * N, bad, badfr, sorted(set(cyc))))
    print('RESULT', 'PASS' if bad == 0 else 'FAIL')
    return bad == 0


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd'); ap.add_argument('design', nargs='?')
    ap.add_argument('--var', default='dsp'); ap.add_argument('-o'); ap.add_argument('--frames', type=int, default=4)
    ap.add_argument('--log', type=int, default=13)
    ap.add_argument('--work', default='sim/w'); ap.add_argument('--top', default='f24'); ap.add_argument('--broken'); ap.add_argument('--skip', type=int, default=0)
    a = ap.parse_args()
    if a.cmd == 'gen':
        open(a.o, 'w').write(gen(a.var, a.log))
    elif a.cmd == 'sim':
        sim(a.design, a.log, a.frames, a.work, top=a.top, broken=a.broken, skip=a.skip)
