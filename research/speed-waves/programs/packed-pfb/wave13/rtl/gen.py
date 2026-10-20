#!/usr/bin/env python3
# RTL generator: one PFB "lane" (one of the 4 parallel sample positions) of a 14-tap polyphase filter branch,
# DSP48E2 primitives instantiated, and the 4-signal x 4-lane slice built from them.
#
#   base     one signal per DSP: a 14-DSP systolic cascade (PCIN/PCOUT), B = W-bit coefficient (W = 18: what
#            SKA-Low builds; W = 11/12: the precision-planned baseline). No fabric logic.
#   packed   two signals per DSP, the whole 14-tap sum stays in ONE cascade:
#              A = a*2^19 + b + 2^18 (27 bits: {a, ~b7, b7 x10, b}; the +2^18 offset keeps a = -128, b < 0 in range)
#              first DSP: P = A*c + C with C = -S*2^18, S = sum of the branch's 14 coefficients (one extra ROM column)
#              P = Ya*2^19 + Yb; low = P[18:0] (signed), high = P[37:19] + P[18]   (one 19-bit fabric adder)
#            Legal when every branch's coefficients satisfy 127*P+128*N <= 2^18-1 and 128*P+127*N <= 2^18
#            (P, N = sums of positive / |negative| taps): then |Ya|, |Yb| < 2^18 for ANY 8-bit samples.
#   pertap   two signals per DSP, each DSP alone (P = A*c - c*2^18 via the C port, -c from the ROM), each tap's two
#            products extracted in fabric and summed by two pipelined adder trees (brief variant (i)).
# Timing convention (all variants): tap t's sample and coefficient arrive t cycles after tap 0's (systolic skew;
# in the real filter this is an address offset into the delay-line and coefficient memories, free);
# the C port (S or -c) arrives one cycle after the A/B ports of its DSP (CREG = 1 aligned with MREG).
import sys, os
T = 14

DSP_PARAMS = dict(A_INPUT='"DIRECT"', B_INPUT='"DIRECT"', AMULTSEL='"A"', BMULTSEL='"B"', USE_MULT='"MULTIPLY"',
                  USE_SIMD='"ONE48"', AREG=1, BREG=1, ACASCREG=1, BCASCREG=1, CREG=1, MREG=1, PREG=1, DREG=1, ADREG=1,
                  OPMODEREG=1, ALUMODEREG=1, INMODEREG=1, CARRYINREG=1, CARRYINSELREG=1)

def dsp(name, A, B, C, pcin, opmode, P, pcout):
    ps = ',\n    '.join(f'.{k}({v})' for k, v in DSP_PARAMS.items())
    return f"""  DSP48E2 #(
    {ps}) {name} (
    .CLK(clk), .A({A}), .B({B}), .C({C}), .D(27'd0), .PCIN({pcin}), .ACIN(30'd0), .BCIN(18'd0),
    .OPMODE(9'b{opmode:09b}), .ALUMODE(4'b0000), .INMODE(5'b00000), .CARRYINSEL(3'b000), .CARRYIN(1'b0),
    .CARRYCASCIN(1'b0), .MULTSIGNIN(1'b0),
    .CEA1(1'b1), .CEA2(1'b1), .CEAD(1'b1), .CEALUMODE(1'b1), .CEB1(1'b1), .CEB2(1'b1), .CEC(1'b1), .CECARRYIN(1'b1),
    .CECTRL(1'b1), .CED(1'b1), .CEINMODE(1'b1), .CEM(1'b1), .CEP(1'b1),
    .RSTA(1'b0), .RSTALLCARRYIN(1'b0), .RSTALUMODE(1'b0), .RSTB(1'b0), .RSTC(1'b0), .RSTCTRL(1'b0), .RSTD(1'b0),
    .RSTINMODE(1'b0), .RSTM(1'b0), .RSTP(1'b0),
    .P({P}), .PCOUT({pcout}));
"""

OP_M = 0b000000101          # X = M, Y = M, Z = 0, W = 0
OP_M_C = 0b000110101        # Z = C
OP_M_PCIN = 0b000010101     # Z = PCIN

def base_lane(W):
    s = [f"// one signal, one lane: 14-DSP cascade, {W}-bit coefficients. y = sum_t c_t x_t (latency {T + 2} from tap 0)",
         f"module base_lane_w{W} (input clk, input [{T*8-1}:0] x, input [{T*W-1}:0] c, output [47:0] y);",
         f"  wire [47:0] pc [0:{T}];", "  assign pc[0] = 48'd0;"]
    for t in range(T):
        A = f"{{{{22{{x[{t*8+7}]}}}}, x[{t*8+7}:{t*8}]}}"
        B = f"{{{{{18-W}{{c[{t*W+W-1}]}}}}, c[{t*W+W-1}:{t*W}]}}" if W < 18 else f"c[{t*W+W-1}:{t*W}]"
        P = "y" if t == T - 1 else ""
        s.append(dsp(f"u{t}", A, B, "48'd0", f"pc[{t}]", OP_M if t == 0 else OP_M_PCIN, P, f"pc[{t+1}]"))
    s.append("endmodule\n")
    return '\n'.join(s)

def packed_lane(W, broken=False):
    # broken: the correction bit is wrong (P[17] instead of P[18])
    cb = 17 if broken else 18
    nm = f"packed_lane_w{W}" + ("_broken" if broken else "")
    s = [f"// two signals (a, b), one lane, ONE 14-DSP cascade. ya = sum a_t c_t, yb = sum b_t c_t (latency {T + 3})",
         f"module {nm} (input clk, input [{T*8-1}:0] a, input [{T*8-1}:0] b, input [{T*W-1}:0] c,",
         f"  input [15:0] negs, output reg [18:0] ya, output reg [18:0] yb);",
         f"  wire [47:0] pc [0:{T}]; wire [47:0] p;", "  assign pc[0] = 48'd0;",
         "  wire [47:0] cin = {{14{negs[15]}}, negs, 18'd0};   // -S * 2^18"]
    for t in range(T):
        at, bt = f"a[{t*8+7}:{t*8}]", f"b[{t*8+7}:{t*8}]"
        A = f"{{{{3{{a[{t*8+7}]}}}}, {at}, ~b[{t*8+7}], {{10{{b[{t*8+7}]}}}}, {bt}}}"
        B = f"{{{{{18-W}{{c[{t*W+W-1}]}}}}, c[{t*W+W-1}:{t*W}]}}"
        s.append(dsp(f"u{t}", A, B, "cin" if t == 0 else "48'd0", f"pc[{t}]", OP_M_C if t == 0 else OP_M_PCIN,
                     "p" if t == T - 1 else "", f"pc[{t+1}]"))
    s.append(f"  always @(posedge clk) begin yb <= p[18:0]; ya <= p[37:19] + p[{cb}]; end")
    s.append("endmodule\n")
    return '\n'.join(s)

def pertap_lane(W):
    s = [f"// two signals, one lane, 14 separate DSPs (P_t = c_t*A_t - c_t*2^18), fabric extraction + adder trees",
         f"module pertap_lane_w{W} (input clk, input [{T*8-1}:0] a, input [{T*8-1}:0] b, input [{T*W-1}:0] c,",
         f"  input [{T*(W+1)-1}:0] negc, output reg [22:0] ya, output reg [22:0] yb);",
         f"  wire [47:0] p [0:{T-1}];"]
    for t in range(T):
        at, bt = f"a[{t*8+7}:{t*8}]", f"b[{t*8+7}:{t*8}]"
        A = f"{{{{3{{a[{t*8+7}]}}}}, {at}, ~b[{t*8+7}], {{10{{b[{t*8+7}]}}}}, {bt}}}"
        B = f"{{{{{18-W}{{c[{t*W+W-1}]}}}}, c[{t*W+W-1}:{t*W}]}}"
        Cn = f"{{{{{48-18-(W+1)}{{negc[{t*(W+1)+W}]}}}}, negc[{t*(W+1)+W}:{t*(W+1)}], 18'd0}}"
        s.append(dsp(f"u{t}", A, B, Cn, "48'd0", OP_M_C, f"p[{t}]", ""))
    # extraction (registered), then adder trees with a register per level
    s.append(f"  reg signed [18:0] lo0 [0:{T-1}]; reg signed [18:0] hi0 [0:{T-1}];")
    s.append("  integer i;")
    s.append(f"  always @(posedge clk) for (i = 0; i < {T}; i = i + 1) begin lo0[i] <= p[i][18:0]; hi0[i] <= p[i][37:19] + p[i][18]; end")
    n, lvl, prev, wdt = T, 0, '0', 19
    while n > 1:
        m = (n + 1) // 2; wdt += 1; lvl += 1
        s.append(f"  reg signed [{wdt-1}:0] lo{lvl} [0:{m-1}]; reg signed [{wdt-1}:0] hi{lvl} [0:{m-1}];")
        body = []
        for j in range(m):
            if 2 * j + 1 < n:
                body.append(f"lo{lvl}[{j}] <= lo{lvl-1}[{2*j}] + lo{lvl-1}[{2*j+1}]; hi{lvl}[{j}] <= hi{lvl-1}[{2*j}] + hi{lvl-1}[{2*j+1}];")
            else:
                body.append(f"lo{lvl}[{j}] <= lo{lvl-1}[{2*j}]; hi{lvl}[{j}] <= hi{lvl-1}[{2*j}];")
        s.append("  always @(posedge clk) begin " + ' '.join(body) + " end")
        n = m
    s.append(f"  always @(posedge clk) begin ya <= hi{lvl}[0]; yb <= lo{lvl}[0]; end")
    s.append("endmodule\n")
    return '\n'.join(s), lvl

def slice_mod(var, W, nsig=4, nlane=4):
    # nsig signals x nlane lanes; coefficients (and S / -c) shared by all signals of a lane
    s = [f"module slice_{var}_w{W} (input clk, input [{nsig*nlane*T*8-1}:0] x, input [{nlane*T*W-1}:0] c,"]
    if var == 'packed': s.append(f"  input [{nlane*16-1}:0] negs,")
    if var == 'pertap': s.append(f"  input [{nlane*T*(W+1)-1}:0] negc,")
    ow = {'base': 48, 'packed': 19, 'pertap': 23}[var]
    s.append(f"  output [{nsig*nlane*ow-1}:0] y);")
    for l in range(nlane):
        cl = f"c[{(l+1)*T*W-1}:{l*T*W}]"
        def xs(sig): k = sig * nlane + l; return f"x[{(k+1)*T*8-1}:{k*T*8}]"
        def ys(sig): k = sig * nlane + l; return f"y[{(k+1)*ow-1}:{k*ow}]"
        if var == 'base':
            for sg in range(nsig):
                s.append(f"  base_lane_w{W} u_l{l}_s{sg} (.clk(clk), .x({xs(sg)}), .c({cl}), .y({ys(sg)}));")
        else:
            for pr in range(nsig // 2):
                extra = f".negs(negs[{l*16+15}:{l*16}])" if var == 'packed' else f".negc(negc[{(l+1)*T*(W+1)-1}:{l*T*(W+1)}])"
                s.append(f"  {var}_lane_w{W} u_l{l}_p{pr} (.clk(clk), .a({xs(2*pr)}), .b({xs(2*pr+1)}), .c({cl}), {extra},"
                         f" .ya({ys(2*pr)}), .yb({ys(2*pr+1)}));")
    s.append("endmodule\n")
    return '\n'.join(s)

if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    for W in (18, 12, 11):
        open(f"{out}/base_lane_w{W}.v", 'w').write(base_lane(W))
        open(f"{out}/slice_base_w{W}.v", 'w').write(slice_mod('base', W))
    for W in (12, 11):
        open(f"{out}/packed_lane_w{W}.v", 'w').write(packed_lane(W))
        open(f"{out}/packed_lane_w{W}_broken.v", 'w').write(packed_lane(W, broken=True))
        open(f"{out}/slice_packed_w{W}.v", 'w').write(slice_mod('packed', W))
        txt, lv = pertap_lane(W)
        open(f"{out}/pertap_lane_w{W}.v", 'w').write(txt)
        open(f"{out}/slice_pertap_w{W}.v", 'w').write(slice_mod('pertap', W))
    print("pertap tree levels:", lv)

def packedg_lane(W, G=7):
    # two signals per DSP, the 14 taps in groups of G taps, one packed cascade per group (each with its own -S_g * 2^18
    # on the first DSP's C port); each group's result extracted, the two fields summed in fabric.
    # Tap t's inputs arrive (t mod G) cycles after the group's first tap (address offsets), so all groups end together.
    ng = (T + G - 1) // G
    ow = 19 + (ng - 1).bit_length()
    s = [f"module packedg{G}_lane_w{W} (input clk, input [{T*8-1}:0] a, input [{T*8-1}:0] b, input [{T*W-1}:0] c,",
         f"  input [{ng*16-1}:0] negs, output reg [{ow-1}:0] ya, output reg [{ow-1}:0] yb);",
         f"  wire [47:0] pc [0:{T}]; wire [47:0] p [0:{ng-1}];"]
    for g in range(ng):
        s.append(f"  wire [47:0] cin{g} = {{{{14{{negs[{g*16+15}]}}}}, negs[{g*16+15}:{g*16}], 18'd0}};")
    for t in range(T):
        g, first, last = t // G, t % G == 0, (t % G == G - 1) or t == T - 1
        A = f"{{{{3{{a[{t*8+7}]}}}}, a[{t*8+7}:{t*8}], ~b[{t*8+7}], {{10{{b[{t*8+7}]}}}}, b[{t*8+7}:{t*8}]}}"
        B = f"{{{{{18-W}{{c[{t*W+W-1}]}}}}, c[{t*W+W-1}:{t*W}]}}"
        s.append(dsp(f"u{t}", A, B, f"cin{g}" if first else "48'd0", "48'd0" if first else f"pc[{t}]",
                     OP_M_C if first else OP_M_PCIN, f"p[{g}]" if last else "", "" if last else f"pc[{t+1}]"))
    s.append(f"  reg signed [18:0] lo [0:{ng-1}]; reg signed [18:0] hi [0:{ng-1}]; integer i;")
    s.append(f"  always @(posedge clk) for (i = 0; i < {ng}; i = i + 1) begin lo[i] <= p[i][18:0]; hi[i] <= p[i][37:19] + p[i][18]; end")
    s.append("  always @(posedge clk) begin ya <= " + ' + '.join(f"hi[{g}]" for g in range(ng)) + "; yb <= "
             + ' + '.join(f"lo[{g}]" for g in range(ng)) + "; end")
    s.append("endmodule\n")
    return '\n'.join(s), ng, ow

def slice_packedg(W, G=7, nsig=4, nlane=4):
    _, ng, ow = packedg_lane(W, G)
    s = [f"module slice_packedg{G}_w{W} (input clk, input [{nsig*nlane*T*8-1}:0] x, input [{nlane*T*W-1}:0] c,",
         f"  input [{nlane*ng*16-1}:0] negs, output [{nsig*nlane*ow-1}:0] y);"]
    for l in range(nlane):
        cl = f"c[{(l+1)*T*W-1}:{l*T*W}]"
        for pr in range(nsig // 2):
            ka, kb = 2 * pr * nlane + l, (2 * pr + 1) * nlane + l
            s.append(f"  packedg{G}_lane_w{W} u_l{l}_p{pr} (.clk(clk), .a(x[{(ka+1)*T*8-1}:{ka*T*8}]), .b(x[{(kb+1)*T*8-1}:{kb*T*8}]),"
                     f" .c({cl}), .negs(negs[{(l+1)*ng*16-1}:{l*ng*16}]), .ya(y[{(ka+1)*ow-1}:{ka*ow}]), .yb(y[{(kb+1)*ow-1}:{kb*ow}]));")
    s.append("endmodule\n")
    return '\n'.join(s)

if __name__ == '__main__':
    for G in (7, 4):
        txt, ng, ow = packedg_lane(12, G)
        open(f"{out}/packedg{G}_lane_w12.v", 'w').write(txt)
        open(f"{out}/slice_packedg{G}_w12.v", 'w').write(slice_packedg(12, G))
