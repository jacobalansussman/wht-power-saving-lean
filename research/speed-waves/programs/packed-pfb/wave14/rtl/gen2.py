#!/usr/bin/env python3
# as-4 RTL generator: packed PFB lane with GRADED GROUPS (taps grouped by tap index, one packed DSP cascade per group,
# its own exact coefficient-sum guard, its own extraction; group results added in fabric with a fixed left shift).
#   FAM = e2 (UltraScale DSP48E2, 27x18, field s = 19) or e1 (7-series DSP48E1, 25x18, field s = 17)
#   A = a*2^s + b + 2^(s-1)  = {sign ext, a, ~b7, {s-9}{b7}, b}; first DSP of group g adds C = -S_g * 2^(s-1)
#   group g: P_g = Ya_g*2^s + Yb_g;  yb_g = P_g[s-1:0] (signed), ya_g = P_g[2s-1:s] + P_g[s-1]
#   ya = sum_g ya_g << sh_g, yb = sum_g yb_g << sh_g     (= sum_t a_t * c_raw_t * 2^sh(t))
# Legal when every group's raw coefficients in every branch satisfy 127P + 128N <= 2^(s-1)-1, 128P + 127N <= 2^(s-1).
# Timing convention (as ap-4): tap at position j of group g gets its sample/coefficient (j + Gmax - n_g) cycles after
# cycle 0 (address offsets into the delay-line / coefficient memories); the C value one cycle after its DSP's A/B.
#   gen2.py FAM SCHEME W [broken]   -> glane_<fam>_<scheme>_w<W>[_broken].v and gslice_<fam>_<scheme>_w<W>.v
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'coef'))
T = 14
def schemes():
    src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'coef', 'gq.py')).read()
    i = src.index('SCHEMES = {'); j = src.index('}\n', i)
    ns = {}; exec(src[i:j + 2], ns); return ns['SCHEMES']

def order(scheme):
    # cascade order: groups in the scheme's order, taps in increasing index inside a group
    return [(g, sorted(e[0]), e[1]) for g, e in enumerate(scheme)]
def unpacked(scheme):
    # groups marked 'U': one product per DSP (a and b in separate cascades), no field, no guard
    return {g for g, e in enumerate(scheme) if len(e) > 2}

E2P = dict(A_INPUT='"DIRECT"', B_INPUT='"DIRECT"', AMULTSEL='"A"', BMULTSEL='"B"', USE_MULT='"MULTIPLY"',
           USE_SIMD='"ONE48"', AREG=1, BREG=1, ACASCREG=1, BCASCREG=1, CREG=1, MREG=1, PREG=1, DREG=1, ADREG=1,
           OPMODEREG=1, ALUMODEREG=1, INMODEREG=1, CARRYINREG=1, CARRYINSELREG=1)
E1P = dict(A_INPUT='"DIRECT"', B_INPUT='"DIRECT"', USE_DPORT='"FALSE"', USE_MULT='"MULTIPLY"', USE_SIMD='"ONE48"',
           USE_PATTERN_DETECT='"NO_PATDET"', AREG=1, BREG=1, ACASCREG=1, BCASCREG=1, CREG=1, MREG=1, PREG=1, DREG=1,
           ADREG=1, OPMODEREG=1, ALUMODEREG=1, INMODEREG=1, CARRYINREG=1, CARRYINSELREG=1)

def dsp(fam, name, A, B, C, pcin, zsel, P, pcout):
    # zsel: 'C' (Z = C), 'PCIN' (Z = PCIN); X = Y = M
    if fam == 'e2':
        op = {'C': 0b000110101, 'PCIN': 0b000010101}[zsel]
        ps = ',\n    '.join(f'.{k}({v})' for k, v in E2P.items())
        return f"""  DSP48E2 #(
    {ps}) {name} (
    .CLK(clk), .A({A}), .B({B}), .C({C}), .D(27'd0), .PCIN({pcin}), .ACIN(30'd0), .BCIN(18'd0),
    .OPMODE(9'b{op:09b}), .ALUMODE(4'b0000), .INMODE(5'b00000), .CARRYINSEL(3'b000), .CARRYIN(1'b0),
    .CARRYCASCIN(1'b0), .MULTSIGNIN(1'b0),
    .CEA1(1'b1), .CEA2(1'b1), .CEAD(1'b1), .CEALUMODE(1'b1), .CEB1(1'b1), .CEB2(1'b1), .CEC(1'b1), .CECARRYIN(1'b1),
    .CECTRL(1'b1), .CED(1'b1), .CEINMODE(1'b1), .CEM(1'b1), .CEP(1'b1),
    .RSTA(1'b0), .RSTALLCARRYIN(1'b0), .RSTALUMODE(1'b0), .RSTB(1'b0), .RSTC(1'b0), .RSTCTRL(1'b0), .RSTD(1'b0),
    .RSTINMODE(1'b0), .RSTM(1'b0), .RSTP(1'b0),
    .P({P}), .PCOUT({pcout}));
"""
    op = {'C': 0b0110101, 'PCIN': 0b0010101}[zsel]
    ps = ',\n    '.join(f'.{k}({v})' for k, v in E1P.items())
    return f"""  DSP48E1 #(
    {ps}) {name} (
    .CLK(clk), .A({A}), .B({B}), .C({C}), .D(25'd0), .PCIN({pcin}), .ACIN(30'd0), .BCIN(18'd0),
    .OPMODE(7'b{op:07b}), .ALUMODE(4'b0000), .INMODE(5'b00000), .CARRYINSEL(3'b000), .CARRYIN(1'b0),
    .CARRYCASCIN(1'b0), .MULTSIGNIN(1'b0),
    .CEA1(1'b1), .CEA2(1'b1), .CEAD(1'b1), .CEALUMODE(1'b1), .CEB1(1'b1), .CEB2(1'b1), .CEC(1'b1), .CECARRYIN(1'b1),
    .CECTRL(1'b1), .CED(1'b1), .CEINMODE(1'b1), .CEM(1'b1), .CEP(1'b1),
    .RSTA(1'b0), .RSTALLCARRYIN(1'b0), .RSTALUMODE(1'b0), .RSTB(1'b0), .RSTC(1'b0), .RSTCTRL(1'b0), .RSTD(1'b0),
    .RSTINMODE(1'b0), .RSTM(1'b0), .RSTP(1'b0),
    .P({P}), .PCOUT({pcout}));
"""

def out_width(fam, scheme, W=12):
    s = 19 if fam == 'e2' else 17
    ng = len(scheme); msh = max(e[1] for e in scheme)
    w = s + msh + (ng - 1).bit_length()
    if unpacked(scheme): w = max(w, W + 8 + 2) + 1
    return w

def lane(fam, sname, W, broken=False):
    sc = schemes()[sname]; s = 19 if fam == 'e2' else 17; aw = 27 if fam == 'e2' else 25
    groups = order(sc); ng = len(groups); ow = out_width(fam, sc, W); UG = unpacked(sc)
    cb = s - 2 if broken else s - 1
    nm = f"glane_{fam}_{sname}_w{W}" + ("_broken" if broken else "")
    L = [f"// as-4 graded-group packed PFB lane: {fam}, scheme {sname} = {sc}, s = {s}, raw coefficients W = {W} bits",
         f"module {nm} (input clk, input [{T*8-1}:0] a, input [{T*8-1}:0] b, input [{T*W-1}:0] c,",
         f"  input [{ng*16-1}:0] negs, output reg signed [{ow-1}:0] ya, output reg signed [{ow-1}:0] yb);",
         f"  wire [47:0] pc [0:{2*T}]; wire [47:0] p [0:{ng-1}];"]
    for g in range(ng):
        L.append(f"  wire [47:0] cin{g} = {{{{{48-16-(s-1)}{{negs[{g*16+15}]}}}}, negs[{g*16+15}:{g*16}], {s-1}'d0}};   // -S_g * 2^{s-1}")
    k = 0
    for g, taps, sh in groups:
        if g in UG:   # one product per DSP: a cascade for a, a cascade for b
            for sig in ('a', 'b'):
                L.append(f"  wire [47:0] pu{g}{sig} [0:{len(taps)}]; wire [47:0] pU{g}{sig};")
                for j, t in enumerate(taps):
                    A = f"{{{{22{{{sig}[{t*8+7}]}}}}, {sig}[{t*8+7}:{t*8}]}}"
                    B = f"{{{{{18-W}{{c[{t*W+W-1}]}}}}, c[{t*W+W-1}:{t*W}]}}" if W < 18 else f"c[{t*W+W-1}:{t*W}]"
                    last = j == len(taps) - 1
                    L.append(dsp(fam, f"u{k}_t{t}{sig}", A, B, "48'd0", "48'd0" if j == 0 else f"pu{g}{sig}[{j}]", 'C' if j == 0 else 'PCIN',
                                 f"pU{g}{sig}" if last else "", "" if last else f"pu{g}{sig}[{j+1}]"))
                    k += 1
            continue
        for j, t in enumerate(taps):
            first, last = j == 0, j == len(taps) - 1
            A = (f"{{{{{30-aw+(0 if fam=='e1' else 0)}{{a[{t*8+7}]}}}}, a[{t*8+7}:{t*8}], ~b[{t*8+7}], {{{s-9}{{b[{t*8+7}]}}}}, b[{t*8+7}:{t*8}]}}")
            B = f"{{{{{18-W}{{c[{t*W+W-1}]}}}}, c[{t*W+W-1}:{t*W}]}}" if W < 18 else f"c[{t*W+W-1}:{t*W}]"
            L.append(dsp(fam, f"u{k}_t{t}", A, B, f"cin{g}" if first else "48'd0", "48'd0" if first else f"pc[{k}]",
                         'C' if first else 'PCIN', f"p[{g}]" if last else "", "" if last else f"pc[{k+1}]"))
            k += 1
    L.append(f"  reg signed [{s-1}:0] lo [0:{ng-1}]; reg signed [{s-1}:0] hi [0:{ng-1}]; integer i;")
    PG = [g for g, _, _ in groups if g not in UG]
    L.append("  always @(posedge clk) begin " + ' '.join(f"lo[{i}] <= p[{i}][{s-1}:0]; hi[{i}] <= p[{i}][{2*s-1}:{s}] + p[{i}][{cb}];" for i in PG) + " end")
    for g in UG:
        L.append(f"  reg signed [{ow-1}:0] uA{g}, uB{g}; always @(posedge clk) begin uA{g} <= pU{g}a[{ow-1}:0]; uB{g} <= pU{g}b[{ow-1}:0]; end")
    def term(arr, g, sh):
        if g in UG: return f"uA{g}" if arr == 'hi' else f"uB{g}"
        return f"({{{{{ow-s}{{{arr}[{g}][{s-1}]}}}}, {arr}[{g}]}} << {sh})" if sh else f"{{{{{ow-s}{{{arr}[{g}][{s-1}]}}}}, {arr}[{g}]}}"
    if ng <= 2:
        L.append("  always @(posedge clk) begin ya <= " + ' + '.join(term('hi', g, sh) for g, _, sh in groups) + ";")
        L.append("    yb <= " + ' + '.join(term('lo', g, sh) for g, _, sh in groups) + "; end")
    else:   # pipelined pairwise adder tree, one register per level (one 2-input adder per level)
        for arr, out in (('hi', 'ya'), ('lo', 'yb')):
            cur = [term(arr, g, sh) for g, _, sh in groups]; lv = 0
            while len(cur) > 2:
                lv += 1; nxt = []
                for j in range(0, len(cur), 2):
                    nm_ = f"{arr}_l{lv}_{j // 2}"
                    L.append(f"  reg signed [{ow-1}:0] {nm_}; always @(posedge clk) {nm_} <= " + (' + '.join(cur[j:j + 2])) + ";")
                    nxt.append(nm_)
                cur = nxt
            L.append(f"  always @(posedge clk) {out} <= " + ' + '.join(cur) + ";")
    L.append("endmodule\n")
    return '\n'.join(L)

def slice_mod(fam, sname, W, nsig=4, nlane=4):
    sc = schemes()[sname]; ng = len(sc); ow = out_width(fam, sc, W)
    L = [f"module gslice_{fam}_{sname}_w{W} (input clk, input [{nsig*nlane*T*8-1}:0] x, input [{nlane*T*W-1}:0] c,",
         f"  input [{nlane*ng*16-1}:0] negs, output [{nsig*nlane*ow-1}:0] y);"]
    for l in range(nlane):
        for pr in range(nsig // 2):
            ka, kb = 2 * pr * nlane + l, (2 * pr + 1) * nlane + l
            L.append(f"  glane_{fam}_{sname}_w{W} u_l{l}_p{pr} (.clk(clk), .a(x[{(ka+1)*T*8-1}:{ka*T*8}]), .b(x[{(kb+1)*T*8-1}:{kb*T*8}]),"
                     f" .c(c[{(l+1)*T*W-1}:{l*T*W}]), .negs(negs[{(l+1)*ng*16-1}:{l*ng*16}]), .ya(y[{(ka+1)*ow-1}:{ka*ow}]), .yb(y[{(kb+1)*ow-1}:{kb*ow}]));")
    L.append("endmodule\n")
    return '\n'.join(L)

def base_lane(fam, W):
    # baseline: one signal per DSP, a 14-DSP cascade, W-bit coefficients (what SKA-Low / CASPER build)
    aw = 27 if fam == 'e2' else 25
    L = [f"module gbase_{fam}_w{W} (input clk, input [{T*8-1}:0] x, input [{T*W-1}:0] c, output [47:0] y);",
         f"  wire [47:0] pc [0:{T}];"]
    for t in range(T):
        A = f"{{{{22{{x[{t*8+7}]}}}}, x[{t*8+7}:{t*8}]}}"
        B = f"{{{{{18-W}{{c[{t*W+W-1}]}}}}, c[{t*W+W-1}:{t*W}]}}" if W < 18 else f"c[{t*W+W-1}:{t*W}]"
        L.append(dsp(fam, f"u{t}", A, B, "48'd0", "48'd0" if t == 0 else f"pc[{t}]", 'C' if t == 0 else 'PCIN',
                     "y" if t == T - 1 else "", "" if t == T - 1 else f"pc[{t+1}]"))
    L.append("endmodule\n"); return '\n'.join(L)
def base_slice(fam, W, nsig=4, nlane=4):
    L = [f"module gbslice_{fam}_w{W} (input clk, input [{nsig*nlane*T*8-1}:0] x, input [{nlane*T*W-1}:0] c, output [{nsig*nlane*48-1}:0] y);"]
    for l in range(nlane):
        for sg in range(nsig):
            k = sg * nlane + l
            L.append(f"  gbase_{fam}_w{W} u_l{l}_s{sg} (.clk(clk), .x(x[{(k+1)*T*8-1}:{k*T*8}]), .c(c[{(l+1)*T*W-1}:{l*T*W}]), .y(y[{(k+1)*48-1}:{k*48}]));")
    L.append("endmodule\n"); return '\n'.join(L)

if __name__ == '__main__' and sys.argv[1] == 'base':
    fam, W = sys.argv[2], int(sys.argv[3]); d = os.path.dirname(os.path.abspath(__file__))
    open(f"{d}/gbase_{fam}_w{W}.v", 'w').write(base_lane(fam, W)); open(f"{d}/gbslice_{fam}_w{W}.v", 'w').write(base_slice(fam, W))
    print("written base", fam, W); sys.exit(0)
if __name__ == '__main__':
    fam, sname, W = sys.argv[1], sys.argv[2], int(sys.argv[3])
    d = os.path.dirname(os.path.abspath(__file__))
    open(f"{d}/glane_{fam}_{sname}_w{W}.v", 'w').write(lane(fam, sname, W))
    open(f"{d}/glane_{fam}_{sname}_w{W}_broken.v", 'w').write(lane(fam, sname, W, True))
    open(f"{d}/gslice_{fam}_{sname}_w{W}.v", 'w').write(slice_mod(fam, sname, W))
    print(f"written glane_{fam}_{sname}_w{W}(.v, _broken.v), gslice; output width {out_width(fam, schemes()[sname], W)}")
