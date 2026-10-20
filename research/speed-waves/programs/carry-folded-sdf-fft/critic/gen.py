# gen.py (as-5, wave 14): one generator, two streaming radix-2^2 DIF FFT architectures that share everything
# except the butterfly stage:
#   sdf : R2^2SDF (single-path delay feedback; the vendor-core architecture): per stage a feedback delay of
#         h = 2^(n-s) words, 2 complex adders (sum, difference), 2 complex muxes.
#   sc  : serial commutator (Garrido et al. 2016 style): per stage an elementary bit exchange (time bits 0 and
#         n-s; delay 2^(n-s) - 1 words, 2 complex muxes) and a serial butterfly with ONE complex adder/subtractor
#         used every clock (sum on even, difference on odd slots) and 3 complex registers.
# Common: N = 2^n (n even), 16-bit complex in, one sample per clock, continuous frames, natural-order input.
# Every butterfly halves with ap-6's V6 rounding: S = (a>>1)+(b>>1)+b0, D = (a>>1)-(b>>1) (b = later sample).
# Trivial -j before every even stage is folded into the butterfly (swap + S/D role swap; exact, no negation).
# General rotator after stages 2,4,..,n-2: 18-bit twiddles (scale 2^16, so 1.0 exact), 4 multiplies,
# round half up (>>16), saturate to 16 bits; twiddle ROM addressed by the slot counter bits it needs.
# Output order: sdf = bit-reversed positions; sc = a fixed permutation (both unordered; the corner turn
# of a radar chain absorbs either by its write address).
import sys, argparse, json
import numpy as np

W = 16; TW = 18; TSH = 16


def s16(v):
    v = int(v) & 0xFFFF
    return v - 0x10000 if v & 0x8000 else v


# ------------------------------------------------------------------ bit-true model (positions, in place)
def S(a, b):
    return (a >> 1) + (b >> 1) + (b & 1)


def D(a, b):
    return (a >> 1) - (b >> 1)


def twiddle(e, Np):
    ang = -2 * np.pi * e / Np
    return int(np.round(np.cos(ang) * (1 << TSH))), int(np.round(np.sin(ang) * (1 << TSH)))


def rot_exp(pos, p, n):
    """exponent and modulus of the twiddle after stage 2p at position pos (p = 1..n/2-1)"""
    Np = 1 << (n - 2 * (p - 1))
    r = pos % (Np >> 2)
    b1 = (pos >> (n - 2 * p + 1)) & 1
    b2 = (pos >> (n - 2 * p)) & 1
    return (r * (b1 + 2 * b2)) % Np, Np


def rotate(zr, zi, wr, wi):
    pr = zr * wr - zi * wi
    pi = zr * wi + zi * wr
    rr = (pr + (1 << (TSH - 1))) >> TSH
    ri = (pi + (1 << (TSH - 1))) >> TSH
    return max(-32768, min(32767, rr)), max(-32768, min(32767, ri))


def model(xr, xi, n, broken=False):
    """xr, xi: int arrays of length N (natural order). Returns arrays by POSITION (pos i = X[bitrev(i)]/N)."""
    N = 1 << n
    ar = [int(v) for v in xr]; ai = [int(v) for v in xi]
    for s in range(1, n + 1):
        h = 1 << (n - s)
        for i in range(N):
            if i & h:
                continue
            j = i + h
            pr, pi_ = ar[i], ai[i]
            qr, qi = ar[j], ai[j]
            if s % 2 == 0 and (j >> (n - s + 1)) & 1:
                # a +/- (-j) b
                ar[i], ai[i], ar[j], ai[j] = S(pr, qi), D(pi_, qr), D(pr, qi), S(pi_, qr)
            else:
                ar[i], ai[i], ar[j], ai[j] = S(pr, qr), S(pi_, qi), D(pr, qr), D(pi_, qi)
        if s % 2 == 0 and s < n:
            p = s // 2
            for i in range(N):
                e, Np = rot_exp(i, p, n)
                wr, wi = twiddle(e, Np)
                if broken and p == 1 and i == (N >> 2) + 3:
                    wr += 4096
                ar[i], ai[i] = rotate(ar[i], ai[i], wr, wi)
    return ar, ai


# ------------------------------------------------------------------ Verilog emission
class G:
    def __init__(s, n, arch, broken=False):
        s.n = n; s.N = 1 << n; s.arch = arch; s.broken = broken
        s.L = []; s.k = 0; s.mods = []

    def w(s, t):
        s.L.append(t)

    def uid(s, p):
        s.k += 1
        return '%s%d' % (p, s.k)

    def pos(s, lat):
        """register holding the slot (0..N-1) of the sample at a unit input whose stream lags the top input by lat"""
        nm = s.uid('pos')
        s.w('  reg [%d:0] %s; always @(posedge clk) %s <= gcnt - %d\'d%d;' % (s.n - 1, nm, nm, s.n, (lat - 1) % s.N))
        return nm

    def delay(s, din_r, din_i, M):
        """complex delay of exactly M clocks (M>=1); returns (out_r, out_i)"""
        nm = s.uid('dl')
        if M == 1:
            s.w('  reg signed [15:0] %s_r, %s_i; always @(posedge clk) begin %s_r <= %s; %s_i <= %s; end' % (nm, nm, nm, din_r, nm, din_i))
            s.w('`ifndef SYNTHESIS\n  initial begin %s_r = 0; %s_i = 0; end\n`endif' % (nm, nm))
            return nm + '_r', nm + '_i'
        aw = max(1, (M - 1).bit_length())
        s.w('  reg [31:0] %s_m [0:%d]; reg [%d:0] %s_wp; reg [31:0] %s_q;' % (nm, M - 1, aw - 1, nm, nm))
        s.w('  wire [%d:0] %s_rp = (%s_wp == %d\'d%d) ? %d\'d0 : %s_wp + %d\'d1;' % (aw - 1, nm, nm, aw, M - 1, aw, nm, aw))
        s.w('  always @(posedge clk) begin %s_wp <= rst ? %d\'d0 : %s_rp; %s_m[%s_wp] <= {%s, %s}; %s_q <= %s_m[%s_rp]; end'
            % (nm, aw, nm, nm, nm, din_r, din_i, nm, nm, nm))
        s.w('  wire signed [15:0] %s_r = %s_q[31:16]; wire signed [15:0] %s_i = %s_q[15:0];' % (nm, nm, nm, nm))
        # simulation only: start the memory at 0 so that X does not leak through LUT-indexed cell models
        s.w('`ifndef SYNTHESIS\n  integer %s_k; initial begin for (%s_k = 0; %s_k < %d; %s_k = %s_k + 1) %s_m[%s_k] = 0; %s_q = 0; end\n`endif' % ((nm,) * 3 + (M,) + (nm,) * 5))
        return nm + '_r', nm + '_i'

    def rotator(s, zr, zi, posreg, p, P):
        """after stage 2p; posreg = slot counter at rotator input; P[time bit] = position bit. latency 4."""
        n, N = s.n, s.N
        Np = 1 << (n - 2 * (p - 1))
        need = list(range(0, n - 2 * p + 2))  # position bits that the exponent depends on
        tbits = sorted([tb for tb in range(n) if P[tb] in need])
        nm = s.uid('rot')
        aw = len(tbits)
        addr = '{' + ','.join('%s[%d]' % (posreg, tb) for tb in reversed(tbits)) + '}'
        s.w('  wire [%d:0] %s_a = %s;' % (aw - 1, nm, addr))
        s.w('  reg signed [17:0] %s_wr, %s_wi; reg signed [15:0] %s_zr, %s_zi;' % (nm, nm, nm, nm))
        s.w('  always @(posedge clk) begin %s_zr <= %s; %s_zi <= %s; end' % (nm, zr, nm, zi))
        s.w('  always @(posedge clk) case (%s_a)' % nm)
        for a in range(1 << aw):
            posv = 0
            for j, tb in enumerate(tbits):
                if (a >> j) & 1:
                    posv |= 1 << P[tb]
            e, _ = rot_exp(posv, p, n)
            wr, wi = twiddle(e, Np)
            # broken copy: the same entry the model's broken flag perturbs (p=1, position N/4+3)
            if s.broken and p == 1 and posv == (N >> 2) + 3:
                wr += 4096
            s.w('    %d\'d%d: begin %s_wr <= %d; %s_wi <= %d; end' % (aw, a, nm, wr, nm, wi))
        s.w('  endcase')
        s.w('  reg signed [33:0] %s_p0, %s_p1, %s_p2, %s_p3;' % (nm, nm, nm, nm))
        s.w('  always @(posedge clk) begin %s_p0 <= %s_zr * %s_wr; %s_p1 <= %s_zi * %s_wi; %s_p2 <= %s_zr * %s_wi; %s_p3 <= %s_zi * %s_wr; end'
            % ((nm,) * 12))
        s.w('  reg signed [18:0] %s_sr, %s_si;' % (nm, nm))
        s.w('  wire signed [34:0] %s_ur = %s_p0 - %s_p1 + 35\'sd%d; wire signed [34:0] %s_ui = %s_p2 + %s_p3 + 35\'sd%d;'
            % (nm, nm, nm, 1 << (TSH - 1), nm, nm, nm, 1 << (TSH - 1)))
        s.w('  always @(posedge clk) begin %s_sr <= %s_ur[34:16]; %s_si <= %s_ui[34:16]; end' % (nm, nm, nm, nm))
        for c in 'ri':
            s.w('  reg signed [15:0] %s_o%s; always @(posedge clk) %s_o%s <= (%s_s%s > 19\'sd32767) ? 16\'sd32767 : (%s_s%s < -19\'sd32768) ? -16\'sd32768 : %s_s%s[15:0];'
                % (nm, c, nm, c, nm, c, nm, c, nm, c))
        return nm + '_or', nm + '_oi', 4

    def sd(s, nm, a, b, sub):
        """V6 add/sub: sub ? (a>>1)-(b>>1) : (a>>1)+(b>>1)+b0 ; sub may be a constant '0'/'1' or a signal"""
        s.w('  wire signed [15:0] %s_ah = %s >>> 1; wire signed [15:0] %s_bh = %s >>> 1; wire %s_b0 = %s[0];' % (nm, a, nm, b, nm, b))
        if sub == '0':
            s.w('  wire signed [15:0] %s = %s_ah + %s_bh + {15\'d0, %s_b0};' % (nm, nm, nm, nm))
        elif sub == '1':
            s.w('  wire signed [15:0] %s = %s_ah - %s_bh;' % (nm, nm, nm))
        else:
            s.w('  wire signed [15:0] %s = %s_ah + (%s ? ~%s_bh : %s_bh) + {15\'d0, (%s | %s_b0)};' % (nm, nm, sub, nm, nm, sub, nm))
        return nm

    def sdf_stage(s, xr, xi, st, lat, P):
        n = s.n; h = 1 << (n - st)
        nm = s.uid('sdf')
        pr = s.pos(lat)
        ph = '%s[%d]' % (pr, n - st)
        s.w('  wire signed [15:0] %s_dr, %s_di;' % (nm, nm))
        mr, mi = s.delay(nm + '_dr', nm + '_di', h)
        if st % 2 == 0:
            jf = '%s_jf' % nm
            s.w('  wire %s = %s[%d];' % (jf, pr, n - st + 1))
            s.w('  wire signed [15:0] %s_br = %s ? %s : %s; wire signed [15:0] %s_bi = %s ? %s : %s;' % (nm, jf, xi, xr, nm, jf, xr, xi))
            Sr = s.sd(nm + '_Sr', mr, nm + '_br', '0'); Dr = s.sd(nm + '_Dr', mr, nm + '_br', '1')
            Pi = s.sd(nm + '_Pi', mi, nm + '_bi', '0'); Qi = s.sd(nm + '_Qi', mi, nm + '_bi', '1')
            s.w('  wire signed [15:0] %s_Si = %s ? %s : %s; wire signed [15:0] %s_Di = %s ? %s : %s;' % (nm, jf, Qi, Pi, nm, jf, Pi, Qi))
            Si, Di = nm + '_Si', nm + '_Di'
        else:
            Sr = s.sd(nm + '_Sr', mr, xr, '0'); Dr = s.sd(nm + '_Dr', mr, xr, '1')
            Si = s.sd(nm + '_Si', mi, xi, '0'); Di = s.sd(nm + '_Di', mi, xi, '1')
        s.w('  assign %s_dr = %s ? %s : %s; assign %s_di = %s ? %s : %s;' % (nm, ph, Dr, xr, nm, ph, Di, xi))
        s.w('  reg signed [15:0] %s_or, %s_oi; always @(posedge clk) begin %s_or <= %s ? %s : %s; %s_oi <= %s ? %s : %s; end'
            % (nm, nm, nm, ph, Sr, mr, nm, ph, Si, mi))
        return nm + '_or', nm + '_oi', lat + h + 1, P

    def sc_stage(s, xr, xi, st, lat, P):
        n = s.n
        P = list(P)
        k = P.index(n - st)  # time bit that holds the position bit this stage pairs on
        if k > 0:
            # elementary bit exchange of time bits 0 and k; delay L = 2^k - 1
            Lk = (1 << k) - 1
            nm = s.uid('bx')
            pe = s.pos(lat)
            s.w('  wire %s_byp = ~%s[0] & %s[%d];' % (nm, pe, pe, k))
            s.w('  wire signed [15:0] %s_ir, %s_ii;' % (nm, nm))
            dr, di = s.delay(nm + '_ir', nm + '_ii', Lk)
            s.w('  assign %s_ir = %s_byp ? %s : %s; assign %s_ii = %s_byp ? %s : %s;' % (nm, nm, dr, xr, nm, nm, di, xi))
            s.w('  wire signed [15:0] %s_xr = %s_byp ? %s : %s; wire signed [15:0] %s_xi = %s_byp ? %s : %s;' % (nm, nm, xr, dr, nm, nm, xi, di))
            xr, xi = nm + '_xr', nm + '_xi'
            lat += Lk
            P[0], P[k] = P[k], P[0]
        assert P[0] == n - st
        nm = s.uid('bf')
        pb = s.pos(lat)
        odd = '%s[0]' % pb
        s.w('  reg signed [15:0] %s_ar, %s_ai, %s_br, %s_bi, %s_cr, %s_ci; reg %s_jf;' % ((nm,) * 7))
        if st % 2 == 0:
            jt = P.index(n - st + 1)
            jsrc = '%s[%d]' % (pb, jt)
            s.w('  always @(posedge clk) if (%s) begin %s_br <= %s ? %s : %s; %s_bi <= %s ? %s : %s; %s_jf <= %s; end'
                % (odd, nm, jsrc, xi, xr, nm, jsrc, xr, xi, nm, jsrc))
            subi = '(%s ^ %s_jf)' % (odd, nm)
        else:
            s.w('  always @(posedge clk) if (%s) begin %s_br <= %s; %s_bi <= %s; end' % (odd, nm, xr, nm, xi))
            subi = odd
        s.w('  always @(posedge clk) if (!%s) begin %s_ar <= %s; %s_ai <= %s; end' % (odd, nm, xr, nm, xi))
        s.w('  always @(posedge clk) if (%s) begin %s_cr <= %s_ar; %s_ci <= %s_ai; end' % (odd, nm, nm, nm, nm))
        rr = s.sd(nm + '_yr', nm + '_cr', nm + '_br', odd)
        ri = s.sd(nm + '_yi', nm + '_ci', nm + '_bi', subi)
        s.w('  reg signed [15:0] %s_or, %s_oi; always @(posedge clk) begin %s_or <= %s; %s_oi <= %s; end' % (nm, nm, nm, rr, nm, ri))
        return nm + '_or', nm + '_oi', lat + 3, P

    # ---------------- fused SDF (ECP5): the stage's muxes absorbed into the CCU2C carry cells
    # ECP5 CCU2C half: S = P ^ cin, cout = P ? cin : G, P = LUT4(A,B,C,D) (INIT), G = LUT2(A,B) = INIT[3:0].
    CELLS = {
        # name: (P(A,B,C,D), G(A,B) or None if D is data (then G = P(A,B,0,0)))
        'OUT':  (lambda A, B, C, D: A ^ (B & C), lambda A, B: A),            # A=m  B=xh  C=ph  D=1
        'MEM':  (lambda A, B, C, D: (C ^ (1 - A)) if B else A, lambda A, B: B & (1 - A)),  # A=xh B=ph C=m D=1
        'OUTJ': (lambda A, B, C, D: A ^ (C & (B ^ D)), None),               # A=m  B=bh  C=ph  D=jf
        'MEMJ': (lambda A, B, C, D: (A ^ C ^ D) if B else C, None),         # A=m  B=ph  C=bh  D=njf
    }

    def cell_init(s, kind):
        P, G = s.CELLS[kind]
        v = 0
        for idx in range(16):
            A, B, C, D = idx & 1, (idx >> 1) & 1, (idx >> 2) & 1, (idx >> 3) & 1
            if G is None:
                bit = P(A, B, C, D)
            else:
                bit = P(A, B, C, 1) if D else G(A, B)
            v |= bit << idx
        if G is None:
            assert all(((v >> (a | (b << 1))) & 1) == P(a, b, 0, 0) for a in (0, 1) for b in (0, 1))
        return v

    def fused(s, nm, kind, A, B, C, D, cin):
        """16-bit fused adder: A,B,C,D lists of 16 bit expressions; returns wire name [15:0]"""
        init = s.cell_init(kind)
        s.w('  wire [15:0] %s; wire [19:0] %s_c;' % (nm, nm))
        halves = [('%s' % cin, '1\'b0', '1\'b0', '1\'b1', 0x00AA, None)] + [(A[i], B[i], C[i], D[i], init, i) for i in range(16)] + [('1\'b0', '1\'b0', '1\'b0', '1\'b1', 0x0000, None)]
        for j in range(0, 18, 2):
            h0, h1 = halves[j], halves[j + 1]
            s0 = '%s[%d]' % (nm, h0[5]) if h0[5] is not None else ''
            s1 = '%s[%d]' % (nm, h1[5]) if h1[5] is not None else ''
            cinj = '1\'b0' if j == 0 else '%s_c[%d]' % (nm, j)
            s.w('  CCU2C #(.INIT0(16\'h%04x), .INIT1(16\'h%04x), .INJECT1_0("NO"), .INJECT1_1("NO")) %s_k%d (.CIN(%s), .A0(%s), .B0(%s), .C0(%s), .D0(%s), .A1(%s), .B1(%s), .C1(%s), .D1(%s), .S0(%s), .S1(%s), .COUT(%s_c[%d]));'
                % (h0[4], h1[4], nm, j, cinj, h0[0], h0[1], h0[2], h0[3], h1[0], h1[1], h1[2], h1[3], s0, s1, nm, j + 2))
        return nm

    CELLS_N = {
        'OUT':   (lambda A, B, C, D: A ^ (B & C), lambda A, B, C: A),                 # A=m  B=xh  C=ph  D=1
        'MEM':   (lambda A, B, C, D: (C ^ (1 - A)) if B else A, lambda A, B, C: B & (1 - A)),  # A=xh B=ph C=m D=1
        'OUTJN': (lambda A, B, C, D: A ^ (D & (B ^ C)), None),                       # A=m  B=bh  C=jf  D=ph
        'MEMJN': (lambda A, B, C, D: (D ^ B ^ C) if A else B, None),                 # A=ph B=bh  C=njf D=m
    }

    def cell_init_n(s, kind):
        P, G = s.CELLS_N[kind]
        v = 0
        for idx in range(16):
            A, B, C, D = idx & 1, (idx >> 1) & 1, (idx >> 2) & 1, (idx >> 3) & 1
            bit = (P(A, B, C, 1) if D else G(A, B, C)) if G is not None else P(A, B, C, D)
            v |= bit << idx
        return v

    def fused_n(s, nm, kind, A, B, C, D, cin):
        init = s.cell_init_n(kind)
        s.w('  wire [15:0] %s; wire [19:0] %s_c;' % (nm, nm))
        # feed-in cell as Yosys' arith_map_nexus.v: half 0 passes (P=1), half 1 has P=0 and G=A1=cin
        s.w('  CCU2 #(.INIT0("0xFFFF"), .INIT1("0x00AA"), .INJECT("NO")) %s_fi (.A0(1\'b1), .B0(1\'b1), .C0(1\'b1), .D0(1\'b1), .A1(%s), .B1(1\'b1), .C1(1\'b1), .D1(1\'b1), .COUT(%s_c[0]));' % (nm, cin, nm))
        for j in range(0, 16, 2):
            s.w('  CCU2 #(.INIT0("0x%04X"), .INIT1("0x%04X"), .INJECT("NO")) %s_k%d (.CIN(%s_c[%d]), .A0(%s), .B0(%s), .C0(%s), .D0(%s), .A1(%s), .B1(%s), .C1(%s), .D1(%s), .S0(%s[%d]), .S1(%s[%d]), .COUT(%s_c[%d]));'
                % (init, init, nm, j, nm, j, A[j], B[j], C[j], D[j], A[j + 1], B[j + 1], C[j + 1], D[j + 1], nm, j, nm, j + 1, nm, j + 2))
        return nm

    def sdfx_stage(s, xr, xi, st, lat, P):
        n = s.n; h = 1 << (n - st)
        nm = s.uid('sdx')
        pr = s.pos(lat)
        ph = '%s_ph' % nm
        s.w('  wire %s = %s[%d];' % (ph, pr, n - st))
        s.w('  wire signed [15:0] %s_dr, %s_di;' % (nm, nm))
        mr, mi = s.delay(nm + '_dr', nm + '_di', h)   # memory holds a>>1 (phase 0) or D (phase 1)
        if st % 2 == 0:
            s.w('  wire %s_jf = %s & %s[%d]; wire %s_nj = ~%s_jf;' % (nm, ph, pr, n - st + 1, nm, nm))
            s.w('  wire signed [15:0] %s_br = %s_jf ? %s : %s; wire signed [15:0] %s_bi = %s_jf ? %s : %s;' % (nm, nm, xi, xr, nm, nm, xr, xi))
            br, bi = nm + '_br', nm + '_bi'
        else:
            br, bi = xr, xi
        def hb(x, i):  # bit i of x>>>1
            return '%s[%d]' % (x, min(i + 1, 15))
        one = "1'b1"
        # real lane (and imag lane of odd stages): OUT / MEM
        lanes = [('r', mr, br, False), ('i', mi, bi, st % 2 == 0)]
        outs = {}
        for c, m, b, jl in lanes:
            if s.arch == 'sdfn':
                M = ['%s[%d]' % (m, i) for i in range(16)]; Bh = [hb(b, i) for i in range(16)]
                if not jl:
                    o = s.fused_n(nm + '_fo' + c, 'OUT', M, Bh, [ph] * 16, [one] * 16, '(%s[0] & %s)' % (b, ph))
                    d = s.fused_n(nm + '_fw' + c, 'MEM', Bh, [ph] * 16, M, [one] * 16, ph)
                else:
                    o = s.fused_n(nm + '_fo' + c, 'OUTJN', M, Bh, [nm + '_jf'] * 16, [ph] * 16, '(%s & (%s_jf | %s[0]))' % (ph, nm, b))
                    d = s.fused_n(nm + '_fw' + c, 'MEMJN', [ph] * 16, Bh, [nm + '_nj'] * 16, M, '(%s & (%s_nj | %s[0]))' % (ph, nm, b))
                outs[c] = (o, d)
                continue
            if not jl:
                o = s.fused(nm + '_fo' + c, 'OUT', ['%s[%d]' % (m, i) for i in range(16)], [hb(b, i) for i in range(16)], [ph] * 16, [one] * 16, '(%s[0] & %s)' % (b, ph))
                d = s.fused(nm + '_fw' + c, 'MEM', [hb(b, i) for i in range(16)], [ph] * 16, ['%s[%d]' % (m, i) for i in range(16)], [one] * 16, ph)
            else:
                o = s.fused(nm + '_fo' + c, 'OUTJ', ['%s[%d]' % (m, i) for i in range(16)], [hb(b, i) for i in range(16)], [ph] * 16, [nm + '_jf'] * 16, '(%s & (%s_jf | %s[0]))' % (ph, nm, b))
                d = s.fused(nm + '_fw' + c, 'MEMJ', ['%s[%d]' % (m, i) for i in range(16)], [ph] * 16, [hb(b, i) for i in range(16)], [nm + '_nj'] * 16, '(%s & (%s_nj | %s[0]))' % (ph, nm, b))
            outs[c] = (o, d)
        s.w('  assign %s_dr = %s; assign %s_di = %s;' % (nm, outs['r'][1], nm, outs['i'][1]))
        s.w('  reg signed [15:0] %s_or, %s_oi; always @(posedge clk) begin %s_or <= %s; %s_oi <= %s; end' % (nm, nm, nm, outs['r'][0], nm, outs['i'][0]))
        return nm + '_or', nm + '_oi', lat + h + 1, P

    # ---------------- au-2: the SAME gated formulation as sdfx/sdfn, but written as plain generic RTL
    # (memory holds a>>1 in phase 0); left to Yosys/ABC to map. The obvious thing a practitioner tries.
    def sdfg_stage(s, xr, xi, st, lat, P):
        n = s.n; h = 1 << (n - st)
        nm = s.uid('sdg')
        pr = s.pos(lat)
        ph = '%s_ph' % nm
        s.w('  wire %s = %s[%d];' % (ph, pr, n - st))
        s.w('  wire signed [15:0] %s_dr, %s_di;' % (nm, nm))
        mr, mi = s.delay(nm + '_dr', nm + '_di', h)
        if st % 2 == 0:
            s.w('  wire %s_jf = %s & %s[%d];' % (nm, ph, pr, n - st + 1))
            s.w('  wire signed [15:0] %s_br = %s_jf ? %s : %s; wire signed [15:0] %s_bi = %s_jf ? %s : %s;' % (nm, nm, xi, xr, nm, nm, xr, xi))
            br, bi = nm + '_br', nm + '_bi'
        else:
            br, bi = xr, xi
        outs = {}
        for c, m, b, jl in [('r', mr, br, False), ('i', mi, bi, st % 2 == 0)]:
            xh = '%s_h%s' % (nm, c); x0 = '%s[0]' % b
            s.w('  wire signed [15:0] %s = %s >>> 1;' % (xh, b))
            if not jl:
                s.w('  wire signed [15:0] %s_y%s = %s + (%s ? %s : 16\'sd0) + {15\'d0, %s & %s};' % (nm, c, m, ph, xh, ph, x0))
                s.w('  wire signed [15:0] %s_w%s = %s ? (%s - %s) : %s;' % (nm, c, ph, m, xh, xh))
            else:
                jf = nm + '_jf'
                s.w('  wire signed [15:0] %s_y%s = %s + (%s ? (%s ? ~%s : %s) : 16\'sd0) + {15\'d0, %s & (%s | %s)};' % (nm, c, m, ph, jf, xh, xh, ph, jf, x0))
                s.w('  wire signed [15:0] %s_w%s = %s ? (%s + (%s ? %s : ~%s) + {15\'d0, (%s ? %s : 1\'b1)}) : %s;' % (nm, c, ph, m, jf, xh, xh, jf, x0, xh))
            outs[c] = ('%s_y%s' % (nm, c), '%s_w%s' % (nm, c))
        s.w('  assign %s_dr = %s; assign %s_di = %s;' % (nm, outs['r'][1], nm, outs['i'][1]))
        s.w('  reg signed [15:0] %s_or, %s_oi; always @(posedge clk) begin %s_or <= %s; %s_oi <= %s; end' % (nm, nm, nm, outs['r'][0], nm, outs['i'][0]))
        return nm + '_or', nm + '_oi', lat + h + 1, P

    def build(s, name):
        n = s.n
        s.w('module %s (input clk, input rst, input signed [15:0] in_re, input signed [15:0] in_im,' % name)
        s.w('  output signed [15:0] out_re, output signed [15:0] out_im);')
        s.w('  reg [%d:0] gcnt; always @(posedge clk) gcnt <= rst ? %d\'d0 : gcnt + %d\'d1;' % (n - 1, n, n))
        s.w('  reg signed [15:0] x0r, x0i; always @(posedge clk) begin x0r <= in_re; x0i <= in_im; end')
        # sample entering at the cycle with gcnt = c is at x0 one cycle later
        xr, xi, lat = 'x0r', 'x0i', 1
        P = list(range(n))  # P[time bit] = position bit
        for st in range(1, n + 1):
            if s.arch == 'sdf':
                xr, xi, lat, P = s.sdf_stage(xr, xi, st, lat, P)
            elif s.arch in ('sdfx', 'sdfn'):
                xr, xi, lat, P = s.sdfx_stage(xr, xi, st, lat, P)
            elif s.arch == 'sdfg':
                xr, xi, lat, P = s.sdfg_stage(xr, xi, st, lat, P)
            else:
                xr, xi, lat, P = s.sc_stage(xr, xi, st, lat, P)
            if st % 2 == 0 and st < n:
                pr = s.pos(lat)
                xr, xi, dl = s.rotator(xr, xi, pr, st // 2, P)
                lat += dl
        s.w('  assign out_re = %s; assign out_im = %s;' % (xr, xi))
        s.w('endmodule')
        s.lat = lat; s.P = P
        return '\n'.join(s.L) + '\n'


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=8)
    ap.add_argument('--arch', default='sdf')
    ap.add_argument('--broken', action='store_true')
    ap.add_argument('-o', required=True)
    a = ap.parse_args()
    g = G(a.n, a.arch, a.broken)
    name = 'fft_%s_%d%s' % (a.arch, 1 << a.n, '_broken' if a.broken else '')
    txt = g.build(name)
    open(a.o, 'w').write(txt)
    json.dump({'name': name, 'lat': g.lat, 'P': g.P, 'n': a.n, 'arch': a.arch}, open(a.o + '.json', 'w'))
    print(name, 'latency', g.lat, 'P', g.P)
