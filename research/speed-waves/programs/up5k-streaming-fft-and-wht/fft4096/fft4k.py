#!/usr/bin/env python3
"""fft4k.py (ap-6, wave 13) -- a 4,096-point streaming FFT for the iCE40-UP5K: bit-true numpy model, Verilog
generator, and simulation check.

Architecture: radix-2^4 single-path delay-feedback (SDF), decimation in frequency, 12 stages, one complex 16-bit
sample per clock, output in bit-reversed order (slot u holds X[bitrev12(u)] / N).
  stage s (1..12): butterfly with a feedback delay of D = 2^(12-s) complex words; every stage halves (scale 1/N).
  rotators: after stages 1,3,5,7,9,11: trivial -j; after 2,6,10: W16 constant rotator (shift-add, CSD);
            after 4 and 8: general rotator (4 SB_MAC16 multipliers each, twiddle octant ROM).
  memories: D=2048 in the 4 SPRAMs (two 32-bit banks, ping-pong), D=8..1024 in block RAM, D<=4 in flip-flops.

Rounding modes of the butterfly (env/arg ROUND):
  V0  halve both inputs first (floor), then add/sub (am-5's WHT arithmetic)
  V7  sum: (a+b)>>1 on odd stages (floor), (a+b+1)>>1 on even stages (round half up), i.e. 17-bit adder with a
      constant carry-in; difference: halve first, a' - b' (unbiased).  [chosen: one carry cell more per sum]
  V4  both sum and difference rounded after the add, alternating floor / round half up (reference only: the
      difference's round-half-up needs a carry-in of 2, not built)
Everything wraps to 16 bits (two's complement), in the model exactly as in the hardware.

Usage:
  fft4k.py model [--frames F]                 precision report (SQNR, worst spur) for V0/V7/V4
  fft4k.py gen -o fft.v [--round V7] [--mem ebr|spram] [--rot mac|lc] [--break twiddle|round]
  fft4k.py sim fft.v --samples 110000 --work DIR [--round V7]
"""
import argparse, math, os, subprocess, sys
import numpy as np

LOGN = 12
N = 1 << LOGN
G = int(os.environ.get('G', '2'))
W16CHAIN = os.environ.get('W16CHAIN', '1') == '1'
SWAP = os.environ.get('SWAP', '1') == '1'
W16PIPE = os.environ.get('W16PIPE', '1') == '1'   # pipelined W16 products (ap-6 v4)   # swapped butterfly (output a-b, store a+b): every k bit is inverted in the slot label  # guard bits in the W16 constant products (v1: 4)

# CSD constants (value = sum s * 2^-k), searched in ap-6/fft (fft4k.py csd): errors -85, -95, -79 dB
CSD = {
    'c1': [(1, 0), (-1, 4), (-1, 6), (1, 9)],          # cos(pi/8)  0.9238281
    's1': [(1, 1), (-1, 3), (1, 7), (-1, 13)],         # sin(pi/8)  0.3826904
    'c2': [(1, 0), (-1, 2), (-1, 4), (1, 6), (1, 8)],  # cos(pi/4)  0.7070312
}

def wrap16(v):
    return ((v + 32768) & 0xFFFF) - 32768

def bitrev(x, bits):
    r = 0
    for i in range(bits):
        r |= ((x >> i) & 1) << (bits - 1 - i)
    return r

BR12 = np.array([bitrev(u ^ (0xFFF if SWAP else 0), 12) for u in range(N)])   # slot u holds X[BR12[u]]

# ---------- twiddles ----------
def octant_rom():
    """T[r] for r = 0..512: (C, S) = (round(32768 cos), round(32768 sin)) of 2*pi*r/4096, C clipped to 32767."""
    T = []
    for r in range(513):
        c = min(int(round(32768 * math.cos(2 * math.pi * r / N))), 32767)
        s = min(int(round(32768 * math.sin(2 * math.pi * r / N))), 32767)
        T.append((c, s))
    return T
TROM = octant_rom()

def twiddle_cs(e):
    """quadrant q, (C, S) for W_N^e = (-j)^q * (C - jS)/32768."""
    q = (e >> 10) & 3
    r = e & 1023
    if r <= 512:
        c, s = TROM[r]
    else:
        s, c = TROM[1024 - r]
    return q, c, s

def rot1_exp(u):  # after stage 4
    K = bitrev((u >> 8) & 15, 4)
    return K * (u & 255)

def rot2_exp(u):  # after stage 8
    K = bitrev((u >> 4) & 15, 4)
    return 16 * (K * (u & 15))

def w16_exp(u, s):  # W16 exponent after stage s in {2, 6, 10}
    p = 12 - s  # bits p+1, p are k bits; p-1, p-2 are input bits
    K = ((u >> (p + 1)) & 1) + 2 * ((u >> p) & 1)
    M = 2 * ((u >> (p - 1)) & 1) + ((u >> (p - 2)) & 1)
    return K * M

def trivial_on(u, s):  # -j after stage s in {1,3,5,7,9,11}
    p = 12 - s
    return ((u >> p) & 1) & ((u >> (p - 1)) & 1)

# ---------- bit-true model (frames x N arrays, slot order) ----------
def prerot(xr, xi, q):
    """multiply by (-j)^q, wrap: q=1: (xi, -xr); q=2: (-xr, -xi); q=3: (-xi, xr)"""
    r0 = np.where(q == 0, xr, np.where(q == 1, xi, np.where(q == 2, -xr, -xi)))
    i0 = np.where(q == 0, xi, np.where(q == 1, -xr, np.where(q == 2, -xi, xr)))
    return wrap16(r0), wrap16(i0)

def csd_prod(x, name):
    xs = x << G
    p = np.zeros_like(x)
    for s, k in CSD[name]:
        p = p + s * (xs >> k)
    return p

def butterfly(a, b, s, mode):
    if mode == 'V0':
        return wrap16((a >> 1) + (b >> 1)), wrap16((a >> 1) - (b >> 1))
    if mode == 'V6':  # halve first, carry-in = LSB of the later input into the sum (stored halved: muxes merge into adders)
        return wrap16((a >> 1) + (b >> 1) + (b & 1)), wrap16((a >> 1) - (b >> 1))
    if mode == 'V7':
        r = 0 if s % 2 == 1 else 1
        return wrap16((a + b + r) >> 1), wrap16((a >> 1) - (b >> 1))
    if mode == 'V4':
        r = 0 if s % 2 == 1 else 1
        return wrap16((a + b + r) >> 1), wrap16((a - b + r) >> 1)
    raise SystemExit('mode ' + mode)

U = np.arange(N)
def km(s):  # slot bits holding the k bits after stage s (inverted when SWAP)
    return (((1 << s) - 1) << (12 - s)) if SWAP else 0
E1 = np.array([rot1_exp(u ^ km(4)) for u in U]); E2 = np.array([rot2_exp(u ^ km(8)) for u in U])
TW1 = np.array([twiddle_cs(e) for e in E1]); TW2 = np.array([twiddle_cs(e) for e in E2])
W16E = {s: np.array([w16_exp(u ^ km(s), s) for u in U]) for s in (2, 6, 10)}
TRIV = {s: np.array([trivial_on(u ^ km(s), s) for u in U]) for s in (1, 3, 5, 7, 9, 11)}

def model(xr, xi, mode='V7', brk=None):
    """xr, xi: int arrays (F, N) of 16-bit values in natural input order. Returns (yr, yi) in slot order."""
    xr = xr.astype(np.int64).copy(); xi = xi.astype(np.int64).copy()
    F = xr.shape[0]
    for s in range(1, LOGN + 1):
        D = N >> s
        ar = xr.reshape(F, -1, 2, D); ai = xi.reshape(F, -1, 2, D)
        sr, dr = butterfly(ar[:, :, 0, :], ar[:, :, 1, :], s, mode)
        si, di = butterfly(ai[:, :, 0, :], ai[:, :, 1, :], s, mode)
        if SWAP:
            xr = np.stack([dr, sr], axis=2).reshape(F, N); xi = np.stack([di, si], axis=2).reshape(F, N)
        else:
            xr = np.stack([sr, dr], axis=2).reshape(F, N); xi = np.stack([si, di], axis=2).reshape(F, N)
        if s in TRIV:
            on = TRIV[s].astype(bool)
            nr, ni = wrap16(xi), wrap16(-xr)
            xr = np.where(on, nr, xr); xi = np.where(on, ni, xi)
        elif s in W16E:
            e = W16E[s]
            q, r = e >> 2, e & 3
            ur, ui = prerot(xr, xi, np.broadcast_to(q, xr.shape))
            P = {nm: (csd_prod(ur, nm), csd_prod(ui, nm)) for nm in CSD}
            P0 = (ur << G, ui << G)
            Cn = {1: 'c1', 2: 'c2', 3: 's1'}; Sn = {1: 's1', 2: 'c2', 3: 'c1'}
            cr = np.select([r == 0, r == 1, r == 2, r == 3], [P0[0], P['c1'][0], P['c2'][0], P['s1'][0]])
            ci = np.select([r == 0, r == 1, r == 2, r == 3], [P0[1], P['c1'][1], P['c2'][1], P['s1'][1]])
            sr_ = np.select([r == 0, r == 1, r == 2, r == 3], [0 * ur, P['s1'][0], P['c2'][0], P['c1'][0]])
            si_ = np.select([r == 0, r == 1, r == 2, r == 3], [0 * ui, P['s1'][1], P['c2'][1], P['c1'][1]])
            half = 1 << (G - 1)
            xr = wrap16((cr + si_ + half) >> G)
            xi = wrap16((ci - sr_ + half) >> G)
        elif s in (4, 8):
            T = TW1 if s == 4 else TW2
            q, C, S = T[:, 0], T[:, 1], T[:, 2]
            if brk == 'twiddle' and s == 4:
                C = C.copy(); C[1234] ^= 1
            ur, ui = prerot(xr, xi, np.broadcast_to(q, xr.shape))
            re = ur * C + 16384 + ui * S
            im = ui * C + 16384 - ur * S
            xr = wrap16(re >> 15); xi = wrap16(im >> 15)
    return xr, xi

def reference(xr, xi):
    X = np.fft.fft(xr.astype(np.float64) + 1j * xi.astype(np.float64), axis=1) / N
    return X[:, BR12]

# ---------- precision report ----------
def precision(frames=16, modes=('V0', 'V6', 'V7', 'V4')):
    rng = np.random.default_rng(5)
    n = np.arange(N)
    cases = []
    # complex Gaussian noise at -12 dBFS per component
    cases.append(('noise -12 dBFS (complex)', np.round(rng.normal(0, 8192, (frames, N))), np.round(rng.normal(0, 8192, (frames, N)))))
    # real Gaussian noise -12 dBFS (vibration / audio)
    cases.append(('noise -12 dBFS (real)', np.round(rng.normal(0, 8192, (frames, N))), np.zeros((frames, N))))
    # real tone off-bin, -1 dBFS, plus a little noise
    f0 = 0.1234567
    tr = np.round(0.891 * 32767 * np.cos(2 * np.pi * f0 * n + rng.uniform(0, 6.28, (frames, 1))))
    cases.append(('real tone -1 dBFS, off-bin', tr, np.zeros((frames, N))))
    # real tone at a bin centre (bin 333), -1 dBFS: the spur figure
    cases.append(('real tone -1 dBFS, bin 333', np.round(0.891 * 32767 * np.cos(2 * np.pi * 333 * n / N + rng.uniform(0, 6.28, (frames, 1)))), np.zeros((frames, N))))
    # complex tone at a bin centre, -1 dBFS magnitude
    ph = rng.uniform(0, 6.28, (frames, 1))
    cases.append(('complex tone -1 dBFS, bin 333', np.round(0.891 * 32767 * np.cos(2 * np.pi * 333 * n / N + ph)),
                  np.round(0.891 * 32767 * np.sin(2 * np.pi * 333 * n / N + ph))))
    out = []
    for name, xr, xi in cases:
        xr = np.clip(xr, -32767, 32767).astype(np.int64); xi = np.clip(xi, -32767, 32767).astype(np.int64)
        ref = reference(xr, xi)
        for m in modes:
            yr, yi = model(xr, xi, m)
            y = yr + 1j * yi
            err = y - ref
            sqnr = 10 * np.log10(np.mean(np.abs(ref) ** 2) / np.mean(np.abs(err) ** 2))
            line = '%-32s %s  SQNR %6.2f dB  max|err| %.2f LSB  mean err (re) %+.3f' % (name, m, sqnr, np.abs(err).max(), err.real.mean())
            if 'tone' in name:
                # worst spur: largest output bin outside the tone's +-3 bins (in natural order), relative to the peak (dBc)
                Y = np.zeros_like(y); Y[:, BR12] = y
                P = np.abs(Y) ** 2
                pk = P.argmax(axis=1)
                mask = np.ones_like(P, dtype=bool)
                for f in range(P.shape[0]):
                    for d in range(-3, 4):
                        mask[f, (pk[f] + d) % N] = False
                        mask[f, (N - pk[f] + d) % N] = False  # image of a real tone
                spur = 10 * np.log10((np.where(mask, P, 0).max(axis=1) / P.max(axis=1)).max())
                Rf = np.zeros_like(ref); Rf[:, BR12] = ref
                Pr = np.abs(Rf) ** 2
                spur_ref = 10 * np.log10((np.where(mask, Pr, 0).max(axis=1) / Pr.max(axis=1)).max())
                line += '  worst spur %.1f dBc (input itself %.1f dBc)' % (spur, spur_ref)
            print(line); out.append(line); sys.stdout.flush()
    return out

# ---------- Verilog generator ----------
def gen(mode='V7', mem='spram', rot='mac', brk=None, name='fft4k'):
    L = []; w = L.append
    w('// generated by fft4k.py: N=4096 radix-2^4 SDF DIF, 16-bit, round=%s mem=%s rot=%s break=%s' % (mode, mem, rot, brk))
    w('module %s (input clk, input rst, input [15:0] xr, input [15:0] xi, output [15:0] yr, output [15:0] yi, output y_first);' % name)
    w('  reg [11:0] n; wire [11:0] n1 = n + 12\'d1;')
    w('  always @(posedge clk) n <= rst ? 12\'d0 : n1;')
    off = 0  # data at the current point carries slot (t - off) mod N
    cur = ('xr', 'xi')
    for s in range(1, LOGN + 1):
        k = LOGN - s; D = 1 << k
        ar, ai = cur
        w('  // ===== stage %d: D=%d' % (s, D))
        # phase counter: bit k of (t - off)
        w('  reg [%d:0] c%d; always @(posedge clk) c%d <= rst ? %d\'d%d : c%d + 1\'b1;' % (k, s, s, k + 1, (-off) % (2 * D), s))
        w('  wire ph%d = c%d[%d];' % (s, s, k))
        w('  wire signed [15:0] br%d = %s, bi%d = %s;' % (s, ar, s, ai))
        w('  wire [31:0] din%d; wire [31:0] dout%d;' % (s, s))
        w('  wire signed [15:0] dr%d = dout%d[31:16], di%d = dout%d[15:0];' % (s, s, s, s))
        # butterfly: a = delayed (dr), b = current (br)
        if mode == 'V6':
            for c in ('r', 'i'):
                a = 'd%s%d' % (c, s); b = 'b%s%d' % (c, s)
                w('  wire signed [15:0] h%s%d = {%s[15], %s[15:1]};' % (c, s, b, b))
                cin = '%s[0]' % b
                if brk == 'round' and s == 7 and c == 'r':
                    cin = "1'b0"
                w('  wire signed [15:0] su%s%d = %s + h%s%d + {15\'d0, %s};' % (c, s, a, c, s, cin))
                w('  wire signed [15:0] df%s%d = %s - h%s%d;' % (c, s, a, c, s))
            if SWAP:
                w('  assign din%d = ph%d ? {sur%d, sui%d} : {hr%d, hi%d};' % (s, s, s, s, s, s))
            else:
                w('  assign din%d = ph%d ? {dfr%d, dfi%d} : {hr%d, hi%d};' % (s, s, s, s, s, s))
        for c in (('r', 'i') if mode != 'V6' else ()):
            a = 'd%s%d' % (c, s); b = 'b%s%d' % (c, s)
            if mode == 'V0':
                w('  wire signed [15:0] su%s%d = {%s[15], %s[15:1]} + {%s[15], %s[15:1]};' % (c, s, a, a, b, b))
            elif mode in ('V7', 'V4'):
                r = 0 if s % 2 == 1 else 1
                if brk == 'round' and s == 7 and c == 'r':
                    r = 1 - r
                w('  wire signed [16:0] sw%s%d = {%s[15], %s} + {%s[15], %s} + 17\'d%d;' % (c, s, a, a, b, b, r))
                w('  wire signed [15:0] su%s%d = sw%s%d[16:1];' % (c, s, c, s))
            w('  wire signed [15:0] df%s%d = {%s[15], %s[15:1]} - {%s[15], %s[15:1]};' % (c, s, a, a, b, b))
        if mode != 'V6':
            w('  assign din%d = ph%d ? {dfr%d, dfi%d} : {br%d, bi%d};' % (s, s, s, s, s, s))
        w('  reg signed [15:0] or%d, oi%d;' % (s, s))
        if SWAP:
            assert mode == 'V6'
            w('  always @(posedge clk) begin or%d <= ph%d ? dfr%d : dr%d; oi%d <= ph%d ? dfi%d : di%d; end' % (s, s, s, s, s, s, s, s))
        else:
            w('  always @(posedge clk) begin or%d <= ph%d ? sur%d : dr%d; oi%d <= ph%d ? sui%d : di%d; end' % (s, s, s, s, s, s, s, s))
        # delay line
        if D == 2048 and mem == 'spram':
            w('  wire [15:0] dAh%d, dAl%d, dBh%d, dBl%d;' % (s, s, s, s))
            w('  wire [13:0] aa%d = {4\'d0, n1[%d:1]}; wire [13:0] ab%d = {4\'d0, n[%d:1]};' % (s, k - 1, s, k - 1))
            for bank, addr, wren in (('A', 'aa%d' % s, '~n[0]'), ('B', 'ab%d' % s, 'n[0]')):
                for half, sl in (('h', '31:16'), ('l', '15:0')):
                    w('  SB_SPRAM256KA sp%s%s%d (.ADDRESS(%s), .DATAIN(din%d[%s]), .MASKWREN(4\'b1111), .WREN(%s), .CHIPSELECT(1\'b1), .CLOCK(clk), .STANDBY(1\'b0), .SLEEP(1\'b0), .POWEROFF(1\'b1), .DATAOUT(d%s%s%d));'
                      % (bank, half, s, addr, s, sl, wren, bank, half, s))
            w('  assign dout%d = n[0] ? {dBh%d, dBl%d} : {dAh%d, dAl%d};' % (s, s, s, s, s))
        elif D >= 8:
            w('  (* ram_style = "block" *) reg [31:0] mem%d [0:%d];' % (s, D - 1))
            w('  reg [31:0] rd%d;' % s)
            w('  always @(posedge clk) begin mem%d[n[%d:0]] <= din%d; rd%d <= mem%d[n1[%d:0]]; end' % (s, k - 1, s, s, s, k - 1))
            w('  assign dout%d = rd%d;' % (s, s))
        else:
            w('  reg [31:0] sr%d [0:%d]; integer j%d;' % (s, D - 1, s))
            w('  always @(posedge clk) begin sr%d[0] <= din%d; for (j%d = 1; j%d < %d; j%d = j%d + 1) sr%d[j%d] <= sr%d[j%d-1]; end' % (s, s, s, s, D, s, s, s, s, s, s))
            w('  assign dout%d = sr%d[%d];' % (s, s, D - 1))
        off += D + 1
        cur = ('or%d' % s, 'oi%d' % s)
        # rotator after stage s
        if s in TRIV:
            p = 12 - s
            w('  // trivial -j after stage %d (bits %d,%d of the slot)' % (s, p, p - 1))
            w('  reg [11:0] t%d; always @(posedge clk) t%d <= rst ? 12\'d%d : t%d + 1\'b1;' % (s, s, (-off) % N, s))
            w('  wire tj%d = %st%d[%d] & t%d[%d];' % (s, '~' if SWAP else '', s, p, s, p - 1))
            w('  reg signed [15:0] tr%d, ti%d;' % (s, s))
            w('  always @(posedge clk) begin tr%d <= tj%d ? %s : %s; ti%d <= tj%d ? -%s : %s; end' % (s, s, cur[1], cur[0], s, s, cur[0], cur[1]))
            off += 1
            cur = ('tr%d' % s, 'ti%d' % s)
        elif s in W16E:
            p = 12 - s
            w('  // W16 constant rotator after stage %d' % s)
            w('  reg [11:0] t%d; always @(posedge clk) t%d <= rst ? 12\'d%d : t%d + 1\'b1;' % (s, s, (-off) % N, s))
            w('  wire [1:0] K%d = %s{t%d[%d], t%d[%d]}; wire [1:0] M%d = {t%d[%d], t%d[%d]};' % (s, '~' if SWAP else '', s, p, s, p + 1, s, s, p - 1, s, p - 2))
            w('  wire [3:0] e%d = K%d * M%d;' % (s, s, s))
            # stage A: pre-rotate by (-j)^q, register r
            xr_, xi_ = cur
            w('  reg signed [15:0] ur%d, ui%d; reg [1:0] rr%d;' % (s, s, s))
            w('  always @(posedge clk) begin rr%d <= e%d[1:0];' % (s, s))
            w('    case (e%d[3:2]) 2\'d0: begin ur%d <= %s; ui%d <= %s; end 2\'d1: begin ur%d <= %s; ui%d <= -%s; end' % (s, s, xr_, s, xi_, s, xi_, s, xr_))
            w('      2\'d2: begin ur%d <= -%s; ui%d <= -%s; end default: begin ur%d <= -%s; ui%d <= %s; end endcase end' % (s, xr_, s, xi_, s, xi_, s, xr_))
            # stage B: CSD products (registered)
            W = 16 + G + int(os.environ.get('WX', '0'))  # product width (v1: G=4, +2 -> 22 bits; v6c: 20); final sums one bit wider
            for c in ('r', 'i'):
                x = 'u%s%d' % (c, s)
                if W - 16 - G > 0:
                    w('  wire signed [%d:0] xs%s%d = {{%d{%s[15]}}, %s, %d\'d0};' % (W - 1, c, s, W - 16 - G, x, x, G))
                else:
                    w('  wire signed [%d:0] xs%s%d = {%s, %d\'d0};' % (W - 1, c, s, x, G))
                if W16PIPE:   # balanced tree: level-1 pair sums registered, level-2 sums after the register (latency +1)
                    w('  wire signed [%d:0] nxs%s%d = ~xs%s%d;' % (W - 1, c, s, c, s))
                    w('  reg signed [%d:0] xsd%s%d; always @(posedge clk) xsd%s%d <= xs%s%d;' % (W - 1, c, s, c, s, c, s))
                    for nm, digs in CSD.items():
                        terms = [('(xs%s%d >>> %d)' if sg > 0 else '(nxs%s%d >>> %d)') % (c, s, kk) for sg, kk in digs]
                        negs = [1 if sg < 0 else 0 for sg, kk in digs]
                        parts = []; owed = 0
                        for j in range(0, len(terms), 2):
                            pn = 'L%s%s%d_%d' % (nm, c, s, j // 2)
                            if j + 1 < len(terms):
                                o = negs[j] + negs[j + 1]
                                cin = min(o, 1); owed += o - cin
                                w('  reg signed [%d:0] %s; always @(posedge clk) %s <= %s + %s%s;' % (W - 1, pn, pn, terms[j], terms[j + 1], (" + %d'sd1" % W) if cin else ''))
                            else:
                                owed += negs[j]
                                w('  reg signed [%d:0] %s; always @(posedge clk) %s <= %s;' % (W - 1, pn, pn, terms[j]))
                            parts.append(pn)
                        prev = parts[0]
                        for j, pn in enumerate(parts[1:]):
                            nmj = 'P%s%s%d_%d' % (nm, c, s, j)
                            cin = 1 if owed > 0 else 0; owed -= cin
                            w('  (* keep *) wire signed [%d:0] %s = %s + %s%s;' % (W - 1, nmj, prev, pn, (" + %d'sd1" % W) if cin else ''))
                            prev = nmj
                        assert owed == 0, (nm, owed)
                        w('  wire signed [%d:0] P%s%s%d = %s;' % (W - 1, nm, c, s, prev))
                elif W16CHAIN:   # explicit chain of 2-input adders; a subtracted term uses the shared ~x (+1 as carry-in)
                    w('  wire signed [%d:0] nxs%s%d = ~xs%s%d;' % (W - 1, c, s, c, s))
                    for nm, digs in CSD.items():
                        sg0, k0 = digs[0]
                        prev = '(xs%s%d >>> %d)' % (c, s, k0)
                        for j, (sg, kk) in enumerate(digs[1:]):
                            nmj = 'P%s%s%d_%d' % (nm, c, s, j)
                            if sg > 0:
                                w('  (* keep *) wire signed [%d:0] %s = %s + (xs%s%d >>> %d);' % (W - 1, nmj, prev, c, s, kk))
                            else:
                                w('  (* keep *) wire signed [%d:0] %s = %s + (nxs%s%d >>> %d) + %d\'sd1;' % (W - 1, nmj, prev, c, s, kk, W))
                            prev = nmj
                        w('  wire signed [%d:0] P%s%s%d = %s;' % (W - 1, nm, c, s, prev))
                else:
                  for nm, digs in CSD.items():
                    terms = []
                    for sg, kk in digs:
                        t = '(xs%s%d >>> %d)' % (c, s, kk)
                        terms.append(('+ ' if sg > 0 else '- ') + t)
                    expr = ' '.join(terms)
                    if expr.startswith('+ '): expr = expr[2:]
                    w('  wire signed [%d:0] P%s%s%d = %s;' % (W - 1, nm, c, s, expr))
            w('  reg signed [%d:0] Cr%d, Ci%d, Sr%d, Si%d;' % (W - 1, s, s, s, s))
            if W16PIPE:
                w('  reg [1:0] rp%d; always @(posedge clk) rp%d <= rr%d;' % (s, s, s))
                w('  always @(posedge clk) case (rp%d)' % s)
                w('    2\'d0: begin Cr%d <= xsdr%d; Ci%d <= xsdi%d; Sr%d <= 0; Si%d <= 0; end' % (s, s, s, s, s, s))
            else:
                w('  always @(posedge clk) case (rr%d)' % s)
                w('    2\'d0: begin Cr%d <= xsr%d; Ci%d <= xsi%d; Sr%d <= 0; Si%d <= 0; end' % (s, s, s, s, s, s))
            w('    2\'d1: begin Cr%d <= Pc1r%d; Ci%d <= Pc1i%d; Sr%d <= Ps1r%d; Si%d <= Ps1i%d; end' % (s, s, s, s, s, s, s, s))
            w('    2\'d2: begin Cr%d <= Pc2r%d; Ci%d <= Pc2i%d; Sr%d <= Pc2r%d; Si%d <= Pc2i%d; end' % (s, s, s, s, s, s, s, s))
            w('    default: begin Cr%d <= Ps1r%d; Ci%d <= Ps1i%d; Sr%d <= Pc1r%d; Si%d <= Pc1i%d; end endcase' % (s, s, s, s, s, s, s, s))
            # stage C: re = (Cr + Si + half) >> G ; im = (Ci - Sr + half) >> G
            w('  wire signed [%d:0] zr%d = Cr%d + Si%d + %d; wire signed [%d:0] zi%d = Ci%d - Sr%d + %d;' % (W, s, s, s, 1 << (G - 1), W, s, s, s, 1 << (G - 1)))
            w('  reg signed [15:0] wr%d, wi%d; always @(posedge clk) begin wr%d <= zr%d[%d:%d]; wi%d <= zi%d[%d:%d]; end' % (s, s, s, s, G + 15, G, s, s, G + 15, G))
            off += 4 if W16PIPE else 3
            cur = ('wr%d' % s, 'wi%d' % s)
        elif s in (4, 8):
            w('  // general rotator after stage %d' % s)
            # exponent pipeline: counter leads the data by 2 cycles
            w('  reg [11:0] t%d; always @(posedge clk) t%d <= rst ? 12\'d%d : t%d + 1\'b1;' % (s, s, (2 - off) % N, s))
            if s == 4:
                w('  wire [3:0] K%d = %s{t%d[8], t%d[9], t%d[10], t%d[11]};' % (s, '~' if SWAP else '', s, s, s, s))
                w('  wire [11:0] e%d = K%d * t%d[7:0];' % (s, s, s))
            else:
                w('  wire [3:0] K%d = %s{t%d[4], t%d[5], t%d[6], t%d[7]};' % (s, '~' if SWAP else '', s, s, s, s))
                w('  wire [7:0] f%d = K%d * t%d[3:0];' % (s, s, s))
                w('  wire [11:0] e%d = {f%d, 4\'d0};' % (s, s))
            w('  reg [11:0] er%d; always @(posedge clk) er%d <= e%d;' % (s, s, s))
            # cycle t-1: er holds e(slot of data at t). ROM address, flags
            w('  wire [9:0] rr%d = er%d[9:0];' % (s, s))
            w('  wire sw%d = rr%d > 10\'d512;' % (s, s))
            w('  wire [9:0] ix%d = sw%d ? (10\'d0 - rr%d) : rr%d;' % (s, s, s, s))  # 1024 - r (mod 1024)
            w('  wire mid%d = (ix%d == 10\'d512);' % (s, s))
            w('  reg [1:0] q%d; reg sw%dr, mid%dr; always @(posedge clk) begin q%d <= er%d[11:10]; sw%dr <= sw%d; mid%dr <= mid%d; end' % (s, s, s, s, s, s, s, s, s))
            # ROM: 512 entries (index 0..511), registered read
            if s == 4:
                w('  (* ram_style = "block" *) reg [15:0] romc%d [0:511]; (* ram_style = "block" *) reg [15:0] roms%d [0:511];' % (s, s))
                w('  initial begin')
                for i in range(512):
                    c_, s_ = TROM[i]
                    if brk == 'twiddle' and i == 300:
                        c_ ^= 1
                    w('    romc%d[%d] = 16\'h%04x; roms%d[%d] = 16\'h%04x;' % (s, i, c_ & 0xFFFF, s, i, s_ & 0xFFFF))
                w('  end')
                w('  reg [15:0] tc%d, ts%d; always @(posedge clk) begin tc%d <= romc%d[ix%d[8:0]]; ts%d <= roms%d[ix%d[8:0]]; end' % (s, s, s, s, s, s, s, s))
            else:  # W256 twiddles: 32 entries (multiples of 16), in logic
                w('  reg [15:0] tc%d, ts%d;' % (s, s))
                w('  always @(posedge clk) case (ix%d[8:4])' % s)
                for i in range(32):
                    c_, s_ = TROM[16 * i]
                    w("    5'd%d: begin tc%d <= 16'h%04x; ts%d <= 16'h%04x; end" % (i, s, c_ & 0xFFFF, s, s_ & 0xFFFF))
                w('  endcase')
            # cycle t: ROM out; swap / mid; register -> cycle t+1 at the MAC inputs
            w('  reg [15:0] Cw%d, Sw%d;' % (s, s))
            w('  always @(posedge clk) begin Cw%d <= mid%dr ? 16\'d23170 : (sw%dr ? ts%d : tc%d); Sw%d <= mid%dr ? 16\'d23170 : (sw%dr ? tc%d : ts%d); end' % (s, s, s, s, s, s, s, s, s, s))
            # data: pre-rotate by (-j)^q at cycle t -> registered at t+1
            xr_, xi_ = cur
            w('  reg signed [15:0] ur%d, ui%d;' % (s, s))
            w('  always @(posedge clk) case (q%d) 2\'d0: begin ur%d <= %s; ui%d <= %s; end 2\'d1: begin ur%d <= %s; ui%d <= -%s; end' % (s, s, xr_, s, xi_, s, xi_, s, xr_))
            w('    2\'d2: begin ur%d <= -%s; ui%d <= -%s; end default: begin ur%d <= -%s; ui%d <= %s; end endcase' % (s, xr_, s, xi_, s, xi_, s, xr_))
            if rot == 'mac':
                def mac(nm, A, B, Cin, Din, sub, pipe):
                    w("  wire [31:0] %s;" % nm)
                    w("  SB_MAC16 #(.C_REG(1'b0), .A_REG(1'b1), .B_REG(1'b1), .D_REG(1'b0), .TOP_8x8_MULT_REG(1'b0), .BOT_8x8_MULT_REG(1'b0),"
                      " .PIPELINE_16x16_MULT_REG1(1'b0), .PIPELINE_16x16_MULT_REG2(1'b%d), .TOPOUTPUT_SELECT(2'd1), .TOPADDSUB_LOWERINPUT(2'd2),"
                      " .TOPADDSUB_UPPERINPUT(1'b1), .TOPADDSUB_CARRYSELECT(2'd2), .BOTOUTPUT_SELECT(2'd1), .BOTADDSUB_LOWERINPUT(2'd2),"
                      " .BOTADDSUB_UPPERINPUT(1'b1), .BOTADDSUB_CARRYSELECT(2'd0), .MODE_8x8(1'b0), .A_SIGNED(1'b1), .B_SIGNED(1'b1))" % pipe)
                    w("    u_%s (.CLK(clk), .CE(1'b1), .C(%s), .A(%s), .B(%s), .D(%s), .AHOLD(1'b0), .BHOLD(1'b0), .CHOLD(1'b0), .DHOLD(1'b0),"
                      " .IRSTTOP(1'b0), .IRSTBOT(1'b0), .ORSTTOP(1'b0), .ORSTBOT(1'b0), .OLOADTOP(1'b0), .OLOADBOT(1'b0),"
                      " .ADDSUBTOP(1'b%d), .ADDSUBBOT(1'b%d), .OHOLDTOP(1'b0), .OHOLDBOT(1'b0), .CI(1'b0), .ACCUMCI(1'b0), .SIGNEXTIN(1'b0), .O(%s),"
                      " .CO(), .ACCUMCO(), .SIGNEXTOUT());" % (nm, Cin, A, B, Din, sub, sub, nm))
                mac('ma%d' % s, 'ur%d' % s, 'Cw%d' % s, "16'd0", "16'd16384", 0, 0)
                mac('mb%d' % s, 'ui%d' % s, 'Sw%d' % s, 'ma%d[31:16]' % s, 'ma%d[15:0]' % s, 0, 1)
                mac('mc%d' % s, 'ui%d' % s, 'Cw%d' % s, "16'd0", "16'd16384", 0, 0)
                mac('md%d' % s, 'ur%d' % s, 'Sw%d' % s, 'mc%d[31:16]' % s, 'mc%d[15:0]' % s, 1, 1)
                w('  wire signed [15:0] gr%d = mb%d[30:15], gi%d = md%d[30:15];' % (s, s, s, s))
            else:  # logic-cell multipliers, same arithmetic and latency (3 cycles)
                w('  reg signed [15:0] ur%dd, ui%dd, Cw%dd, Sw%dd; always @(posedge clk) begin ur%dd <= ur%d; ui%dd <= ui%d; Cw%dd <= Cw%d; Sw%dd <= Sw%d; end' % ((s,) * 4 + (s, s) * 4))
                w('  reg signed [31:0] pa%d, pc%d, pb%d, pd%d; always @(posedge clk) begin pa%d <= ur%dd * $signed(Cw%dd) + 32\'sd16384; pc%d <= ui%dd * $signed(Cw%dd) + 32\'sd16384; pb%d <= ui%dd * $signed(Sw%dd); pd%d <= ur%dd * $signed(Sw%dd); end' % ((s,) * 4 + (s,) * 12))
                w('  reg signed [31:0] qa%d, qc%d; always @(posedge clk) begin qa%d <= pa%d + pb%d; qc%d <= pc%d - pd%d; end' % ((s,) * 8))
                w('  wire signed [15:0] gr%d = qa%d[30:15], gi%d = qc%d[30:15];' % (s, s, s, s))
            off += 4
            cur = ('gr%d' % s, 'gi%d' % s)
    w('  assign yr = %s; assign yi = %s;' % cur)
    w("  assign y_first = (n == 12'd%d);" % (off % N))
    w('endmodule')
    return '\n'.join(L) + '\n', off

# ---------- simulation ----------
CELLS = '$OSS_CAD_SUITE/share/yosys/ice40/cells_sim.v'
BIN = '$OSS_CAD_SUITE/bin'

def frames(total, seed):
    rng = np.random.default_rng(seed)
    n = np.arange(N)
    F = []
    F.append((np.full(N, 32767), np.zeros(N)))                       # full-scale DC (real)
    F.append((np.full(N, -32768), np.full(N, -32768)))               # most negative, both components
    F.append((np.where(n % 2 == 0, 32767, -32768), np.zeros(N)))     # alternating (Nyquist)
    F.append((np.round(32767 * np.cos(2 * np.pi * 17 * n / N)), np.round(32767 * np.sin(2 * np.pi * 17 * n / N))))  # bin-centre tone
    F.append((np.round(32000 * np.cos(2 * np.pi * 2047.5 * n / N)), np.zeros(N)))  # tone at a bin edge near Nyquist
    e = np.zeros(N); e[0] = 32767; F.append((e, np.zeros(N)))           # impulse
    F.append((rng.choice([-32768, 32767, -1, 0, 1], N), rng.choice([-32768, 32767, -1, 0, 1], N)))
    while len(F) * N < total:
        if len(F) % 2:
            F.append((rng.integers(-32768, 32768, N), rng.integers(-32768, 32768, N)))
        else:
            F.append((np.round(rng.normal(0, 8192, N)), np.round(rng.normal(0, 8192, N))))
    return [(np.clip(a, -32768, 32767).astype(np.int64), np.clip(b, -32768, 32767).astype(np.int64)) for a, b in F]

def sim(design, total, work, mode, brk_model=None, netlist=False):
    os.makedirs(work, exist_ok=True)
    F = frames(total, 1)
    xr = np.concatenate([f[0] for f in F]); xi = np.concatenate([f[1] for f in F])
    lat = 13000
    with open(os.path.join(work, 'in.hex'), 'w') as fh:
        for a, b in zip(xr, xi):
            fh.write('%04x%04x\n' % (a & 0xFFFF, b & 0xFFFF))
        for _ in range(lat):
            fh.write('00000000\n')
    tot = len(xr) + lat
    tb = """`timescale 1ns/1ps
module tb;
  reg clk = 0, rst = 1; reg [31:0] mem [0:%d]; reg [15:0] xr = 0, xi = 0; wire [15:0] yr, yi; wire y_first;
  integer t, fo;
  fft4k dut (.clk(clk), .rst(rst), .xr(xr), .xi(xi), .yr(yr), .yi(yi), .y_first(y_first));
  initial begin
    $readmemh("in.hex", mem); fo = $fopen("out.txt", "w");
    rst = 1; #5 clk = 1; #5 clk = 0; #5 clk = 1; #5 clk = 0; rst = 0;
    for (t = 0; t < %d; t = t + 1) begin
      xr = mem[t][31:16]; xi = mem[t][15:0];
      #1 $fwrite(fo, "%%h %%h %%b\\n", yr, yi, y_first);
      #4 clk = 1; #5 clk = 0;
    end
    $fclose(fo); $finish;
  end
endmodule
""" % (tot - 1, tot)
    open(os.path.join(work, 'tb.v'), 'w').write(tb)
    exe = os.path.join(work, 'sim.vvp')
    srcs = [os.path.join(work, 'tb.v'), os.path.abspath(design), CELLS]
    r = subprocess.run(['/usr/sbin/taskpolicy', '-b', BIN + '/iverilog', '-g2012', '-o', exe, '-DNO_ICE40_DEFAULT_ASSIGNMENTS'] + srcs, capture_output=True, text=True, cwd=work)
    if r.returncode:
        print(r.stdout[-3000:], r.stderr[-3000:]); sys.exit(2)
    r = subprocess.run(['/usr/sbin/taskpolicy', '-b', BIN + '/vvp', '-n', exe], capture_output=True, text=True, cwd=work)
    lines = [ln.split() for ln in open(os.path.join(work, 'out.txt')).read().split('\n') if ln.strip()]
    def sv(h):
        if 'x' in h or 'z' in h: return None
        v = int(h, 16); return v - 65536 if v >= 32768 else v
    ys = [(sv(a), sv(b)) for a, b, c in lines]; fs = [c == '1' for a, b, c in lines]
    marks = [i for i, f in enumerate(fs) if f]
    FR = np.stack([f[0] for f in F]); FI = np.stack([f[1] for f in F])
    er, ei = model(FR, FI, mode, brk_model)
    base = None
    for m in marks[:4]:
        seg = ys[m:m + N]
        if len(seg) == N and all(v[0] is not None for v in seg):
            if np.array_equal(np.array([v[0] for v in seg]), er[0]) and np.array_equal(np.array([v[1] for v in seg]), ei[0]):
                base = m; break
    if base is None:
        print('RESULT FAIL: no marker aligns frame 0; markers', marks[:4])
        if marks:
            m = marks[1] if len(marks) > 1 else marks[0]
            seg = ys[m:m + 8]
            print('got', seg, 'expected', list(zip(er[0][:8].tolist(), ei[0][:8].tolist())))
            for mm in marks[:3]:
                seg = ys[mm:mm + N]
                got = np.array([v[0] if v[0] is not None else -99999 for v in seg])
                if len(got) == N:
                    print('marker', mm, 'mismatches re', int((got != er[0]).sum()))
        sys.exit(1)
    ok = bad = 0; first_bad = None
    for f in range(len(F)):
        seg = ys[base + f * N: base + (f + 1) * N]
        if len(seg) < N: break
        gr = np.array([-99999 if v[0] is None else v[0] for v in seg]); gi = np.array([-99999 if v[1] is None else v[1] for v in seg])
        if np.array_equal(gr, er[f]) and np.array_equal(gi, ei[f]) and fs[base + f * N]:
            ok += 1
        else:
            bad += 1
            if first_bad is None:
                j = int(np.nonzero((gr != er[f]) | (gi != ei[f]))[0][0]) if not (np.array_equal(gr, er[f]) and np.array_equal(gi, ei[f])) else -1
                first_bad = (f, j, (int(gr[j]), int(gi[j])) if j >= 0 else None, (int(er[f][j]), int(ei[f][j])) if j >= 0 else None)
    print('frames %d, input samples %d, ok %d, bad %d, first_bad %s, base line %d' % (len(F), len(xr), ok, bad, first_bad, base))
    print('RESULT', 'PASS' if bad == 0 and ok == len(F) else 'FAIL')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd')
    ap.add_argument('design', nargs='?')
    ap.add_argument('-o')
    ap.add_argument('--round', default='V6')
    ap.add_argument('--mem', default='spram')
    ap.add_argument('--rot', default='mac')
    ap.add_argument('--break', dest='brk', default=None)
    ap.add_argument('--frames', type=int, default=16)
    ap.add_argument('--samples', type=int, default=110000)
    ap.add_argument('--work', default='simwork')
    ap.add_argument('--model-break', default=None)
    a = ap.parse_args()
    if a.cmd == 'model':
        precision(a.frames)
    elif a.cmd == 'gen':
        v, off = gen(a.round, a.mem, a.rot, a.brk)
        open(a.o, 'w').write(v)
        print('wrote', a.o, 'latency', off, 'lines', v.count('\n'))
    elif a.cmd == 'sim':
        sim(a.design, a.samples, a.work, a.round, a.model_break)

if __name__ == '__main__':
    main()
