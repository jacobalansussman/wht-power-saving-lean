#!/usr/bin/env python3
"""up5k_sta.py -- static timing for iCE40-UP5K designs that use SB_MAC16 (DSP) and SB_SPRAM256KA blocks.

WHY: nextpnr-ice40 (0.11.1-59, ice40/arch.cc) treats EVERY non-clock pin of ICESTORM_DSP and ICESTORM_SPRAM
as a register pin: inputs get a 0.1 ns setup, DSP outputs a 0.1 ns clock-to-out (the chipdb has no MAC16
timing), and no combinational arc crosses the block, whatever A_REG/B_REG/.../TOPOUTPUT_SELECT say. A DSP used
as a combinational adder or multiplier is therefore cut into two short half-paths, and a registered DSP output
starts 1.2 to 1.9 ns too early. SPRAM is better: its CLOCK->DATAOUT (1.82 ns) IS in the chipdb; only its setups
(0.14-0.40 ns real, 0.1 used) are placeholders.

WHAT: re-time the routed design. nextpnr writes an SDF file (--sdf) holding every routed net delay
(INTERCONNECT, one per driver->sink pin) and every cell arc it uses; this script keeps all of that (so a design
without DSP/SPRAM gives nextpnr's own figure), drops the placeholder arcs of DSP and SPRAM cells, and inserts
arcs built per cell from its parameters (packed JSON, --write) and from Lattice's characterised delays in
icestorm's timings_up5k.txt (iCEcube2-derived; worst corner, max of rise/fall, exactly as nextpnr's chipdb.py
reads that file for logic cells). Then it runs an ordinary longest-path analysis (one clock, rising edge).

USAGE
  nextpnr-ice40 --up5k ... --sdf D.sdf --write D.packed.json        # (common/fpga_measure.py: --sdf flag, or run yourself)
  python3 -I up5k_sta.py D.sdf D.packed.json [--paths 3] [--json out.json]
Prints: nextpnr-model Fmax (recomputed from the SDF as nextpnr has it; must equal nextpnr's log),
corrected Fmax, the corrected critical paths with every arc, and the internal DSP register-to-register
estimates. Delays in ns.

LIMITS (see the README of this folder): one clock domain, rising edge; clock-network
skew ignored (nextpnr ignores it too here); the DSP model is exact (per-bit characterised arcs) for the
characterised configurations (adder bypass, accumulator, multiplier bypass/pipelined) and DERIVED for the
mixed ones (marked 'derived' in the path listing): an unregistered C/D pin into a registered adder output,
a product feeding the adder combinationally, internal registers driving a combinational output.
"""
import os, sys, re, json, argparse, collections

TDB_PATH = os.path.join(os.environ.get('OSS_CAD_SUITE', '/usr/local'), 'share', 'icebox', 'timings_up5k.txt')   # icestorm's database
T_INT_CQ = 0.60   # [guess] clock-to-q of a DSP-internal register before its output mux (only for derived arcs)

# ---------------------------------------------------------------- Lattice/icestorm timing database
def _bitname(s):
    return re.sub(r'\[(\d+)\]', r'_\1', s.split(':')[-1])

def load_tdb(path=TDB_PATH):
    db = collections.defaultdict(lambda: {'IOPATH': {}, 'SETUP': {}, 'RECOVERY': {}})
    cur = None
    for ln in open(path):
        t = ln.split()
        if not t:
            continue
        if t[0] == 'CELL':
            cur = t[1]; continue
        if cur is None or not (cur.startswith('SB_MAC16') or cur in ('SB_SPRAM256KA', 'SB_RAM40_4K')):
            continue
        if t[0] in ('IOPATH', 'SETUP', 'RECOVERY'):
            vals = []
            for f in t[3:5]:
                p = f.split(':')
                if len(p) == 3 and p[2] not in ('*', ''):
                    vals.append(float(p[2]) / 1000.0)
            if not vals:
                continue
            v = max(vals)
            a, b = _bitname(t[1]), _bitname(t[2])
            if t[0] == 'IOPATH':
                k = (a, b); d = db[cur]['IOPATH']
            else:
                k = a; d = db[cur][t[0]]
            d[k] = max(d.get(k, -1e9), v)
    return db

# ---------------------------------------------------------------- SDF (nextpnr's own output)
_tok = re.compile(r'\(|\)|(?:\\.|[^\s()\\])+')

def _parse_sexp(text):
    stack = [[]]
    for m in _tok.finditer(text):
        s = m.group(0)
        if s == '(':
            stack.append([])
        elif s == ')':
            x = stack.pop(); stack[-1].append(x)
        else:
            stack[-1].append(s)
    return stack[0]

def _unesc(s):
    return re.sub(r'\\(.)', r'\1', s)

def _split_pin(s):
    # "cell/port" with escaped characters; split at the last unescaped '/'
    i = len(s) - 1
    while i >= 0:
        if s[i] == '/' and (i == 0 or s[i - 1] != '\\'):
            return _unesc(s[:i]), _unesc(s[i + 1:])
        i -= 1
    raise ValueError(s)

def _dly(lst):
    # [(a:b:c)] lists -> max of third fields, ns
    v = 0.0
    for x in lst:
        if isinstance(x, list) and x:
            p = x[0].split(':')
            v = max(v, float(p[-1]) / 1000.0)
    return v

def _port(x):
    # 'PORT' or ['posedge','PORT']
    if isinstance(x, list):
        return x[1], x[0]
    return x, None

# nextpnr writes a register's clock-to-out as a plain IOPATH from its clock port
CLOCK_PORTS = {'CLK', 'RCLK', 'WCLK', 'CLOCK', 'INPUT_CLK', 'OUTPUT_CLK', 'SBCLKI'}

def load_sdf(path):
    tree = _parse_sexp(open(path).read())
    df = tree[0]
    cells = []
    for e in df:
        if isinstance(e, list) and e and e[0] == 'CELL':
            ctype = inst = None; arcs = []; setups = []
            for f in e[1:]:
                if f[0] == 'CELLTYPE':
                    ctype = f[1].strip('"')
                elif f[0] == 'INSTANCE':
                    inst = _unesc(f[1]) if len(f) > 1 else ''
                elif f[0] == 'DELAY':
                    for g in f[1:]:
                        for h in g[1:]:
                            if h[0] == 'INTERCONNECT':
                                arcs.append(('NET', _split_pin(h[1]), _split_pin(h[2]), _dly(h[3:])))
                            elif h[0] == 'IOPATH':
                                fp, edge = _port(h[1])
                                arcs.append(('CQ' if (edge or fp in CLOCK_PORTS) else 'COMB', fp, h[2], _dly(h[3:])))
                elif f[0] == 'TIMINGCHECK':
                    for g in f[1:]:
                        if g[0] in ('SETUPHOLD', 'SETUP'):
                            p, _ = _port(g[1])
                            ck, _ = _port(g[2])
                            setups.append((p, _dly(g[3:4]), ck))
            cells.append((ctype, inst, arcs, setups))
    return cells

# ---------------------------------------------------------------- the SB_MAC16 model
def _pint(v):
    if isinstance(v, int):
        return v
    s = str(v).strip()
    if s and set(s) <= set('01xz'):
        return int(s.replace('x', '0').replace('z', '0'), 2)
    try:
        return int(s)
    except ValueError:
        return 0

class DSPModel:
    def __init__(self, tdb):
        self.t = tdb
        g = lambda m: tdb['SB_MAC16_' + m]
        self.ADS16, self.ADS32 = g('ADS_U_16P16_BYPASS'), g('ADS_U_32P32_BYPASS')
        self.ACC16, self.ACC32 = g('ACC_U_16P16_BYPASS'), g('ACC_U_32P32_BYPASS')
        self.PIPE = g('ADS_U_32P32_ALL_PIPELINE')          # registered pins, registered outputs
        self.MAC16B, self.MAC8B = g('MAC_U_16X16_BYPASS'), g('MAC_U_8X8_BYPASS')
        # Calibration [derived]: how much shorter "pin -> adder -> accumulator register (setup)" is than
        # "pin -> adder -> O pin (combinational)". Smallest difference over the lower-input bits (conservative).
        ds = []
        for half, pin, osel in (('top', 'A', range(16, 32)), ('bot', 'B', range(0, 16))):
            for i in range(16):
                c = max(self.ADS16['IOPATH'].get((f'{pin}_{i}', f'O_{j}'), 0) for j in osel)
                s = self.ACC16['SETUP'].get(f'{pin}_{i}')
                if c and s is not None:
                    ds.append(c - s)
        self.DELTA_ADDER = min(ds)
        # Output-mux part of a combinational multiplier arc [derived]: MUL bypass A->O minus MAC bypass A->CLK
        # would mix two paths; use the adder figure as the mux estimate (same output mux).
        self.DELTA_OMUX = self.DELTA_ADDER

    def mul_mode(self, p, eight):
        s = 'S' if (p['A_SIGNED'] or p['B_SIGNED']) else 'U'
        return self.t['SB_MAC16_MUL_%s_%s_BYPASS' % (s, '8X8' if eight else '16X16')], \
               self.t['SB_MAC16_MUL_%s_%s_ALL_PIPELINE' % (s, '8X8' if eight else '16X16')]

    def arcs(self, params):
        """Return comb arcs [(from,to,d,tag)], starts {out: (cq,tag)}, ends {in: (setup,tag)}, internal [(desc,d)]."""
        P = {k: _pint(params.get(k, 0)) for k in (
            'A_REG', 'B_REG', 'C_REG', 'D_REG', 'TOP_8x8_MULT_REG', 'BOT_8x8_MULT_REG', 'PIPELINE_16x16_MULT_REG1',
            'PIPELINE_16x16_MULT_REG2', 'TOPOUTPUT_SELECT', 'TOPADDSUB_LOWERINPUT', 'TOPADDSUB_UPPERINPUT',
            'TOPADDSUB_CARRYSELECT', 'BOTOUTPUT_SELECT', 'BOTADDSUB_LOWERINPUT', 'BOTADDSUB_UPPERINPUT',
            'BOTADDSUB_CARRYSELECT', 'MODE_8x8', 'A_SIGNED', 'B_SIGNED', 'NEG_TRIGGER')}
        comb, starts, ends, internal = [], {}, {}, []
        cascade = P['TOPADDSUB_CARRYSELECT'] in (2, 3)
        ADS = self.ADS32 if cascade else self.ADS16
        ACC = self.ACC32 if cascade else self.ACC16
        reg = {'A': P['A_REG'], 'B': P['B_REG'], 'C': P['C_REG'], 'D': P['D_REG']}
        PIPE = self.PIPE

        def setup_pin(pin, i, val, tag):
            k = f'{pin}_{i}' if i is not None else pin
            if k not in ends or ends[k][0] < val:
                ends[k] = (val, tag)

        # input-register setups and control pins (characterised) [read]
        for pin in 'ABCD':
            if reg[pin]:
                for i in range(16):
                    setup_pin(pin, i, max(0.0, PIPE['SETUP'].get(f'{pin}_{i}', 0.01)), 'read:input-reg setup')
        for c in ('AHOLD', 'BHOLD', 'CHOLD', 'DHOLD', 'OHOLDTOP', 'OHOLDBOT'):
            setup_pin(c, None, PIPE['SETUP'].get(c, 0.45), 'read')
        for c in ('IRSTTOP', 'IRSTBOT', 'ORSTTOP', 'ORSTBOT'):
            setup_pin(c, None, PIPE['RECOVERY'].get(c, 1.2), 'read:recovery')
        setup_pin('CE', None, 0.46, 'guess:CE setup (not in database)')

        for half in ('top', 'bot'):
            T = half == 'top'
            sel = P['TOPOUTPUT_SELECT' if T else 'BOTOUTPUT_SELECT']
            lower = P['TOPADDSUB_LOWERINPUT' if T else 'BOTADDSUB_LOWERINPUT']
            upper_is_pin = P['TOPADDSUB_UPPERINPUT' if T else 'BOTADDSUB_UPPERINPUT']
            upin, lpin = ('C', 'A') if T else ('D', 'B')
            outs = [f'O_{j}' for j in (range(16, 32) if T else range(0, 16))]
            oload, addsub = ('OLOADTOP', 'ADDSUBTOP') if T else ('OLOADBOT', 'ADDSUBBOT')
            # which pins feed this half's adder, combinationally, with the database arcs that describe them
            add_src = []                                  # (pin, bits, regd)
            if upper_is_pin:
                add_src.append((upin, 16, reg[upin]))
            if lower == 0:
                add_src.append((lpin, 16, reg[lpin]))
            mult_into_adder = lower in (1, 2)
            if T and cascade:                             # bottom adder's inputs reach the top half via the carry
                if P['BOTADDSUB_UPPERINPUT']:
                    add_src.append(('D', 16, reg['D']))
                if P['BOTADDSUB_LOWERINPUT'] == 0:
                    add_src.append(('B', 16, reg['B']))
            eight = bool(P['MODE_8x8']) or lower == 1
            mregs = (P['TOP_8x8_MULT_REG'] if T else P['BOT_8x8_MULT_REG']) if lower == 1 else \
                    (P['PIPELINE_16x16_MULT_REG1'] or P['PIPELINE_16x16_MULT_REG2'])
            mulB, mulP = self.mul_mode(P, lower == 1 or (sel == 2))

            if sel == 0:                                  # combinational adder output
                for pin, n, r in add_src:
                    for i in range(n):
                        for o in outs:
                            d = ADS['IOPATH'].get((f'{pin}_{i}', o))
                            if d is None:
                                continue
                            if r:
                                starts[o] = max(starts.get(o, (0, '')), (T_INT_CQ + d, 'derived:input reg -> adder -> O'))
                            else:
                                comb.append((f'{pin}_{i}', o, d, 'read:ADS bypass'))
                for c in (oload, addsub) + (('ADDSUBBOT',) if T and cascade else ()):
                    for o in outs:
                        d = ADS['IOPATH'].get((c, o))
                        if d is not None:
                            comb.append((c, o, d, 'read:ADS bypass'))
                if P['BOTADDSUB_CARRYSELECT'] in (2, 3) and (not T or cascade):
                    cpin = 'ACCUMCI' if P['BOTADDSUB_CARRYSELECT'] == 2 else 'CI'
                    for o in outs:
                        d = self.ADS32['IOPATH'].get(('CI', o))
                        if d is not None:
                            comb.append((cpin, o, d, 'read:ADS32 CI'))
                if not upper_is_pin:                      # accumulator feedback rQ -> adder -> O
                    for o in outs:
                        starts[o] = max(starts.get(o, (0, '')), (T_INT_CQ + max(ADS['IOPATH'].get((f'{upin}_{i}', o), 0) for i in range(16)), 'derived:acc reg -> adder -> O'))
                if mult_into_adder:
                    mmax = {pin: max(v for (a, b), v in mulB['IOPATH'].items() if a.startswith(pin + '_')) for pin in 'AB'}
                    lmax = {o: max(ADS['IOPATH'].get((f'{lpin}_{i}', o), 0) for i in range(16)) for o in outs}
                    for o in outs:
                        if reg['A'] or reg['B'] or mregs:
                            starts[o] = max(starts.get(o, (0, '')), (T_INT_CQ + max(mmax.values()) - self.DELTA_OMUX + lmax[o], 'derived:reg -> mult -> adder -> O'))
                        else:
                            for pin in 'AB':
                                for i in range(16):
                                    comb.append((f'{pin}_{i}', o, mmax[pin] - self.DELTA_OMUX + lmax[o], 'derived:mult -> adder -> O'))
            elif sel == 1:                                # registered accumulator output
                for o in outs:
                    starts[o] = (PIPE['IOPATH'].get(('CLK', o), 1.98), 'read:acc reg CLK->O')
                for pin, n, r in add_src:
                    if r:
                        internal.append((f'{half}: {pin} input reg -> adder -> acc reg', T_INT_CQ + max(ACC['SETUP'].get(f'A_{i}', 0) for i in range(16))))
                        continue
                    for i in range(n):
                        if pin == lpin and not (T and cascade and pin in 'BD'):
                            s = ACC['SETUP'].get(f'{pin}_{i}')
                            tag = 'read:ACC bypass setup'
                        else:
                            s = None
                        if s is None:
                            s = max(ADS['IOPATH'].get((f'{pin}_{i}', o), 0) for o in outs) - self.DELTA_ADDER
                            tag = 'derived:ADS comb - %.2f' % self.DELTA_ADDER
                        setup_pin(pin, i, s, tag)
                setup_pin(oload, None, PIPE['SETUP'].get(oload, 0.46), 'read')
                setup_pin(addsub, None, PIPE['SETUP'].get(addsub, 2.65), 'read')
                if T and cascade:
                    setup_pin('ADDSUBBOT', None, PIPE['SETUP'].get('ADDSUBBOT', 2.65), 'read')
                if P['BOTADDSUB_CARRYSELECT'] in (2, 3) and (not T or cascade):
                    setup_pin('ACCUMCI' if P['BOTADDSUB_CARRYSELECT'] == 2 else 'CI', None, self.ACC32['SETUP'].get('CI', 1.69), 'read')
                if not upper_is_pin:
                    internal.append((f'{half}: accumulator feedback', T_INT_CQ + max(ACC['SETUP'].get(f'A_{i}', 0) for i in range(16))))
                if mult_into_adder:
                    big = self.MAC8B if eight else self.MAC16B
                    if reg['A'] or reg['B'] or mregs:
                        for pin in 'AB':
                            if not reg[pin] and not mregs:
                                for i in range(16):
                                    setup_pin(pin, i, big['SETUP'].get(f'{pin}_{i}', 6.5), 'read:MAC bypass setup')
                        if mregs:
                            pre = max(v for (a, b), v in mulB['IOPATH'].items()) - self.DELTA_OMUX
                            internal.append((f'{half}: A/B -> multiplier -> pipeline reg', T_INT_CQ + pre))
                            internal.append((f'{half}: pipeline reg -> adder -> acc reg', T_INT_CQ + max(ACC['SETUP'].get(f'A_{i}', 0) for i in range(16))))
                            for pin in 'AB':
                                if not reg[pin]:
                                    for i in range(16):
                                        setup_pin(pin, i, pre, 'derived:MUL comb - %.2f' % self.DELTA_OMUX)
                        else:
                            internal.append((f'{half}: A/B input reg -> multiplier -> adder -> acc reg', T_INT_CQ + max(big['SETUP'].get(f'A_{i}', 0) for i in range(16))))
                    else:
                        for pin in 'AB':
                            for i in range(16):
                                setup_pin(pin, i, big['SETUP'].get(f'{pin}_{i}', 6.5), 'read:MAC bypass setup')
            else:                                         # sel 2 / 3: product straight out
                registered = (P['TOP_8x8_MULT_REG'] if T else P['BOT_8x8_MULT_REG']) if sel == 2 else P['PIPELINE_16x16_MULT_REG2']
                if registered:
                    for o in outs:
                        starts[o] = (mulP['IOPATH'].get(('CLK', o), 2.12), 'read:product reg CLK->O')
                    pre = max(v for (a, b), v in mulB['IOPATH'].items()) - self.DELTA_OMUX
                    for pin in 'AB':
                        if reg[pin]:
                            internal.append((f'{half}: {pin} input reg -> multiplier -> product reg', T_INT_CQ + pre))
                        else:
                            for i in range(16):
                                setup_pin(pin, i, pre, 'derived:MUL comb - %.2f' % self.DELTA_OMUX)
                else:
                    for pin in 'AB':
                        for i in range(16):
                            for o in outs:
                                d = mulB['IOPATH'].get((f'{pin}_{i}', o))
                                if d is None:
                                    continue
                                if reg[pin] or (sel == 3 and P['PIPELINE_16x16_MULT_REG1']):
                                    starts[o] = max(starts.get(o, (0, '')), (T_INT_CQ + d, 'derived:reg -> multiplier -> O'))
                                else:
                                    comb.append((f'{pin}_{i}', o, d, 'read:MUL bypass'))
        return comb, starts, ends, internal

SPRAM_SETUP = {'ADDRESS': 'ADDRESS', 'DATAIN': 'DATAIN', 'MASKWREN': 'MASKWREN', 'WREN': 'WREN',
               'CHIPSELECT': 'CHIPSELECT', 'STANDBY': 'STANDBY', 'SLEEP': 'SLEEP', 'POWEROFF': 'POWEROFF'}

def spram_model(tdb):
    sp = tdb['SB_SPRAM256KA']
    starts = {f'DATAOUT_{i}': (sp['IOPATH'].get(('CLOCK', f'DATAOUT_{i}'), 1.821), 'read:SPRAM CLOCK->DATAOUT') for i in range(16)}
    ends = {}
    for k, v in sp['SETUP'].items():
        ends[k] = (v, 'read:SPRAM setup')
    ends.setdefault('POWEROFF', (0.1, 'guess'))
    return starts, ends

# ---------------------------------------------------------------- the analysis
def analyse(sdf_path, packed_path, tdb_path=TDB_PATH, correct=True, npaths=3, equal_clock=False, ebr_setups=False):
    tdb = load_tdb(tdb_path)
    dsp = DSPModel(tdb)
    cells = load_sdf(sdf_path)
    params = {}
    if packed_path:
        m = json.load(open(packed_path))['modules']
        top = [v for v in m.values() if v.get('attributes', {}).get('top')] or list(m.values())
        for n, c in top[0]['cells'].items():
            if c['type'] in ('ICESTORM_DSP', 'ICESTORM_SPRAM'):
                params[n] = c.get('parameters', {})
    edges = collections.defaultdict(list)
    starts, ends = {}, {}
    internal = []
    unused = set()
    hard_insts = set()
    modelled = collections.Counter()
    # clock insertion delay at each clock pin = routed delay of the clock net into it (nextpnr's clock-skew model)
    clk_arr = {}
    for ctype, inst, arcs, setups in cells:
        for a in arcs:
            if a[0] == 'NET' and a[2][1] in CLOCK_PORTS:
                clk_arr[a[2]] = a[3]
    lc_clk = sorted(v for (c, p), v in clk_arr.items() if p == 'CLK' and v > 0) or [0.0]
    lc_clk = lc_clk[len(lc_clk) // 2]
    hard_ck = {}
    for ctype, inst, arcs, setups in cells:
        if ctype in ('ICESTORM_DSP', 'ICESTORM_SPRAM'):
            ck = 'CLOCK' if ctype == 'ICESTORM_SPRAM' else 'CLK'
            hard_ck[inst] = ck
            if equal_clock:
                clk_arr[(inst, ck)] = lc_clk
    def ca(inst, ck):
        return clk_arr.get((inst, ck), 0.0)
    for ctype, inst, arcs, setups in cells:
        special = correct and ctype in ('ICESTORM_DSP', 'ICESTORM_SPRAM')
        for a in arcs:
            if a[0] == 'NET':
                (c1, p1), (c2, p2), d = a[1], a[2], a[3]
                edges[(c1, p1)].append(((c2, p2), d, 'net'))
            elif special:
                continue
            elif a[0] == 'COMB':
                edges[(inst, a[1])].append(((inst, a[2]), a[3], ctype))
            else:
                k = (inst, a[2])
                starts[k] = max(starts.get(k, (-1e9, '')), (ca(inst, a[1]) + a[3], ctype + ' clk->out (+clock net %.3f)' % ca(inst, a[1])))
        if ctype in ('ICESTORM_DSP', 'ICESTORM_SPRAM'):
            hard_insts.add(inst)
        if not special:
            for p, s, ck in setups:
                k = (inst, p)
                tag = ctype + ' setup (-clock net %.3f)' % ca(inst, ck)
                if ebr_setups and ctype == 'ICESTORM_RAM':
                    real = tdb['SB_RAM40_4K']['SETUP'].get(p, tdb['SB_RAM40_4K']['SETUP'].get(re.sub(r'_\d+$', '', p)))
                    if real is not None:
                        s = real; tag = 'read:EBR setup (-clock net %.3f)' % ca(inst, ck)
                ends[k] = max(ends.get(k, (-1e9, '')), (s - ca(inst, ck), tag))
            continue
        hc = ca(inst, hard_ck[inst])
        modelled[ctype] += 1
        if ctype == 'ICESTORM_SPRAM':
            st, en = spram_model(tdb)
            for p, v in st.items():
                starts[(inst, p)] = (v[0] + hc, v[1])
            sp_ports = {p for p, s, ck in setups}
            for p in sp_ports:
                base = re.sub(r'_\d+$', '', p)
                v = en.get(p, en.get(base, (0.1, 'nextpnr placeholder (not in database)')))
                ends[(inst, p)] = (v[0] - hc, v[1])
        else:
            if inst not in params:
                raise SystemExit('DSP %s not found in packed JSON' % inst)
            comb, st, en, intl = dsp.arcs(params[inst])
            for f, t, d, tag in comb:
                edges[(inst, f)].append(((inst, t), d, 'DSP ' + tag))
            for p, v in st.items():
                starts[(inst, p)] = (v[0] + hc, v[1])
            comb_from = {f for f, t, d, tag in comb}
            for p, s, ck in setups:
                if p in en:
                    ends[(inst, p)] = (en[p][0] - hc, en[p][1])
                elif p not in comb_from:
                    unused.add((inst, p))       # pin feeds nothing in this configuration (or the model misses it)
            for desc, d in intl:
                internal.append((inst, desc, d))
    # longest path, Kahn order
    indeg = collections.Counter()
    nodes = set(edges) | set(starts) | set(ends)
    for u, lst in edges.items():
        for v, d, k in lst:
            indeg[v] += 1; nodes.add(v)
    arr, pred = {}, {}
    arr_h, pred_h = {}, {}     # arrival over paths that start in, or cross, a DSP/SPRAM block
    for k, (d, tag) in starts.items():
        arr[k] = d; pred[k] = (None, d, tag)
        if k[0] in hard_insts:
            arr_h[k] = d; pred_h[k] = (None, d, tag)
    q = collections.deque([n for n in nodes if indeg[n] == 0])
    seen = 0
    while q:
        u = q.popleft(); seen += 1
        for v, d, kind in edges.get(u, ()):
            if u in arr and (v not in arr or arr[u] + d > arr[v]) and not (v in starts and pred.get(v, (1,))[0] is None):
                arr[v] = arr[u] + d; pred[v] = (u, d, kind)
            if not (v in starts and pred.get(v, (1,))[0] is None):
                cand = None
                if u in arr_h:
                    cand = (arr_h[u] + d, 'h')
                if kind.startswith('DSP') and u in arr and (cand is None or arr[u] + d > cand[0]):
                    cand = (arr[u] + d, 'a')
                if cand and (v not in arr_h or cand[0] > arr_h[v]):
                    arr_h[v] = cand[0]; pred_h[v] = (u, d, kind, cand[1])
            indeg[v] -= 1
            if indeg[v] == 0:
                q.append(v)
    loops = len(nodes) - seen
    res = []
    for k, (s, tag) in ends.items():
        if k in arr:
            res.append((arr[k] + s, k, s, tag))
    for inst, desc, d in internal:
        res.append((d, (inst, 'internal: ' + desc), 0.0, 'derived:internal DSP path'))
    res.sort(key=lambda x: -x[0])
    def trace(k):
        path = []
        while k is not None and k in pred:
            u, d, kind = pred[k]
            path.append((k, d, kind)); k = u
        return path[::-1]
    out = {'loops_unresolved': loops, 'modelled': dict(modelled), 'paths': [],
           'unmodelled_driven_pins': sorted('%s/%s' % k for k in unused if k in arr)}
    for tot, k, s, tag in res[:npaths]:
        out['paths'].append({'total_ns': tot, 'fmax_mhz': 1000.0 / tot if tot > 0 else None, 'end': k, 'setup': s,
                             'setup_tag': tag, 'arcs': trace(k) if 'internal' not in k[1] else []})
    resh = []
    for k, (s_, tag) in ends.items():
        a_ = arr.get(k) if k[0] in hard_insts else arr_h.get(k)
        if a_ is not None:
            resh.append((a_ + s_, k, s_, tag))
    resh.sort(key=lambda x: -x[0])
    def trace_h(k):
        path = []; mode = 'h' if k[0] not in hard_insts else 'a'
        while k is not None:
            if mode == 'h':
                if k not in pred_h: break
                e = pred_h[k]
                if e[0] is None:
                    path.append((k, e[1], e[2])); break
                u, d, kind, m = e
                path.append((k, d, kind)); k = u; mode = m
            else:
                if k not in pred: break
                u, d, kind = pred[k]
                path.append((k, d, kind)); k = u
        return path[::-1]
    if resh:
        tot, k, s_, tag = resh[0]
        out['hard_path'] = {'total_ns': tot, 'fmax_mhz': 1000.0 / tot, 'end': k, 'setup': s_, 'setup_tag': tag, 'arcs': trace_h(k)}
    out['worst_ns'] = res[0][0] if res else None
    out['fmax_mhz'] = 1000.0 / res[0][0] if res else None
    out['internal'] = internal
    return out

def PIPE_DEFAULT_SETUP(p, dsp):
    return dsp.PIPE['SETUP'].get(p, 0.1)

def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('sdf'); ap.add_argument('packed', nargs='?')
    ap.add_argument('--paths', type=int, default=3)
    ap.add_argument('--tdb', default=TDB_PATH)
    ap.add_argument('--json'); ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args()
    base = analyse(a.sdf, a.packed, a.tdb, correct=False, npaths=1)
    cor = analyse(a.sdf, a.packed, a.tdb, correct=True, npaths=a.paths)
    eq = analyse(a.sdf, a.packed, a.tdb, correct=True, npaths=1, equal_clock=True)
    eb = analyse(a.sdf, a.packed, a.tdb, correct=True, npaths=1, ebr_setups=True)
    eqeb = analyse(a.sdf, a.packed, a.tdb, correct=True, npaths=1, ebr_setups=True, equal_clock=True)
    cons = min(eb['fmax_mhz'], eqeb['fmax_mhz'])
    print('nextpnr model (recomputed from SDF): %.2f MHz (%.3f ns)' % (base['fmax_mhz'], base['worst_ns']))
    print('corrected (DSP/SPRAM arcs from timings_up5k.txt): %.2f MHz (%.3f ns)   cells re-modelled: %s   unresolved loop nodes: %d'
          % (cor['fmax_mhz'], cor['worst_ns'], cor['modelled'], cor['loops_unresolved']))
    print('corrected, DSP/SPRAM clock pins given the logic cells\' clock delay (no 0.93 ns hard-block skew): %.2f MHz (%.3f ns)' % (eq['fmax_mhz'], eq['worst_ns']))
    print('corrected, and block-RAM (EBR) setups from the database instead of nextpnr\'s 0.1 ns: %.2f MHz (%.3f ns)' % (eb['fmax_mhz'], eb['worst_ns']))
    print('CONSERVATIVE (DSP/SPRAM + EBR setups, the worse of the two hard-block clock models): %.2f MHz' % cons)
    hp = cor.get('hard_path')
    if hp:
        print('worst path through a DSP/SPRAM block (corrected): %.2f MHz (%.3f ns), ends at %s/%s' % (hp['fmax_mhz'], hp['total_ns'], hp['end'][0], hp['end'][1]))
    if cor['unmodelled_driven_pins']:
        print('WARNING: DSP pins driven by logic but feeding nothing in the model (check the configuration): %d, e.g. %s'
              % (len(cor['unmodelled_driven_pins']), cor['unmodelled_driven_pins'][:6]))
    if not a.quiet:
        for i, p in enumerate(cor['paths']):
            print('\n-- path %d: %.3f ns (%.2f MHz), ends at %s/%s (setup %.3f, %s)' % (i + 1, p['total_ns'], p['fmax_mhz'], p['end'][0], p['end'][1], p['setup'], p['setup_tag']))
            for (n, d, kind) in p['arcs']:
                print('   %+7.3f  %-60s %s' % (d, (n[0] + '/' + n[1])[-60:], kind))
        if hp:
            print('\n-- worst DSP/SPRAM path: %.3f ns (%.2f MHz)' % (hp['total_ns'], hp['fmax_mhz']))
            for (n, d, kind) in hp['arcs']:
                print('   %+7.3f  %-60s %s' % (d, (n[0] + '/' + n[1])[-60:], kind))
        intl = sorted(cor['internal'], key=lambda x: -x[2])[:4]
        if intl:
            print('\ninternal DSP register-to-register estimates (worst 4): ' + '; '.join('%s %s %.2f ns' % (i, d, v) for i, d, v in intl))
    if a.json:
        json.dump({'nextpnr_model_mhz': base['fmax_mhz'], 'corrected_mhz': cor['fmax_mhz'], 'corrected_equal_clock_mhz': eq['fmax_mhz'], 'corrected_ebr_mhz': eb['fmax_mhz'], 'conservative_mhz': cons, 'hard_path_mhz': (cor.get('hard_path') or {}).get('fmax_mhz'), 'worst_ns': cor['worst_ns'],
                   'paths': [{k: v for k, v in p.items() if k != 'arcs'} | {'arcs': [[n[0], n[1], d, kind] for n, d, kind in p['arcs']]} for p in cor['paths']],
                   'modelled': cor['modelled']}, open(a.json, 'w'), indent=1, default=str)

if __name__ == '__main__':
    main()
