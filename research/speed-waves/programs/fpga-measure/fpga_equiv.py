#!/usr/bin/env python3
"""fpga_equiv.py -- is a new design equal to a reference?  Proof (SymbiYosys) or simulation (Icarus
Verilog), and a deliberately broken copy to show that the check can fail.  Standard library only.

    python3 -B fpga_equiv.py prove  --new NEW.v [--new MORE.v] --ref REF.v [--new-top A] [--ref-top B]
                                    [--lat-new 0] [--lat-ref 0] [--engine bitwuzla|z3|yices|abc-pdr]
                                    [--timeout 300] [--work DIR] [--bg]
    python3 -B fpga_equiv.py sim    --new NEW.v --ref REF.v [...tops, latencies...] [--vectors 100000]
                                    [--field 8] [--rand-seed 1] [--work DIR] [--heavy | --bg]
    python3 -B fpga_equiv.py mutate --new NEW.v -o BROKEN.v [--line N]
    python3 -B fpga_equiv.py selftest [--work DIR]      (the recipe run on the circuit study's 8-point designs;
                                                         set WHT41_RTL to their folder, see below)

THE TWO MODULES must have the same ports (names, directions, widths; read by Yosys).  Every input
except the clock (--clk, default clk/clock/...) is driven with the same value in both; ALL outputs
(data and valid alike) are compared.

LATENCY.  --lat-new / --lat-ref: clocks from an input to its output (0 = combinational).  The output
of the design with the smaller latency is delayed by the difference, and the comparison starts once
both pipelines are full.  So a pipelined design can be checked against a combinational reference.

prove   Writes a miter (both modules, the delay registers, a saturating cycle counter and
        "assert (out_new == out_ref_delayed)" once the counter reaches max latency) and runs
        SymbiYosys in mode prove (bounded check from the initial state + k-induction) with depth
        max latency + 2 and the chosen engine.  PASS = equal for every input sequence.  k-induction
        is complete for feed-forward pipelines (no register feeds back) of depth <= max latency;
        a design with feedback (accumulators, counters, folded designs with a state machine) can
        make the induction fail with an unreachable state: then the BMC part is still a check of
        the first steps, and the report says "induction failed".
        WARNING (lesson of the earlier circuit study): two DIFFERENT networks of adders compared bit by bit (e.g. the
        19-adder circuit against the textbook) time out in every engine, even at 8 points.  Prove
        through a helper instead: (A) design = the same lines at full width (bitwuzla; this checks
        the exact widths), (B) full-width lines = the definition (z3: it normalises linear sums).
        gen_rtl.py --variant wide writes such helpers; selftest shows the recipe.
sim     Writes the input vectors (corner cases first: all 0, all 1, every field at its signed
        minimum and maximum, alternating bits; then a third of "extreme" vectors whose fields are
        drawn from {min, max, -1, 0, 1}; the rest uniform random) and a test bench that compares the
        two outputs every cycle.  PASS = no mismatch.  --field W: the width of one data field
        (inputs whose width is a multiple of W are cut in W-bit fields for the corner cases).
mutate  Writes a broken copy: the first (or --line N) "a + b" on an assignment line becomes
        "a - b" (or "-" becomes "+").  Run prove and sim on it: both must FAIL.

Exit status: 0 PASS, 1 FAIL, 2 TIMEOUT/UNKNOWN, 3 error.  Last line of output: one JSON object.

RUNNING.  As in fpga_measure.py: when executable wrappers heavy.sh and bg.sh sit next to this file,
Yosys, SymbiYosys and the solvers run through heavy.sh and the simulation through bg.sh (--heavy:
everything through heavy.sh; --bg: everything through bg.sh); without them every tool runs directly.
At most 4 simulations run at once on the machine (lock files in the system temporary folder).
SELFTEST.  It reads the 8-point designs of the circuit study published in the repository
wht-fewer-operations (tools/circuits/fpga/rtl): set WHT41_RTL to that folder.  It also needs the
helper formal/helpers/8_ours_wide_w8.v there (written by that folder's gen_rtl.py --variant wide).
"""
import argparse
import contextlib
import fcntl
import json
import os
import random
import re
import shutil
import signal
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.dont_write_bytecode = True
import fpga_measure as FM  # noqa: E402

BIN = FM.BIN
WHT41 = os.environ.get('WHT41_RTL', os.path.join(HERE, 'wht41-rtl'))   # the circuit study's FPGA rtl folder (selftest only)
LOCKS = os.path.join(tempfile.gettempdir(), 'fpga_equiv_locks')
os.makedirs(LOCKS, exist_ok=True)
SIM_SLOTS = 4      # simulations at once on the machine 


def runner(a, sim=False):
    """heavy.sh for Yosys on whole designs and for the solvers (the default); bg.sh for the simulation
    (unless --heavy) and for everything with --bg"""
    return FM.runner(not a.bg and (a.heavy or not sim))


def run(r, cmd, log, timeout, cwd=None):
    """cmd under the runner r; the seconds and the timeout count from the tool's own start (see FM.sh);
    a timeout kills every descendant (sby starts its engines in process groups of their own)"""
    rc, dt = FM.sh(r, cmd, log, timeout, cwd=cwd)
    return rc, round(dt, 2)


def guard(res, rc):
    """res, with a note when the memory guard stopped the tool (exit code 75, 97 or 98)"""
    if rc in FM.GUARD:
        res['guard'] = FM.GUARD[rc]
    return res


@contextlib.contextmanager
def sim_slot():
    """hold one of the SIM_SLOTS lock files LOCKS/sim-N.lock while simulating; waits for a free one"""
    said = False
    while True:
        for i in range(1, SIM_SLOTS + 1):
            f = open(os.path.join(LOCKS, 'sim-%d.lock' % i), 'a')
            try:
                fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError:
                f.close()
                continue
            try:
                yield i
            finally:
                f.close()
            return
        if not said:
            print('[fpga_equiv] waiting for one of the %d simulation slots' % SIM_SLOTS, file=sys.stderr, flush=True)
            said = True
        time.sleep(2)


def ports_of(files, top, wd, a):
    ns = argparse.Namespace(top=top, param=[], incdir=a.incdir, define=a.define, sv=a.sv, heavy=not a.bg)
    mod, ports, seq = FM.probe(ns, files, wd)
    return mod, ports, seq


def setup(a, wd):
    os.makedirs(wd, exist_ok=True)
    nf = [os.path.abspath(f) for f in a.new]
    rf = [os.path.abspath(f) for f in a.ref]
    mn, pn, _ = ports_of(nf, a.new_top, os.path.join(wd, 'probe_new'), a)
    mr, pr, _ = ports_of(rf, a.ref_top, os.path.join(wd, 'probe_ref'), a)
    if pn != pr:
        raise SystemExit('the two modules have different ports:\n new %s\n ref %s' % (pn, pr))
    if mn == mr:
        raise SystemExit('both modules are called %s: rename one (e.g. the broken copy from mutate is renamed)' % mn)
    names = [p for p, d, w in pn if d == 'input' and w == 1]
    clk = a.clk or next((n for n in FM.CLOCK_NAMES if n in names), '__noclock__')
    ins = [(p, w) for p, d, w in pn if d == 'input' and p != clk]
    outs = [(p, w) for p, d, w in pn if d == 'output']
    return nf, rf, mn, mr, clk, ins, outs, any(p == clk for p, _, _ in pn)


def inst(mod, name, clk, has_clk, ins, outs, outbus):
    c, lo = [], 0
    if has_clk:
        c.append('.%s(clk)' % clk)
    for p, w in ins:
        c.append('.%s(hin[%d:%d])' % (p, lo + w - 1, lo))
        lo += w
    lo = 0
    for p, w in outs:
        c.append('.%s(%s[%d:%d])' % (p, outbus, lo + w - 1, lo))
        lo += w
    return '    %s %s (%s);' % (mod, name, ', '.join(c))


def delay_lines(bus, n, NO, prefix):
    """n delay registers behind bus; returns (lines, name of the delayed bus)"""
    if n == 0:
        return [], bus
    L = ['    reg [%d:0] %s1%s;' % (NO - 1, prefix, ''.join(', %s%d' % (prefix, i) for i in range(2, n + 1)))]
    body = '%s1 <= %s;' % (prefix, bus) + ''.join(' %s%d <= %s%d;' % (prefix, i, prefix, i - 1) for i in range(2, n + 1))
    L.append('    always @(posedge clk) begin %s end' % body)
    return L, '%s%d' % (prefix, n)


# ------------------------------------------------------------------------------------------------
def prove(a):
    wd = os.path.abspath(a.work)
    nf, rf, mn, mr, clk, ins, outs, has_clk = setup(a, wd)
    NI, NO = sum(w for _, w in ins), sum(w for _, w in outs)
    LM = max(a.lat_new, a.lat_ref)
    s = ['module miter (input wire clk, input wire [%d:0] hin);' % (NI - 1),
         '    wire [%d:0] on, orf;' % (NO - 1),
         inst(mn, 'u_new', clk, has_clk, ins, outs, 'on'),
         inst(mr, 'u_ref', clk, has_clk, ins, outs, 'orf')]
    dl, bn = delay_lines('on', LM - a.lat_new, NO, 'dn')
    s += dl
    dl, br = delay_lines('orf', LM - a.lat_ref, NO, 'dr')
    s += dl
    s += ['    reg [15:0] cnt = 16\'d0;',
          '    always @(posedge clk) if (cnt != 16\'hffff) cnt <= cnt + 16\'d1;',
          '    always @(*) if (cnt >= 16\'d%d) assert (%s == %s);' % (LM, bn, br),
          'endmodule', '']
    open(os.path.join(wd, 'miter.v'), 'w').write('\n'.join(s))
    files = nf + rf
    if LM == 0 and a.engine != 'abc-pdr':
        return prove_comb(a, wd, files, mn, mr, ins, outs, clk, has_clk, NI, NO)
    eng = {'bitwuzla': 'smtbmc bitwuzla', 'z3': 'smtbmc z3', 'yices': 'smtbmc yices', 'abc-pdr': 'abc pdr'}[a.engine]
    rd = '\n'.join('read_verilog %s %s' % ('-sv' if a.sv else '', os.path.basename(f)) for f in files)
    sby = ('[options]\nmode prove\ndepth %d\ntimeout %d\n\n[engines]\n%s\n\n[script]\nread_verilog -formal miter.v\n%s\n'
           'prep -top miter\nflatten\nasync2sync\n\n[files]\nmiter.v\n%s\n') % (
        LM + 2, a.timeout, eng, rd, '\n'.join(files))
    # (two files of the same base name would clash in sby's [files]; copy them under unique names if so)
    if len({os.path.basename(f) for f in files}) != len(files):
        raise SystemExit('two input files share a base name; rename one')
    open(os.path.join(wd, 'eq.sby'), 'w').write(sby)
    shutil.rmtree(os.path.join(wd, 'eq'), ignore_errors=True)
    cmd = [os.path.join(BIN, 'sby'), '-f', 'eq.sby']
    rc, dt = run(runner(a), cmd, os.path.join(wd, 'sby.log'), a.timeout + 30, cwd=wd)
    log = open(os.path.join(wd, 'sby.log')).read()
    st = 'UNKNOWN'
    if rc == 'timeout' or re.search(r'DONE \(TIMEOUT', log):
        st = 'TIMEOUT'
    elif re.search(r'DONE \(PASS', log):
        st = 'PASS'
    elif re.search(r'DONE \(FAIL', log):
        st = 'FAIL'
        if re.search(r'induction', log, re.I) and not re.search(r'basecase.*(failed|FAIL)', log, re.I) and re.search(r'Assert failed in miter.*\n.*induction|induction.*failed', log, re.I):
            st = 'FAIL (induction failed: may be an unreachable state; see log)'
    res = {'mode': 'prove', 'result': st, 'seconds': dt, 'engine': a.engine, 'new': mn, 'ref': mr,
           'lat_new': a.lat_new, 'lat_ref': a.lat_ref, 'depth': LM + 2, 'NI': NI, 'NO': NO, 'log': os.path.join(wd, 'sby.log')}
    m = re.findall(r'(Assert failed in miter.*|counterexample trace.*|Writing trace to VCD file: .*)', log)
    if m:
        res['counterexample'] = m[:3]
    print(json.dumps(guard(res, rc)), flush=True)
    return {'PASS': 0, 'FAIL': 1}.get(st.split(' ')[0], 2)


def prove_comb(a, wd, files, mn, mr, ins, outs, clk, has_clk, NI, NO):
    """both combinational: one time step, no counter (the light path of wht41/fpga/rtl/formal/run_formal.py):
    yosys writes the miter as SMT-LIB, yosys-smtbmc -t 1 asks the solver for an input that breaks it"""
    s = ['module miter (input wire clk, input wire [%d:0] hin);' % (NI - 1),
         '    wire [%d:0] on, orf;' % (NO - 1),
         inst(mn, 'u_new', clk, has_clk, ins, outs, 'on'),
         inst(mr, 'u_ref', clk, has_clk, ins, outs, 'orf'),
         '    always @(*) assert (on == orf);', 'endmodule', '']
    open(os.path.join(wd, 'miter_comb.v'), 'w').write('\n'.join(s))
    ys = ['read_verilog -formal miter_comb.v'] + ['read_verilog %s %s' % ('-sv' if a.sv else '', f) for f in files] + [
          'hierarchy -top miter', 'proc', 'flatten', 'async2sync', 'opt_clean', 'write_smt2 -wires miter.smt2']
    open(os.path.join(wd, 'miter.ys'), 'w').write('\n'.join(ys) + '\n')
    rc, dt1 = run(runner(a), [os.path.join(BIN, 'yosys'), '-q', '-l', 'yosys.log', 'miter.ys'],
                  os.path.join(wd, 'yosys.out'), 600, cwd=wd)
    if rc != 0:
        print(json.dumps(guard({'mode': 'prove', 'result': 'ERROR', 'log': os.path.join(wd, 'yosys.log')}, rc)))
        return 3
    rc, dt2 = run(runner(a), [os.path.join(BIN, 'yosys-smtbmc'), '-s', a.engine, '-t', '1', '--noprogress',
                              '--dump-vcd', 'cex.vcd', 'miter.smt2'], os.path.join(wd, 'smtbmc.log'), a.timeout, cwd=wd)
    log = open(os.path.join(wd, 'smtbmc.log')).read()
    st = 'TIMEOUT' if rc == 'timeout' else ('PASS' if 'Status: PASSED' in log else ('FAIL' if 'Status: FAILED' in log else 'UNKNOWN'))
    res = {'mode': 'prove', 'path': 'combinational (smtbmc -t 1)', 'result': st, 'seconds': round(dt1 + dt2, 2),
           'engine': a.engine, 'new': mn, 'ref': mr, 'NI': NI, 'NO': NO, 'log': os.path.join(wd, 'smtbmc.log')}
    if st == 'FAIL':
        res['counterexample'] = os.path.join(wd, 'cex.vcd')
    print(json.dumps(guard(res, rc)), flush=True)
    return {'PASS': 0, 'FAIL': 1}.get(st, 2)


# ------------------------------------------------------------------------------------------------
def vectors(ins, n, field, seed):
    rnd = random.Random(seed)
    V = []

    def fields(p, w):
        if field and w % field == 0 and w > field:
            return [field] * (w // field)
        return [w]

    def pack(choose):
        x, lo = 0, 0
        for p, w in ins:
            for fw in fields(p, w):
                v = choose(fw) & ((1 << fw) - 1)
                x |= v << lo
                lo += fw
        return x

    smin = lambda fw: 1 << (fw - 1) if fw > 1 else 1
    smax = lambda fw: (1 << (fw - 1)) - 1 if fw > 1 else 0
    allw = lambda fw: (1 << fw) - 1
    corners = [lambda fw: 0, allw, smin, smax, lambda fw: int('01' * fw, 2) & allw(fw),
               lambda fw: int('10' * fw, 2) & allw(fw), lambda fw: 1]
    for c in corners:
        V.append(pack(c))
    for i in range(n - len(V)):
        if i % 3 == 0:
            V.append(pack(lambda fw: rnd.choice((smin(fw), smax(fw), allw(fw), 0, 1))))
        else:
            V.append(pack(lambda fw: rnd.getrandbits(fw)))
    return V[:n]


def sim(a):
    wd = os.path.abspath(a.work)
    nf, rf, mn, mr, clk, ins, outs, has_clk = setup(a, wd)
    NI, NO = sum(w for _, w in ins), sum(w for _, w in outs)
    LM = max(a.lat_new, a.lat_ref)
    V = vectors(ins, a.vectors, a.field, a.rand_seed)
    hexw = (NI + 3) // 4
    with open(os.path.join(wd, 'vec.hex'), 'w') as f:
        for v in V:
            f.write('%0*x\n' % (hexw, v))
    dn, dr = LM - a.lat_new, LM - a.lat_ref
    tb = ['`timescale 1ns/1ps', 'module tb;',
          '    reg clk = 1\'b0;', '    reg [%d:0] hin;' % (NI - 1), '    wire [%d:0] on, orf;' % (NO - 1),
          '    reg [%d:0] mem [0:%d];' % (NI - 1, len(V) - 1),
          inst(mn, 'u_new', clk, has_clk, ins, outs, 'on'),
          inst(mr, 'u_ref', clk, has_clk, ins, outs, 'orf'),
          '    reg [%d:0] hn [0:%d];' % (NO - 1, max(dn, 1)), '    reg [%d:0] hr [0:%d];' % (NO - 1, max(dr, 1)),
          '    integer t, k, bad;',
          '    reg [%d:0] cn, cr;' % (NO - 1),
          '    initial begin',
          '        $readmemh("vec.hex", mem);',
          '        bad = 0;',
          '        for (t = 0; t < %d; t = t + 1) begin' % (len(V) + LM),
          '            hin = mem[t %% %d];' % len(V),
          '            #1;',
          # shift the histories, newest at index 0
          '            for (k = %d; k > 0; k = k - 1) begin hn[k] = hn[k-1]; end' % max(dn, 1),
          '            for (k = %d; k > 0; k = k - 1) begin hr[k] = hr[k-1]; end' % max(dr, 1),
          '            hn[0] = on; hr[0] = orf;',
          '            cn = hn[%d]; cr = hr[%d];' % (dn, dr),
          '            if (t >= %d && cn !== cr) begin' % LM,
          '                bad = bad + 1;',
          '                if (bad <= 5) $display("MISMATCH cycle %0d in=%h new=%h ref=%h", t, mem[(t - ' + str(a.lat_new) + ') % ' + str(len(V)) + '], cn, cr);',
          '            end',
          '            #4 clk = 1\'b1; #5 clk = 1\'b0;',
          '        end',
          '        if (bad == 0) $display("RESULT PASS vectors=%0d", ' + str(len(V)) + ');',
          '        else $display("RESULT FAIL mismatches=%0d of %0d", bad, ' + str(len(V)) + ');',
          '        $finish;',
          '    end', 'endmodule', '']
    open(os.path.join(wd, 'tb.v'), 'w').write('\n'.join(tb))
    with sim_slot():
        rc, dt1 = run(runner(a, sim=True), [os.path.join(BIN, 'iverilog'), '-g2012', '-o', 'tb.vvp', '-s', 'tb', 'tb.v'] + nf + rf,
                      os.path.join(wd, 'iverilog.log'), a.timeout, cwd=wd)
        if rc != 0:
            print(json.dumps(guard({'mode': 'sim', 'result': 'ERROR', 'log': os.path.join(wd, 'iverilog.log')}, rc)))
            return 3
        rc, dt2 = run(runner(a, sim=True), [os.path.join(BIN, 'vvp'), '-n', 'tb.vvp'], os.path.join(wd, 'vvp.log'), a.timeout, cwd=wd)
    log = open(os.path.join(wd, 'vvp.log')).read()
    st = 'TIMEOUT' if rc == 'timeout' else ('PASS' if 'RESULT PASS' in log else ('FAIL' if 'RESULT FAIL' in log else 'UNKNOWN'))
    res = {'mode': 'sim', 'result': st, 'seconds': round(dt1 + dt2, 2), 'vectors': len(V), 'new': mn, 'ref': mr,
           'lat_new': a.lat_new, 'lat_ref': a.lat_ref, 'field': a.field, 'NI': NI, 'NO': NO, 'log': os.path.join(wd, 'vvp.log')}
    m = re.findall(r'(MISMATCH.*|RESULT.*)', log)
    if st != 'PASS':
        res['detail'] = [x[:200] for x in m[:3]]
    print(json.dumps(guard(res, rc)), flush=True)
    return {'PASS': 0, 'FAIL': 1}.get(st, 2)


def mutate(a):
    src = open(a.new[0]).read().split('\n')
    done = None
    for i, line in enumerate(src):
        code = line.split('//')[0]
        if '=' not in code or re.match(r'\s*(input|output|parameter|localparam|for\b)', code):
            continue
        rhs_at = code.index('=') + 1
        rhs = code[rhs_at:]
        if a.line and i + 1 != a.line:
            continue
        for op, rep in ((' + ', ' - '), (' - ', ' + ')):
            if op in rhs:
                j = rhs_at + rhs.index(op)
                src[i] = line[:j] + rep + line[j + len(op):]
                done = (i + 1, line.strip(), src[i].strip())
                break
        if done:
            break
    if not done:
        raise SystemExit('no "a + b" or "a - b" found on an assignment line')
    txt = '\n'.join(src)
    # rename every module so that the copy can sit beside the original
    txt = re.sub(r'(^\s*module\s+)(\w+)', lambda m: m.group(1) + 'broken_' + m.group(2), txt, flags=re.M)
    open(a.o, 'w').write(txt)
    print(json.dumps({'mode': 'mutate', 'line': done[0], 'was': done[1], 'now': done[2], 'out': os.path.abspath(a.o),
                      'note': 'modules renamed broken_<name>'}))
    return 0


def selftest(a):
    """the recipe on the wht41 8-point designs, 8-bit inputs"""
    V = os.path.join(WHT41, 'v')
    Hd = os.path.join(WHT41, 'formal', 'helpers')
    w = os.path.abspath(a.work)
    os.makedirs(w, exist_ok=True)
    out = []

    def go(name, fn, **kw):
        ns = argparse.Namespace(new=[], ref=[], new_top=None, ref_top=None, lat_new=0, lat_ref=0, engine='bitwuzla',
                                timeout=a.timeout, work=os.path.join(w, name), clk=None, vectors=a.vectors, field=8,
                                rand_seed=1, heavy=a.heavy, bg=a.bg, incdir=[], define=[], sv=False, o=None, line=None)
        for k, v in kw.items():
            setattr(ns, k, v)
        print('# ' + name, flush=True)
        rc = fn(ns)
        out.append((name, rc))
        return rc

    broken = os.path.join(w, 'broken_8_ours_c_w8.v')
    go('mutate', mutate, new=[os.path.join(V, '8_ours_c_w8.v')], o=broken)
    go('sim_ours_vs_textbook', sim, new=[os.path.join(V, '8_ours_c_w8.v')], ref=[os.path.join(V, '8_textbook_c_w8.v')])
    go('sim_broken_vs_textbook', sim, new=[broken], ref=[os.path.join(V, '8_textbook_c_w8.v')])
    go('sim_pipelined6_vs_textbook', sim, new=[os.path.join(V, '8_ours_p_w8.v')], ref=[os.path.join(V, '8_textbook_c_w8.v')], lat_new=6)
    go('sim_pipelined_wrong_latency', sim, new=[os.path.join(V, '8_ours_p_w8.v')], ref=[os.path.join(V, '8_textbook_c_w8.v')], lat_new=5)
    go('proveA_ours_vs_wide', prove, new=[os.path.join(V, '8_ours_c_w8.v')], ref=[os.path.join(Hd, '8_ours_wide_w8.v')])
    go('proveB_wide_vs_definition', prove, new=[os.path.join(Hd, '8_ours_wide_w8.v')], ref=[os.path.join(V, '8_def_c_w8.v')], engine='z3')
    go('prove_broken_vs_wide', prove, new=[broken], ref=[os.path.join(Hd, '8_ours_wide_w8.v')])
    go('prove_pipelined6_vs_comb', prove, new=[os.path.join(V, '8_ours_p_w8.v')], ref=[os.path.join(V, '8_ours_c_w8.v')], lat_new=6)
    go('prove_pipelined_wrong_latency', prove, new=[os.path.join(V, '8_ours_p_w8.v')], ref=[os.path.join(V, '8_ours_c_w8.v')], lat_new=5)
    want = {'mutate': 0, 'sim_ours_vs_textbook': 0, 'sim_broken_vs_textbook': 1, 'sim_pipelined6_vs_textbook': 0,
            'sim_pipelined_wrong_latency': 1, 'proveA_ours_vs_wide': 0, 'proveB_wide_vs_definition': 0,
            'prove_broken_vs_wide': 1, 'prove_pipelined6_vs_comb': 0, 'prove_pipelined_wrong_latency': 1}
    ok = all(want[n] == rc for n, rc in out)
    print(json.dumps({'mode': 'selftest', 'result': 'PASS' if ok else 'FAIL',
                      'cases': {n: {'rc': rc, 'expected': want[n]} for n, rc in out}}))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('mode', choices=['prove', 'sim', 'mutate', 'selftest'])
    ap.add_argument('--new', action='append', default=[])
    ap.add_argument('--ref', action='append', default=[])
    ap.add_argument('--new-top', default=None)
    ap.add_argument('--ref-top', default=None)
    ap.add_argument('--lat-new', type=int, default=0)
    ap.add_argument('--lat-ref', type=int, default=0)
    ap.add_argument('--clk', default=None)
    ap.add_argument('--engine', default='bitwuzla', choices=['bitwuzla', 'z3', 'yices', 'abc-pdr'])
    ap.add_argument('--timeout', type=int, default=300)
    ap.add_argument('--vectors', type=int, default=100000)
    ap.add_argument('--field', type=int, default=None)
    ap.add_argument('--rand-seed', type=int, default=1)
    ap.add_argument('--work', default=os.path.join(os.getcwd(), 'equiv_work'))
    ap.add_argument('--heavy', action='store_true', help='put the simulation through heavy.sh too (proofs and Yosys go there by default)')
    ap.add_argument('--bg', action='store_true', help='everything through bg.sh: only for tiny designs; still capped by the memory guard')
    ap.add_argument('--sv', action='store_true')
    ap.add_argument('--define', action='append', default=[])
    ap.add_argument('--incdir', action='append', default=[])
    ap.add_argument('-o', default=None)
    ap.add_argument('--line', type=int, default=None)
    a = ap.parse_args()
    if a.bg and a.heavy:
        ap.error('--bg and --heavy together')
    for s in (signal.SIGTERM, signal.SIGHUP, signal.SIGINT):
        signal.signal(s, FM.stop)
    if BIN:
        os.environ['PATH'] = BIN + ':' + os.environ.get('PATH', '')
    if a.mode == 'mutate':
        if not a.o or not a.new:
            ap.error('mutate needs --new FILE and -o OUT')
        return mutate(a)
    if a.mode == 'selftest':
        return selftest(a)
    if not a.new or not a.ref:
        ap.error('--new and --ref are needed')
    return prove(a) if a.mode == 'prove' else sim(a)


if __name__ == '__main__':
    sys.exit(main())
