"""aA-data (wht45/allheads, 19 Oct 2026): the hand-written instance modules of design A (#2979, every head kept;
suffix 2979H), made from x-tables' and x-inst's #2979 modules by a CHECKED renaming, nothing else.

Only the FILL changed (reopt FOUND.md); words, ORDER, DISPATCH, PROVIDERS, TEMPLATE, PLAN and PROFILE are #2979's.
So the instance numbers that change are exactly: the head count H (16,086 -> heads2979H-summary.json), the other
rank-48 bins F1 (12,104 -> 28,190 - H), the per-child splits (prose), and the DFT numerator N of R5 (17556473 ->
rate-2979H.json), with the seam exponent 1 - (100 N + 1)/10^12 and N + 1 (prose).  What must NOT change is checked
from the data before anything is written (the generator refuses, exit 2):
  * the remainder histogram of the new host check (MANIFEST rem_hist) = x-tables' (so `sched2979_hcost`'s numbers
    and the unit profile `unit2979` = PROFILE.json carry over);
  * the residual sums (both providers) = x-tables';
  * the dispatch module regenerated from DISPATCH/ORDER is byte-identical to x-tables' Data/D2979/Disp.lean;
  * the host table's call fields (time, provider, windows, residuals) equal x-tables' table (fill-independence.json);
  * H + F1 = the rank-48 count of PROFILE.json, H = 16,576 = reopt FOUND.md's count.
Renaming: every identifier / module component containing `2979` gets `H` right after `2979` (sched2979_fits ->
sched2979H_fits, HFTest2979W -> HFTest2979HW, S2979HB -> S2979HHB, Data.H2979 -> Data.H2979H, L2979 -> L2979H ...),
EXCEPT the generic or fill-independent names that are reused unchanged:
  I2979 (x-inst's generic table-side head-first namespace), HF.Inst2979 / HF/Inst2979 (its module), Count2979C,
  HFW2979 (the warp heads-first folder), D2979 / D2979W / D2979Net / D2979Geom / Disp2979W (dispatch: fill-independent),
  lever2979 (folder names in comments), #2979 (the hub PR).
After the rewrite every remaining old number and every unrenamed `2979` token outside that list is a refusal.

usage: python3 -I -B mk_inst2979H.py OUT_LEAN_ROOT     (refuses to overwrite an existing file)"""
import sys, os, re, json, hashlib, gzip

HERE = os.path.dirname(os.path.abspath(__file__))
M = os.path.dirname(HERE)
XT = os.environ['GEN_ROOT'] + '/wht45/lever2979/x-tables'
XI = os.environ['GEN_ROOT'] + '/wht45/lever2979/x-inst'
SRC = {  # source (relative to x-tables/lean or x-inst/lean) -> (sha256, output path)
    (XT, 'Work/XHostT/Data/H2979/Heads.lean'): 'Work/XHostT/Data/H2979H/Heads.lean',
    (XT, 'Work/XHostT/Data/H2979/NegHeads.lean'): 'Work/XHostT/Data/H2979H/NegHeads.lean',
    (XT, 'Work/XHostT/Data/H2979/DispLink.lean'): 'Work/XHostT/Data/H2979H/DispLink.lean',
    (XT, 'Work/XHostT/Inst/S2979.lean'): 'Work/XHostT/Inst/S2979H.lean',
    (XT, 'Work/XHostT/Inst/S2979Facts.lean'): 'Work/XHostT/Inst/S2979HFacts.lean',
    (XT, 'Work/XHostT/Inst/S2979Unit.lean'): 'Work/XHostT/Inst/S2979HUnit.lean',
    (XT, 'Work/XHostT/Inst/S2979HB.lean'): 'Work/XHostT/Inst/S2979HHB.lean',
    (XI, 'Work/Mixed/R5Joint2979.lean'): 'Work/Mixed/R5Joint2979H.lean',
    (XI, 'Work/Mixed/Bridge/Inst2979.lean'): 'Work/Mixed/Bridge/Inst2979H.lean',
    (XI, 'Work/Mixed/Bridge/HeadsFirst2979.lean'): 'Work/Mixed/Bridge/HeadsFirst2979H.lean',
    (XI, 'Work/Mixed/Bridge/HFTest2979.lean'): 'Work/Mixed/Bridge/HFTest2979H.lean',
    (XI, 'Work/Mixed/Bridge/HFTest2979W.lean'): 'Work/Mixed/Bridge/HFTest2979HW.lean',
    (XI, 'Work/L2979/EngineOfGeomW.lean'): 'Work/L2979H/EngineOfGeomW.lean',
}
PIN = {  # the sources as compiled by x-tables (tree D, = out/x-tables-src.mac.sha256) and x-inst (T2 on C1, checked)
    'Work/XHostT/Data/H2979/Heads.lean': '6988bed457053a66aa5618e11f58a4832ceb3712aeead3fdec27ce5df68678bb',
    'Work/XHostT/Data/H2979/NegHeads.lean': 'c3897f6b760dd26831255ae7f66f07c4843b7816d175214a5692fdf5bdfa6ca8',
    'Work/XHostT/Data/H2979/DispLink.lean': '5983446e1b4ae5c5e5027e1f6377c46563693d8cd4dd1e8056b87bc100d1608d',
    'Work/XHostT/Inst/S2979.lean': '58d1d0531535ff6bb73daed810f19c876e9e99ab762655c380e434a9f22197f5',
    'Work/XHostT/Inst/S2979Facts.lean': '28f0bc6873178d5b9761a0561c5f68ba7d408a8be1988f770ad222417d6e92ed',
    'Work/XHostT/Inst/S2979Unit.lean': 'ddec4f5459a32a1dd1905042677590811f5797e82e64e61d5a6596bf8793e339',
    'Work/XHostT/Inst/S2979HB.lean': '3666d3a92c28ebe1eb77cd85b387e42b2156d6a6aaeae330ad0a1e548c5a9822',
    'Work/Mixed/R5Joint2979.lean': '1777c3fb02a28509803e84941b0997ae183a1919b82fec724eeb0df0f8f471f4',
    'Work/Mixed/Bridge/Inst2979.lean': '29ab8fcc02cf9862b8f5b2aabcd045cb25b767c69597b589d299cacb925643e8',
    'Work/Mixed/Bridge/HeadsFirst2979.lean': '413faeeeda84b0ee844b2f211a3613c5f7482b7b65539b268b19ba001ee3de83',
    'Work/Mixed/Bridge/HFTest2979.lean': 'f2a7338fc72707ba531912ffe5d113cd0202c9b8659df12f6cd7ae50cf227579',
    'Work/Mixed/Bridge/HFTest2979W.lean': '0414db094c5ebde7de158a813c2f69b6c93af43a0d43610f6946a26aaaa9a0cf',
    'Work/L2979/EngineOfGeomW.lean': '3e7f83d0ad1ba11ea962fa62395cbeed1cac54ca4a4eb570d340e86e8146c944',
}
KEEP = {'I2979', 'Count2979C', 'HFW2979', 'D2979', 'D2979W', 'D2979Net', 'D2979Geom', 'Disp2979W', 'lever2979'}


def die(msg):
    print('REFUSED:', msg)
    sys.exit(2)


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def rd(p):
    b = open(p, 'rb').read()
    return json.loads(gzip.decompress(b) if p.endswith('.gz') else b)


def main():
    out = sys.argv[1]
    # ---- the data: what changes and what must not
    hs = rd(os.path.join(M, 'out/heads2979H-summary.json'))
    aud = rd(os.path.join(M, 'out/H2979H-audit.json'))
    rate = rd(os.path.join(M, 'out/rate-2979H.json'))
    man_new = rd(os.path.join(out, 'Work/XHostT/Data/H2979H/MANIFEST.json'))
    man_old = rd(os.path.join(XT, 'lean/Work/XHostT/Data/H2979/MANIFEST.json'))
    ind = rd(os.path.join(M, 'out/fill-independence.json'))
    hs_old = rd(os.path.join(XT, 'out/heads2979-summary.json'))
    if not (hs['all_ok'] and aud['all_ok']):
        die('heads / audit not all ok')
    H, F1 = hs['heads'], hs['nonheads']
    if H != 16576 or hs['by_child'] != {'8': 8288, '10': 8288}:
        die('head count %d is not reopt FOUND.md\'s 16,576 (8,288 + 8,288)' % H)
    if H + F1 != aud['profile']['48'] or aud['profile'] != aud['profile_hub']:
        die('H + F1 != the rank-48 count of PROFILE.json')
    if man_new['rem_hist'] != man_old['rem_hist'] or man_new['res_sums'] != man_old['res_sums']:
        die('remainder histogram or residual sums differ from x-tables\': hcost / unit profile would change')
    if man_new['sha256'] != aud['hosttab_content_sha256']:
        die('MANIFEST is not of the audited host table')
    for k in ('disp', 'time', 'order_k', 'windows', 'res', 'local_profile', 'child_of'):
        if ind['fields'][k] != 'same':
            die('host-table field %s depends on the fill' % k)
    if open(os.path.join(M, 'rerun/disp/Disp.lean'), 'rb').read() != open(os.path.join(XT, 'lean/Work/XHostT/Data/D2979/Disp.lean'), 'rb').read():
        die('the regenerated dispatch module differs from x-tables\'')
    N = rate['dft_numerator']
    if rate['H'] != H or rate['children'].get('24') != H or rate['children'].get('48') != F1 or rate['refused'] != N + 1:
        die('rate-2979H.json is not of this head data')
    N1s = rate['sharp_numerator']
    H0, F0 = hs_old['heads'], hs_old['nonheads']
    N0 = 17556473
    oc = hs_old['by_child']; on = hs_old['nonheads_by_child']; nc = hs['by_child']; nn = hs['nonheads_by_child']
    com = lambda x: '{:,}'.format(x)
    NUM = [  # (old, new), applied as whole numbers (no digit before or after)
        (str(100 * N0 + 1), str(100 * N + 1)), (str(100 * N0 + 2), str(100 * N + 2)),
        (str(N0 + 1), str(N1s)), (str(N0), str(N)),
        ('1.7556473e-3', ('%.7e' % (N / 1e10)).replace('e-03', 'e-3')),
        (str(H0 + 1), str(H + 1)), (str(H0), str(H)), (str(F0), str(F1)),
        (com(H0), com(H)), (com(F0), com(F1)),
        ('%s of child 8 + %s of child 10' % (com(oc['8']), com(oc['10'])), '%s of child 8 + %s of child 10' % (com(nc['8']), com(nc['10']))),
        ('%s + %s' % (com(on['12']), com(on['13'])), '%s + %s' % (com(nn['12']), com(nn['13']))),
    ]
    OLDNUM = [str(100 * N0 + 1), str(100 * N0 + 2), str(N0), str(N0 + 1), '1.7556473', str(H0), str(H0 + 1), str(F0),
              com(H0), com(F0), com(oc['8']), com(oc['10']), com(on['12']), com(on['13'])]
    tok = re.compile(r'[A-Za-z_][A-Za-z0-9_]*2979[A-Za-z0-9_]*')
    written, report = [], []
    for (base, rel), dst in sorted(SRC.items()):
        src = os.path.join(base, 'lean', rel)
        if sha(src) != PIN[rel]:
            die('source %s is not the compiled one (sha pin)' % rel)
        t = open(src).read()
        nsub = {}

        def ren(mo):
            w = mo.group(0)
            if w in KEEP:
                return w
            # names of the generic namespace I2979 (x-inst's HF/Inst2979: card_live2979, ...) keep theirs
            if t[max(0, mo.start() - 6):mo.start()] == 'I2979.':
                return w
            # x-inst's generic module Work.XHostT.HF.Inst2979 keeps its name
            if w == 'Inst2979' and t[max(0, mo.start() - 3):mo.start()] in ('HF.', 'HF/'):
                return w
            if w.startswith('I2979') or '2979H' in w:
                die('unexpected token %s in %s' % (w, rel))
            i = w.index('2979') + 4
            nsub[w] = nsub.get(w, 0) + 1
            return w[:i] + 'H' + w[i:]
        t2 = tok.sub(ren, t)
        for a, b in NUM:
            t2, k = re.subn(r'(?<![0-9.,])' + re.escape(a) + r'(?![0-9])', b, t2)
            if k:
                nsub[a] = k
        # leftovers
        for a in OLDNUM:
            if re.search(r'(?<![0-9.,])' + re.escape(a) + r'(?![0-9])', t2):
                die('old number %s left in %s' % (a, dst))
        for mo in tok.finditer(t2):
            w = mo.group(0)
            if '2979H' in w or w in KEEP or (w == 'Inst2979' and t2[max(0, mo.start() - 3):mo.start()] in ('HF.', 'HF/')) \
                    or t2[max(0, mo.start() - 6):mo.start()] == 'I2979.':
                continue
            die('unrenamed token %s in %s' % (w, dst))
        # provenance line after the imports
        lines = t2.split('\n')
        k = 0
        while k < len(lines) and lines[k].startswith('import '):
            k += 1
        prov = ('-- aA-data (wht45/allheads, 19 Oct 2026): GENERATED by tools/mk_inst2979H.py from %s/lean/%s (sha256 %s):\n'
                '-- design A = #2979 with every head kept (FILL 764395fc), instance suffix 2979H; heads %d, F1 %d, N %d.'
                % (os.path.basename(base), rel, sha(src)[:16], H, F1, N))
        lines.insert(k, prov)
        t2 = '\n'.join(lines)
        p = os.path.join(out, dst)
        if os.path.exists(p):
            die('refusing to overwrite ' + p)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        open(p, 'w').write(t2)
        written.append(dst)
        report.append(dict(src='%s/lean/%s' % (os.path.basename(base), rel), src_sha256=sha(src), dst=dst,
                           dst_sha256=sha(p), substitutions=nsub))
    json.dump(dict(generator_sha256=sha(os.path.abspath(__file__)), H=H, F1=F1, N=N, sharp=N1s, files=report),
              open(os.path.join(M, 'out/mk_inst2979H-report.json'), 'w'), indent=1)
    print('ok: %d modules written (H %d, F1 %d, N %d)' % (len(written), H, F1, N))
    for r in report:
        print('  %-45s <- %-50s %s' % (r['dst'], r['src'], ', '.join('%s:%d' % kv for kv in sorted(r['substitutions'].items()))))


if __name__ == '__main__':
    main()
