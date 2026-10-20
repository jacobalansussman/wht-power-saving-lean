"""afc-data (wht45/afc, 20 Oct 2026): the hand-written instance modules of A-fc (#2979's words, T74 clock, the
re-solved fill u1-2979sep of round wht45/fullcut, FILL 7fc766c6 / PROFILE 5df79411; suffix 2979Fc), made from
x-tables' #2979 modules by a CHECKED renaming and number substitution: aA-data's mk_inst2979H.py, extended.

What is different from design A (aA-data): the A-fc fill changes the IDLE REMAINDERS, so PROFILE.json changes too.
So, beyond the head count (16,576, as in design A), the numbers that change are:
  * F1 (unhosted rank-48 bins of children 12/13, not heads): 0 (design A 11,614, #2979 12,104);
  * the remainder histogram, hence `sched2979Fc_hcost`'s right side (= the nonzero part of MANIFEST rem_hist) and the
    idle part of `unit2979Fc_arith`;
  * the unit profile `unit2979Fc` (= A-fc's PROFILE.json, checked: 666 calls x the local word + the remainders).
Only cut-INDEPENDENT modules are made here (x-tables' 7 instance modules, plus three small new ones below).  The half-cut
bridge modules of aA-data (Inst2979H, HeadsFirst2979H, HFTest2979H(W), R5Joint2979H, L2979H/EngineOfGeomW) are NOT
made: A-fc uses the FULL cut, whose bridge, R5 and adapter are afc-cut's.
What must NOT change is checked from the data before anything is written (the generator refuses, exit 2):
  * heads 16,576 = 8,288 + 8,288 (FOUND.md), F1 = 0, H + F1 = the rank-48 count of PROFILE.json, audit and heads ALL OK;
  * the residual sums (both providers) = x-tables' (the residual vector is fill-independent);
  * the dispatch module regenerated from DISPATCH/ORDER is byte-identical to x-tables' Data/D2979/Disp.lean;
  * the host table's call fields (time, provider, windows, residuals, ...) equal x-tables' table (fill-independence.json);
  * 666 x the local word (PROVIDERS.json, both providers) + the remainder histogram = PROFILE.json, exactly;
  * MANIFEST is of the audited table.
Renaming: every identifier / module component containing `2979` gets `Fc` right after `2979` (sched2979_fits ->
sched2979Fc_fits, S2979HB -> S2979FcHB, Data.H2979 -> Data.H2979Fc, DH2979 -> DH2979Fc ...), EXCEPT the generic or
fill-independent names: I2979 (and names after `I2979.`), Count2979C, HFW2979, D2979 / D2979W / D2979Net / D2979Geom /
Disp2979W (dispatch), lever2979.  After the rewrite every remaining old number and every unrenamed `2979` token outside
that list is a refusal.
New modules (templates in this file, numbers from the data):
  (`PB.unit2979Fc`, the bridge's copy of the same list, is afc-cut's: Work/Mixed/Bridge/Inst2979Fc.lean, afc-cut NEEDS.md
  section 1; this generator writes `out/unit2979Fc.txt` so the two lists can be compared as text);
  * Work/XHostT/Data/H2979Fc/NegF1.lean: negative control, F1 count 1 (must FAIL);
  * Work/XHostT/Inst/S2979FcHFW.lean: the heads-first certificate of the warp network on `sched2979Fc` (x-inst's
    `xnetGeomW_heads_first_auto` at the A-fc schedule, with the residual link as a hypothesis) and `hLlen2979Fc`.

usage: python3 -I -B mk_inst2979Fc.py OUT_LEAN_ROOT     (refuses to overwrite an existing file)"""
import sys, os, re, json, hashlib, gzip, collections

HERE = os.path.dirname(os.path.abspath(__file__))
M = os.path.dirname(HERE)
XT = os.environ['GEN_ROOT'] + '/wht45/lever2979/x-tables'
SRC = {  # source (relative to x-tables/lean) -> output path
    'Work/XHostT/Data/H2979/Heads.lean': 'Work/XHostT/Data/H2979Fc/Heads.lean',
    'Work/XHostT/Data/H2979/NegHeads.lean': 'Work/XHostT/Data/H2979Fc/NegHeads.lean',
    'Work/XHostT/Data/H2979/DispLink.lean': 'Work/XHostT/Data/H2979Fc/DispLink.lean',
    'Work/XHostT/Inst/S2979.lean': 'Work/XHostT/Inst/S2979Fc.lean',
    'Work/XHostT/Inst/S2979Facts.lean': 'Work/XHostT/Inst/S2979FcFacts.lean',
    'Work/XHostT/Inst/S2979Unit.lean': 'Work/XHostT/Inst/S2979FcUnit.lean',
    'Work/XHostT/Inst/S2979HB.lean': 'Work/XHostT/Inst/S2979FcHB.lean',
}
PIN = {  # = aA-data's pins (x-tables tree D, compiled)
    'Work/XHostT/Data/H2979/Heads.lean': '6988bed457053a66aa5618e11f58a4832ceb3712aeead3fdec27ce5df68678bb',
    'Work/XHostT/Data/H2979/NegHeads.lean': 'c3897f6b760dd26831255ae7f66f07c4843b7816d175214a5692fdf5bdfa6ca8',
    'Work/XHostT/Data/H2979/DispLink.lean': '5983446e1b4ae5c5e5027e1f6377c46563693d8cd4dd1e8056b87bc100d1608d',
    'Work/XHostT/Inst/S2979.lean': '58d1d0531535ff6bb73daed810f19c876e9e99ab762655c380e434a9f22197f5',
    'Work/XHostT/Inst/S2979Facts.lean': '28f0bc6873178d5b9761a0561c5f68ba7d408a8be1988f770ad222417d6e92ed',
    'Work/XHostT/Inst/S2979Unit.lean': 'ddec4f5459a32a1dd1905042677590811f5797e82e64e61d5a6596bf8793e339',
    'Work/XHostT/Inst/S2979HB.lean': '3666d3a92c28ebe1eb77cd85b387e42b2156d6a6aaeae330ad0a1e548c5a9822',
}
KEEP = {'I2979', 'Count2979C', 'HFW2979', 'D2979', 'D2979W', 'D2979Net', 'D2979Geom', 'Disp2979W', 'lever2979'}
SUF = 'Fc'


def die(msg):
    print('REFUSED:', msg)
    sys.exit(2)


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def rd(p):
    b = open(p, 'rb').read()
    return json.loads(gzip.decompress(b) if p.endswith('.gz') else b)


com = lambda x: '{:,}'.format(x)


def main():
    out = sys.argv[1]
    # ---- the data
    hs = rd(os.path.join(M, 'out/heads2979Fc-summary.json'))
    aud = rd(os.path.join(M, 'out/H2979Fc-audit.json'))
    man_new = rd(os.path.join(out, 'Work/XHostT/Data/H2979Fc/MANIFEST.json'))
    man_old = rd(os.path.join(XT, 'lean/Work/XHostT/Data/H2979/MANIFEST.json'))
    ind = rd(os.path.join(M, 'out/fill-independence.json'))
    hs_old = rd(os.path.join(XT, 'out/heads2979-summary.json'))
    prof = {int(a): b for a, b in rd(os.path.join(M, 'in/alloc2979Fc/PROFILE.json')).items()}
    provs = rd(os.path.join(M, 'in/alloc2979Fc/PROVIDERS.json'))
    if sha(os.path.join(M, 'in/alloc2979Fc/PROFILE.json')) != '5df7941120150492c64e62dc3362095eff3dfac9facff041ebb5c25b2e78f75e':
        die('PROFILE.json is not A-fc\'s (5df79411)')
    if not (hs['all_ok'] and aud['all_ok']):
        die('heads / audit not all ok')
    H, F1 = hs['heads'], hs['nonheads']
    if H != 16576 or hs['by_child'] != {'8': 8288, '10': 8288}:
        die('head count %d is not FOUND.md\'s 16,576 (8,288 + 8,288)' % H)
    if F1 != 0:
        die('A-fc has no non-head rank-48 bin (FOUND.md: 48: 16,576 heads only); got %d' % F1)
    if H + F1 != prof[48] or {int(a): b for a, b in aud['profile'].items()} != prof:
        die('H + F1 != the rank-48 count of PROFILE.json, or the table profile is not PROFILE.json')
    if man_new['res_sums'] != man_old['res_sums']:
        die('residual sums differ from x-tables\'')
    if man_new['sha256'] != aud['hosttab_content_sha256']:
        die('MANIFEST is not of the audited host table')
    for k in ('disp', 'time', 'order_k', 'windows', 'res', 'local_profile', 'child_of'):
        if ind['fields'][k] != 'same':
            die('host-table field %s depends on the fill' % k)
    if open(os.path.join(M, 'rerun/disp/Disp.lean'), 'rb').read() != open(os.path.join(XT, 'lean/Work/XHostT/Data/D2979/Disp.lean'), 'rb').read():
        die('the regenerated dispatch module differs from x-tables\'')
    # local word of one call (both providers equal) and the idle part
    lp = {int(a): b for a, b in provs['extra']['actual_paid_profile'].items()}
    if {int(a): b for a, b in provs['gauge44']['actual_paid_profile'].items()} != lp:
        die('the two providers\' local profiles differ')
    rem = {int(a): b for a, b in man_new['rem_hist'].items() if int(a) != 0 and b}
    chk = collections.Counter()
    for a, b in lp.items():
        chk[a] += 666 * b
    for a, b in rem.items():
        chk[a] += b
    if dict(chk) != prof:
        die('666 x the local word + remainders != PROFILE.json')
    rem_old = {int(a): b for a, b in man_old['rem_hist'].items() if int(a) != 0 and b}
    if rem_old != {1: 17310, 3: 4956, 5: 1728, 48: 28190}:
        die('x-tables\' remainder histogram is not the expected one')
    lp_old = [(1, 2445), (2, 706), (3, 232), (4, 108), (5, 18), (8, 8)]
    if sorted(lp.items()) != lp_old:
        die('the local word is not x-tables\' lp2979')
    unit_old = [(1, 1645680), (2, 470196), (3, 159468), (4, 71928), (5, 13716), (8, 5328), (48, 28190)]
    unit_new = sorted(prof.items())

    def lst(u):
        return '[' + ', '.join('(%d, %d)' % ab for ab in u) + ']'

    def hc_lean(r):
        return ' + '.join('%d * φ %d' % (b, a) for a, b in sorted(r.items()))

    def hc_prose(r):
        return ' + '.join('%s φ%d' % (com(b), a) for a, b in sorted(r.items()))

    def hc_prose2(r):
        return ' + '.join('%d φ%d' % (b, a) for a, b in sorted(r.items()))

    H0, F0 = hs_old['heads'], hs_old['nonheads']
    oc = hs_old['by_child']; on = hs_old['nonheads_by_child']; nc = hs['by_child']
    TXT = [  # exact text (old, new), applied first, before the number rules
        (hc_lean(rem_old), hc_lean(rem)),
        (hc_prose(rem_old), hc_prose(rem)),
        (hc_prose2(rem_old), hc_prose2(rem)),
        (lst(unit_old), lst(unit_new)),
        ('`allocation/PROFILE.json`, sha256 94a81974...', 'A-fc\'s `PROFILE.json`, sha256 5df79411...: the idle remainders of the re-solved fill'),
        ('`allocation/PROFILE.json`', 'A-fc\'s `PROFILE.json`'),
        ('the idle part of #2979\'s PROFILE.json', 'the idle part of A-fc\'s PROFILE.json'),
        ('#2979\'s whole unit profile', 'A-fc\'s whole unit profile'),
        ('#2979\'s unit profile', 'A-fc\'s unit profile'),
        ('#2979\'s profile', 'A-fc\'s profile'),
        ('#2979\'s unit price is its PROFILE.json', 'A-fc\'s unit price is its PROFILE.json'),
        ('(unhosted, not heads: %s + %s)' % (com(on['12']), com(on['13'])), '(unhosted, not heads: none in A-fc)'),
        ('the literal heads: %s of child 8 + %s of child 10' % (com(oc['8']), com(oc['10'])),
         'the literal heads: %s of child 8 + %s of child 10' % (com(nc['8']), com(nc['10']))),
        ('data from x-tables\'\nconversion of #2979\'s literal owner fill: mdumitrean, research/rank7-joint-endframes, on Jeremy Kothe\'s #2957 and Mason\nHale\'s #2917 T74 words',
         'data from afc-data\'s\nconversion of A-fc\'s fill (u1-2979sep, 7fc766c6) on #2979\'s words: mdumitrean, research/rank7-joint-endframes, on\nJeremy Kothe\'s #2957 and Mason Hale\'s #2917 T74 words'),
        ('The host table of hub PR #2979 as a Lean', 'The host table of A-fc (hub PR #2979\'s words, re-solved fill) as a Lean'),
        ('the cross-replica host schedule of #2979.', 'the cross-replica host schedule of A-fc (#2979\'s words).'),
        ('the cross-replica host schedule of hub PR #2979', 'the cross-replica host schedule of A-fc (hub PR #2979\'s words, re-solved fill)'),
        ('**the head count of #2979**', '**the head count of A-fc**'),
        ('**the price of the remaining climbs of #2979\'s host chains.**', '**the price of the remaining climbs of A-fc\'s host chains.**'),
        ('#2979\'s head bins', 'A-fc\'s head bins'),
        ('**the head bins of #2979**', '**the head bins of A-fc**'),
        ('per-kind data facts of #2979\'s host table', 'per-kind data facts of A-fc\'s host table'),
        ('The data facts of #2979\'s schedule', 'The data facts of A-fc\'s schedule (#2979\'s words)'),
    ]
    NUM = [  # (old, new), whole numbers only
        (str(H0 + 1), str(H + 1)), (str(H0), str(H)), (str(F0), str(F1)),
        (com(H0), com(H)), (com(F0), com(F1)),
    ]
    OLDNUM = [str(H0), str(H0 + 1), str(F0), com(H0), com(F0), com(oc['8']), com(oc['10']), com(on['12']), com(on['13']),
              '17310', '4956', '1728', '28190', '1645680', '159468', '13716', '17,310', '4,956', '1,728', '28,190',
              '94a81974']
    tok = re.compile(r'[A-Za-z_][A-Za-z0-9_]*2979[A-Za-z0-9_]*')

    def keep(t2, mo):
        w = mo.group(0)
        return (w in KEEP or t2[max(0, mo.start() - 6):mo.start()] == 'I2979.'
                or (w == 'Inst2979' and t2[max(0, mo.start() - 3):mo.start()] in ('HF.', 'HF/')))

    written, report = [], []
    for rel, dst in sorted(SRC.items()):
        src = os.path.join(XT, 'lean', rel)
        if sha(src) != PIN[rel]:
            die('source %s is not the compiled one (sha pin)' % rel)
        t = open(src).read()
        nsub = {}
        for a, b in TXT:
            k = t.count(a)
            if k:
                t = t.replace(a, b)
                nsub['txt:' + a[:30]] = k

        def ren(mo):
            w = mo.group(0)
            if keep(t, mo):
                return w
            if w.startswith('I2979') or ('2979' + SUF) in w or '2979H' in w:
                die('unexpected token %s in %s' % (w, rel))
            i = w.index('2979') + 4
            nsub[w] = nsub.get(w, 0) + 1
            return w[:i] + SUF + w[i:]
        t2 = tok.sub(ren, t)
        for a, b in NUM:
            t2, k = re.subn(r'(?<![0-9.,])' + re.escape(a) + r'(?![0-9])', b, t2)
            if k:
                nsub[a] = k
        for a in OLDNUM:
            if re.search(r'(?<![0-9.,])' + re.escape(a) + r'(?![0-9])', t2):
                die('old number %s left in %s' % (a, dst))
        for mo in tok.finditer(t2):
            w = mo.group(0)
            if ('2979' + SUF) in w or keep(t2, mo):
                continue
            die('unrenamed token %s in %s' % (w, dst))
        lines = t2.split('\n')
        k = 0
        while k < len(lines) and lines[k].startswith('import '):
            k += 1
        prov = ('-- afc-data (wht45/afc, 20 Oct 2026): GENERATED by tools/mk_inst2979Fc.py from x-tables/lean/%s (sha256 %s):\n'
                '-- A-fc = #2979\'s words with the re-solved fill u1-2979sep (FILL 7fc766c6, PROFILE 5df79411), suffix 2979Fc;'
                ' heads %d, F1 %d.' % (rel, sha(src)[:16], H, F1))
        lines.insert(k, prov)
        t2 = '\n'.join(lines)
        write(out, dst, t2, written)
        report.append(dict(src='x-tables/lean/%s' % rel, src_sha256=sha(src), dst=dst,
                           dst_sha256=sha(os.path.join(out, dst)), substitutions=nsub))
    # ---- new modules
    new = {}
    new['Work/XHostT/Data/H2979Fc/NegF1.lean'] = f'''import Work.XHostT.KHeadC
import Work.XHostT.Data.H2979Fc.Calls
/-! afc-data NEGATIVE control: the F1 count (unhosted rank-48 bins of children 12/13) one too high (A-fc has none);
this module must FAIL to compile.  GENERATED by tools/mk_inst2979Fc.py. -/
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace OAI.PowerSaving.XHostT.K.DH2979Fc
open OAI.PowerSaving.XHostT.K
theorem f1_bad : cntK rt 48 6 74 = {F1 + 1} := by decide +kernel
end OAI.PowerSaving.XHostT.K.DH2979Fc
'''
    new['Work/XHostT/Inst/S2979FcHFW.lean'] = f'''import Work.XHostT.Inst.S2979FcHB
import Work.XHostT.HF.Inst2979
import Work.XHostT.HFW2979.GeomInstW

/-!
# (key: heads-first, A-fc) The heads-first certificate of the #2979 warp network on A-fc's schedule

x-inst's warp heads-first certificate (`XHostTW.xnetGeomW_heads_first_auto`, Work/XHostT/HFW2979/GeomInstW.lean,
generic in the schedule) at A-fc's host schedule `sched2979Fc` (this layer): for ANY warp providers whose entrance
residuals are the host table's (`hres`; for the real #2979 provider words this is x-words' `hent` against
`K.DH2979Fc.resT`, the integration's `resIcW` lemma), the word `headWord (heads of headBins sched2979Fc) ++ body` is a
live kernel, `body` proper, with price `card Γ x unitPriceTW Ic sched2979Fc - |heads| x φ 48`, and there are exactly
`card Γ x {H}` head roles (`hLlen2979Fc`; {H} = `headBins2979Fc_card`).  The cut (afc-cut's full cut, `MX.cut_heads_full`)
is applied to this certificate downstream; nothing here depends on the cut.
GENERATED by afc-data's tools/mk_inst2979Fc.py (wht45/afc, 20 Oct 2026); prepared with assistance from Claude.  No `sorry`.
-/

set_option linter.unusedSectionVars false
set_option linter.unusedVariables false
namespace OAI.PowerSaving.XHostT.Inst
open Finset Binary Matrix RAM SS CB RF BR GF GX BC B401 XHost XHostR XHostD XHostT XHostTW
open OAI.PowerSaving.XHostT.K.DH2979Fc

noncomputable section

/-- the head-role count of A-fc: `card Γ` gauge copies of each of the {H} head bins. -/
theorem hLlen2979Fc (Γ : Type) [Fintype Γ] [DecidableEq Γ] :
    (headsL (Γ := Γ) (nb := 0) (headBins sched2979Fc)).length = Fintype.card Γ * {H} :=
  I2979.hLlen_card sched2979Fc (resK ct resT 74 408) sched2979Fc_fits resK2979Fc_pos sched2979Fc_heads

variable {{H C : Type}} [Fintype H] [DecidableEq H] [Fintype C] [DecidableEq C]
variable {{Ok : ∀ {{β : Type}} [Fintype β] [DecidableEq β], GWin H β → Prop}}

/-- **THE HEADS-FIRST CERTIFICATE OF A-fc** on the #2979 warp network (any providers with the table's residuals). -/
theorem heads_first2979Fc (I : GInvW H (Fin 112) (Fin 408) C Ok) (Ic : Fin 74 → Fin 9 → GInvW H (Fin 112) (Fin 408) C Ok)
    (hIc : ∀ ρ l, (Ic ρ l).tv = I.tv) (base : Fin 112 → OBase H) (l : H) (hbase : ∀ t, I.tv t = (base t).v l)
    (hH : Fintype.card H = 9) (hOk : ∀ (n : ℕ) (d : Orth (BC.LC (Fin 8) H l)), Ok (win401N l hH n d))
    (hres : ∀ ρ l' q, resIcW Ic ρ l' q = resK ct resT 74 408 ρ l' q) :
    ∃ body : GWord (BC.LC (Fin 8) H l) (BRoleT (Orth (BC.LC (Fin 8) H l)) (Fin 112) 0 C 74),
      body.Proper (Fintype.card (BC.LC (Fin 8) H l)) ∧
      (∀ φ, body.costR φ = (Fintype.card (Orth (BC.LC (Fin 8) H l)) : ℝ) * unitPriceTW Ic sched2979Fc φ
          - ((Fintype.card (Orth (BC.LC (Fin 8) H l)) * {H} : ℕ) : ℝ) * φ 48) ∧
      (headsL (Γ := Orth (BC.LC (Fin 8) H l)) (nb := 0) (headBins sched2979Fc)).length
          = Fintype.card (Orth (BC.LC (Fin 8) H l)) * {H} ∧
      GLiveKernel (Sum.inl : BLiveT (Orth (BC.LC (Fin 8) H l)) (Fin 112) 0 74
          → BRoleT (Orth (BC.LC (Fin 8) H l)) (Fin 112) 0 C 74)
        (MX.headWord Sum.inl (headsL (headBins sched2979Fc)) (fun _ => 48) (zGeom base l hH) (zGeom_indep base l hH)
          ++ body) := by
  have e : resIcW Ic = resK ct resT 74 408 := by
    funext ρ l' q
    exact hres ρ l' q
  have hfit : sched2979Fc.fits (resIcW Ic) := by
    rw [e]
    exact sched2979Fc_fits
  have hL := hLlen2979Fc (Orth (BC.LC (Fin 8) H l))
  obtain ⟨body, hP, hc, hk⟩ := xnetGeomW_heads_first_auto I Ic hIc base l hbase hH hOk sched2979Fc hfit
  refine ⟨body, hP, fun φ => ?_, hL, hk⟩
  rw [hc φ, hL]

end
end OAI.PowerSaving.XHostT.Inst

#print axioms OAI.PowerSaving.XHostT.Inst.hLlen2979Fc
#print axioms OAI.PowerSaving.XHostT.Inst.heads_first2979Fc
'''
    for dst, txt in sorted(new.items()):
        write(out, dst, txt, written)
        report.append(dict(src='template in mk_inst2979Fc.py', dst=dst, dst_sha256=sha(os.path.join(out, dst))))
    open(os.path.join(M, 'out/unit2979Fc.txt'), 'w').write(lst(unit_new) + '\n')
    json.dump(dict(generator_sha256=sha(os.path.abspath(__file__)), H=H, F1=F1, rem=rem, unit=unit_new, files=report),
              open(os.path.join(M, 'out/mk_inst2979Fc-report.json'), 'w'), indent=1)
    print('ok: %d modules written (H %d, F1 %d, remainders %s)' % (len(written), H, F1, rem))
    for r in report:
        print('  %-45s <- %-50s %s' % (r['dst'], r['src'], ', '.join('%s:%d' % kv for kv in sorted(r.get('substitutions', {}).items()))))


def write(out, dst, txt, written):
    p = os.path.join(out, dst)
    if os.path.exists(p):
        die('refusing to overwrite ' + p)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, 'w').write(txt)
    written.append(dst)


if __name__ == '__main__':
    main()
