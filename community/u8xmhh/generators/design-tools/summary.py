"""fs-price: one JSON with every figure of this run (DFT numerators over 10^10; c = term_le grid point used by the
generators and the Lean route, q = term_le4)."""
import json, os
F = os.path.join(os.environ.get('F', '.'), 'out/')  # as run: the agents' fs-price folder
def J(n):
    p = F + n
    return json.load(open(p)) if os.path.exists(p) else None
ctl = J('controls-price.json')
rows = []
def add(label, status, dft_c, dft_q=None, **kw):
    rows.append(dict(label=label, status=status, dft_c=dft_c, saving_c='%.7e' % (dft_c / 1e10), dft_q=dft_q, **kw))
add('Shea control (#588 atm3 + 112 source copies, dedicated banks)', 'Lean-checked (Shea); reproduced', ctl['shea']['dft'], W=ctl['shea']['W'])
add('C588M control (#588 + mutable sources, annealed table)', 'Lean-checked (m-gen/m-build); reproduced', ctl['c588m']['dft'], W=ctl['c588m']['W'])
add('U8XM control (U8 + exit frames + mutable, 5 centre copies, 3771 banks)', 'Lean-checked (k-combine); reproduced', ctl['u8xm']['dft'], W=ctl['u8xm']['W'])
for tag, label, st in (('h588_c0_samerun.repacked2.json', 'C588MHH control (#588 + mutable + same-run hosting)', 'Lean-checked (h-L3-build); reproduced, tables byte-identical'),
                       ('h869_c0_samerun.repacked2.json', 'U8 + exits, NO copies, same-run hosting (h-L3-design control)', 'Python price; reproduced; needs Cc generalised'),
                       ('h869_c0_samerun.padsa.json', 'U8 + exits, NO copies, same-run hosting, pads consolidated (reference)', 'Python price; needs Cc generalised (200-400 Lean lines)'),
                       ('hu8xm_sr.repacked2.json', 'U8XM + same-run hosting, host_repack2 banks', 'Python price; superseded by the pad-consolidated table'),
                       ('hu8xm_sr.padsa.json', 'U8XMHH = U8XM + same-run hosting, pads consolidated (THE TARGET)', 'Python price + Lean data generated; hostSchedK checked in Lean on the Mac'),
                       ('hu8xm_x_d1.repacked2.json', 'U8XM + cross-run hosting, deferred finals (reference)', 'Python price; no Lean format (cross-run tables)'),
                       ('hu8xm_x_d1.padsa.json', 'U8XM + cross-run hosting, deferred finals, pads consolidated (reference)', 'Python price; no Lean format'),
                       ('hu8xm_x_d0.repacked2.json', 'U8XM + cross-run hosting, strict order (reference)', 'Python price; no Lean format'),
                       ('hu8xm_x_d0.padsa.json', 'U8XM + cross-run hosting, strict order, pads consolidated (reference)', 'Python price; no Lean format')):
    r = J(tag)
    if r:
        add(label, st, r['dft'], r.get('dft_quartic'), W=r['W'], banks=r['banks'], pads=r['pads'], hosted=r['hosted'], unhosted=r['unhosted'],
            defer=r['defer'], samerun=r.get('samerun'), file='out/' + tag)
shea = rows[0]['dft_c']
for r in rows:
    r['vs_shea_pct'] = round(100 * (r['dft_c'] / shea - 1), 2)
gen = J('U8XMHH.json')
json.dump(dict(date='2026-10-13', agent='fs-price', rows=rows, generator_summary=gen), open(F + 'PRICES.json', 'w'), indent=1)
for r in rows:
    print('%-82s %s  q=%s  W=%s  %+.2f%%  [%s]' % (r['label'][:82], r['saving_c'], r['dft_q'], r.get('W'), r['vs_shea_pct'], r['status']))
