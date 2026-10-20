"""x-tables (wht45/lever2979, 19 Oct 2026): our lowering's residuals of #2979's two provider words, in g-heads'
residuals format (data/residuals-3945.json: per provider `residuals` (R widths, helper order), `exits`,
`helper_origin`), read from x-gate's lowered words (gcert/1w, made by p-lower's transport from the hub's BORROW words).
Also records each word's block histogram (our lowering, warps free) for the converter's local-profile check.
usage: python3 -I mk_resid2979.py XGATE_OUT_DIR OUT.json"""
import sys, json, gzip, hashlib
xd, outf = sys.argv[1:3]
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
out = {'producer': 'x-tables mk_resid2979.py', 'source': {}}
for w, name in (('base', 'gauge44'), ('extra', 'extra')):
    gp = f'{xd}/{w}.gcert1w.json.gz'; lp = f'{xd}/lower-{w}.json'
    g = json.loads(gzip.decompress(open(gp, 'rb').read())); lo = json.load(open(lp))
    assert g['provider'] == name and g['R'] == 408 and len(g['residuals']) == 408 and len(g['helper_origin']) == 408
    assert lo['frames_and_price']['hist_equals_hub_local_profile'] is True
    out[name] = dict(residuals=g['residuals'], exits=g['exits'], helper_origin=g['helper_origin'],
                     block_hist=lo['frames_and_price']['hist'])
    out['source'][f'{w}.gcert1w.json.gz'] = sha(gp); out['source'][f'lower-{w}.json'] = sha(lp)
json.dump(out, open(outf, 'w'))
print('wrote', outf, {k: out['source'][k][:16] for k in out['source']})
