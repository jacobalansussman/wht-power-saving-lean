# p-alt: decide whether a hub word file carries warps / twisted adds / zero buffers. Reads JSON data only.
import sys, json, gzip, hashlib, re
def load(p):
    b = open(p, 'rb').read()
    h = hashlib.sha256(b).hexdigest()
    raw = gzip.decompress(b) if b[:2] == b'\x1f\x8b' else b
    return h, raw, json.loads(raw)
def ident(cols):
    # G given as columns of a 0/1 matrix (ints as bit masks or lists); identity test where possible
    try:
        return all((c == (1 << i)) if isinstance(c, int) else (list(c) == [1 if j == i else 0 for j in range(len(c))]) for i, c in enumerate(cols))
    except Exception:
        return None
for p in sys.argv[1:]:
    h, raw, d = load(p)
    s = raw.decode('utf-8', 'replace')
    out = dict(file=p.split('/dl/')[-1], sha256=h, kind=d.get('kind'), base=d.get('base'))
    w = d.get('warps')
    out['warps'] = None if w is None else len(w)
    ic = d.get('initial_corrections')
    out['initial_corrections'] = None if ic is None else len(ic)
    zb = d.get('zero_buffers')
    out['zero_buffers'] = None if zb is None else len(zb)
    out['warp_edges'] = None if d.get('warp_edges') is None else len(d['warp_edges'])
    out['borrow'] = d.get('borrow')
    out['R'] = d.get('R', (d.get('report') or {}).get('R'))
    out['removed'] = None if d.get('removed') is None else len(d['removed'])
    out['str_tadd'] = s.count('"tadd"') + s.count("'tadd'")
    out['str_TADD'] = len(re.findall(r'"(?:T?ADD_?TW|TWADD|tw_add|twisted[^"]*)"', s))
    out['str_Z'] = s.count('["Z",')
    out['str_warp'] = s.lower().count('warp')
    # identity warps?
    if isinstance(w, list) and w:
        nid = 0; nchk = 0
        for x in w:
            cols = x.get('G_columns') if isinstance(x, dict) else None
            if cols is not None:
                nchk += 1; r = ident(cols); nid += 1 if r else 0
        out['warps_G_identity'] = f'{nid}/{nchk}' if nchk else 'n/a (word form)'
        out['warp0'] = json.dumps(w[0])[:160]
    print(json.dumps(out))
