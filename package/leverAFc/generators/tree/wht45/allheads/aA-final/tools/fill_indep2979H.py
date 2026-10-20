"""aA-final (wht45/allheads, 19 Oct 2026): the fill-independence record of design A as a program.  aA-data made
out/fill-independence.json with a short Python snippet (its STATUS, 14:31); this is that comparison written out: every
field of the new host table (design A's FILL) against x-tables' table (#2979's published FILL), every source pin, the
number of roles whose host row changed, and the remainder changes.  mk_inst2979H.py reads the result and refuses
unless the call fields (disp, time, order_k, windows, res, local_profile, child_of) are 'same'.
usage: python3 -I fill_indep2979H.py OLD_HOSTTAB.json.gz NEW_HOSTTAB.json.gz OUT.json"""
import collections, gzip, json, sys
a, b = (json.load(gzip.open(f)) for f in sys.argv[1:3])
out = {
    'fields': {k: ('same' if a[k] == b[k] else 'DIFFERS') for k in sorted(a) if k != 'source'},
    'source': {k: ('same' if a['source'][k] == b['source'][k] else 'DIFFERS') for k in a['source']},
    'roles_rehosted': sum(1 for ra, rb in zip(a['hrows'], b['hrows']) for x, y in zip(ra, rb) if x != y),
    'hrem_changes': dict(sorted(collections.Counter((x, y) for x, y in zip(a['hrem'], b['hrem']) if x != y).items())),
}
out['hrem_changes'] = {f'{x}->{y}': n for (x, y), n in out['hrem_changes'].items()}
assert set(a) == set(b) and len(a['hrows']) == len(b['hrows'])
with open(sys.argv[3], 'w') as f:
    json.dump(out, f, indent=1); f.write('\n')
print('fill independence:', {k: v for k, v in out['fields'].items() if v != 'same'}, 'differ; sources differing:',
      [k for k, v in out['source'].items() if v != 'same'], '; roles rehosted', out['roles_rehosted'])
