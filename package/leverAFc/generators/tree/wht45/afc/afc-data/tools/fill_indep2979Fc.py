"""afc-data (wht45/afc, 20 Oct 2026): fill independence of the A-fc host table, as aA-data's check (STATUS 14:31,
out/fill-independence.json): every call field of the new table equals x-tables' #2979 table; only the role rows,
offsets and remainders may differ.  usage: python3 -I -B fill_indep2979Fc.py NEW.json.gz OLD.json.gz OUT.json"""
import sys, json, gzip, collections
a = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
b = json.loads(gzip.decompress(open(sys.argv[2], 'rb').read()))
fields = {k: ('same' if a[k] == b[k] else 'DIFFERS') for k in sorted(set(a) | set(b)) if k != 'source'}
src = {k: ('same' if a['source'][k] == b['source'][k] else 'DIFFERS') for k in a['source']}
ch = collections.Counter('%d->%d' % (x, y) for x, y in zip(b['hrem'], a['hrem']) if x != y)
reh = sum(1 for ra, rb in zip(a['hrows'], b['hrows']) for x, y in zip(ra, rb) if x != y)
out = dict(fields=fields, source=src, roles_rehosted=reh, hrem_changes=dict(sorted(ch.items())))
json.dump(out, open(sys.argv[3], 'w'), indent=1)
bad = [k for k in ('disp', 'time', 'order_k', 'windows', 'res', 'local_profile', 'child_of', 'T', 'R', 'NH', 'nb', 'kinds',
                   'ranks', 'kz', 'kk', 'hpair', 'providers', 'format') if fields[k] != 'same']
print(json.dumps(out)); print('FILL-INDEPENDENT CALL FIELDS' if not bad else 'DEPENDS: %s' % bad)
sys.exit(0 if not bad else 2)
