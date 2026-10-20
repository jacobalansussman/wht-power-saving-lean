# p-alt: write a minimal 'old-NATIVE' stub holding only helper_origin = range(224, 704) minus the word's removed helpers
# (control: this equals #3945's old-NATIVE helper_origin exactly). heads_fill.py reads nothing else from it.
import sys, json, gzip
w = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
rem = set(w['removed']); assert len(rem) == 72 and w['report']['R'] == 480
ho = [h for h in range(224, 704) if h not in rem]; assert len(ho) == 408
open(sys.argv[2], 'wb').write(gzip.compress(json.dumps({'helper_origin': ho, 'stub': 'p-alt'}).encode()))
print(sys.argv[2], len(ho), 'start_extra' in w, 'helper_final' in w)
