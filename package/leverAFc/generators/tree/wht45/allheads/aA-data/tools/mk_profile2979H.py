"""aA-data (wht45/allheads, 19 Oct 2026): the half-cut R5 profile of design A (#2979, every head kept), read from the
host-table audit (full profile = PROFILE.json) and the head list (literal heads), not from literals:
children = full profile with the H literal heads moved from rank 48 to rank 24 (the other rank-48 bins stay 48),
W = 74 x 896 arrays per gauge copy, kc = 49, m = 73 (x-inst's profile2979-palt.json format).
usage: python3 -I -B mk_profile2979H.py H2979H-audit.json heads2979H-summary.json OUT.json"""
import sys, json
aud, hs, out = sys.argv[1:4]
a = json.load(open(aud)); h = json.load(open(hs))
assert a['all_ok'] and h['all_ok'] and a['profile'] == a['profile_hub']
prof = {int(k): v for k, v in a['profile'].items()}
H = h['heads']; non = h['nonheads']
assert prof[48] == H + non, 'rank-48 bins = heads + non-heads'
ch = dict(prof); ch[48] = non; ch[24] = ch.get(24, 0) + H
W = 74 * 896
p = {'source': 'aA-data mk_profile2979H.py: %s (profile = PROFILE.json) + %s (HEAD-literal %d = %s), half cut: heads -> rank 24, other idle-48 stay 48'
     % (aud.split('/')[-1], hs.split('/')[-1], H, ' + '.join(str(v) for v in h['by_child'].values())),
     'children': {str(k): v for k, v in sorted(ch.items())}, 'W': W, 'H': H, 'kc': 49, 'm': 73}
json.dump(p, open(out, 'w')); print(json.dumps(p))
