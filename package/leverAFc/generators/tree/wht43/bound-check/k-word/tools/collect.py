"""Collect every distinct canonical complex word (CANONICAL-WORD.json[.gz]) from directories and tarballs
(read in memory, nothing extracted), keyed by sha256 of the decompressed JSON. Saves each distinct word as
words/<sha16>.json.gz and an index words/INDEX.json {sha16: {origins:[...], v,h,R, profile, pairs, ...}}.
Read-only towards the inputs; untrusted data is only parsed as JSON."""
import sys, os, json, gzip, hashlib, tarfile, io
from collections import Counter
OUT=sys.argv[1]; roots=sys.argv[2:]
os.makedirs(OUT,exist_ok=True)
idxp=os.path.join(OUT,'INDEX.json')
index=json.load(open(idxp)) if os.path.exists(idxp) else {}
def summar(c):
    v,h,R=c['v'],c['h'],c['R']; F=c['frames']; st,fi=c['start'],c['final']
    pairs=Counter((len(F[st[r]]),len(F[fi[r]])) for r in range(2*v,2*v+R))
    return dict(v=v,h=h,R=R,profile={int(k):n for k,n in c['local_profile'].items()},
        pairs=sorted([list(k)+[n] for k,n in pairs.items()]),
        src=sorted([list(k)+[n] for k,n in Counter((len(F[st[r]]),len(F[fi[r]])) for r in range(v)).items()]),
        tgt=sorted([list(k)+[n] for k,n in Counter((len(F[st[r]]),len(F[fi[r]])) for r in range(v,2*v)).items()]),
        cst=sum(len(F[f]) for f in c['center_frames']))
def add(raw, origin):
    try:
        data=gzip.decompress(raw) if raw[:2]==b'\x1f\x8b' else raw
        c=json.loads(data)
    except Exception as e:
        print('skip',origin,e); return
    if not isinstance(c,dict) or 'local_profile' not in c: return
    k=hashlib.sha256(data).hexdigest()[:16]
    if k not in index:
        try: s=summar(c)
        except Exception as e: print('bad',origin,e); return
        index[k]=dict(origins=[],**s)
        open(os.path.join(OUT,k+'.json.gz'),'wb').write(gzip.compress(data,mtime=0))
    if origin not in index[k]['origins']: index[k]['origins'].append(origin)
def scan_tar(p, label):
    try:
        with tarfile.open(p,'r:*') as t:
            for m in t:
                if m.isfile() and os.path.basename(m.name).startswith('CANONICAL-WORD.json'):
                    add(t.extractfile(m).read(), label+'::'+m.name)
                elif m.isfile() and m.name.endswith('.tar.gz') and m.size<50_000_000:
                    inner=t.extractfile(m).read()
                    tmp=io.BytesIO(inner)
                    try:
                        with tarfile.open(fileobj=tmp,mode='r:*') as t2:
                            for m2 in t2:
                                if m2.isfile() and os.path.basename(m2.name).startswith('CANONICAL-WORD.json'):
                                    add(t2.extractfile(m2).read(), label+'::'+m.name+'::'+m2.name)
                    except Exception as e: print('inner skip',m.name,e)
    except Exception as e: print('tar skip',p,e)
for root in roots:
    if os.path.isfile(root):
        if root.endswith(('.tar.gz','.tgz','.tar')): scan_tar(root, root)
        else: add(open(root,'rb').read(), root)
        continue
    for d,_,fs in os.walk(root):
        for f in fs:
            p=os.path.join(d,f)
            if f.startswith('CANONICAL-WORD.json'): add(open(p,'rb').read(), p)
            elif f.endswith(('.tar.gz','.tgz')) and os.path.getsize(p)<60_000_000: scan_tar(p,p)
json.dump(index,open(idxp,'w'),indent=0)
print('distinct words',len(index))
