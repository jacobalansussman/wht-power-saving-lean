"""Regenerate mathematical FILL/PROFILE payloads into a new explicit directory."""
from pathlib import Path
from collections import Counter
import argparse,json,gzip,hashlib
ROOT=Path(__file__).resolve().parent
def need(ok,why):
    if not ok:raise ValueError(why)
def read(path):
    b=path.read_bytes();return json.loads(gzip.decompress(b) if path.suffix=='.gz' else b)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def rank(rows):
    piv={}
    for x in rows:
        while x:
            b=x.bit_length()
            if b in piv:x^=piv[b]
            else:piv[b]=x;break
    return len(piv)
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir',required=True);args=parser.parse_args()
    pins=read(ROOT/'SOURCES.json')['files']
    for rel,digest in pins.items():need(sha(ROOT/rel)==digest,'source hash mismatch: '+rel)
    I=ROOT/'inputs';order=read(I/'ORDER.json');dispatch={(s,r):pid for s,r,pid in read(I/'DISPATCH.json')}
    fill=read(I/'old-FILL.json.gz');widths={};profiles={};full=[1<<i for i in range(9)]
    for tag,pid in [('base','gauge44'),('extra','extra')]:
        word=read(I/f'WORD-{tag}.json.gz');origin=read(I/f'old-NATIVE-{tag}.json.gz')['helper_origin']
        widths[pid]=[rank(word['helper_final'].get(str(r),full))-rank(word['start_extra'].get(str(r),[])) for r in origin]
        profiles[pid]={int(r):n for r,n in word['report']['local_profile'].items()}
    caps={c['child']:c['rank'] for c in read(I/'TEMPLATE.json.gz')['idle_children']};used={}
    for j,row in enumerate(fill['rows']):
        t,h,o,c,p,z=row;need((t,h)==divmod(j,408),'helper identity');key=(o,c,p);prefix=used.get(key,0)
        width=widths[dispatch[tuple(order[t])]][h];need(prefix+width<=caps[c],'owner capacity')
        row[5]=prefix;used[key]=prefix+width
    profile=Counter()
    for o in range(360):
        for c,cap in caps.items():
            for p in range(112):
                idle=cap-used.get((o,c,p),0)
                if idle:profile[idle]+=1
    for call in order:profile.update(profiles[dispatch[tuple(call)]])
    need(fill==read(ROOT/'FILL.json.gz'),'regenerated role payload differs')
    need(dict(profile)=={int(r):n for r,n in read(ROOT/'PROFILE.json').items()},'regenerated profile differs')
    out=Path(args.output_dir);need(not out.exists(),'output directory must be fresh');out.mkdir(parents=True)
    payload=json.dumps(fill,separators=(',',':')).encode();(out/'FILL.json.gz').write_bytes(gzip.compress(payload,compresslevel=6,mtime=0))
    (out/'PROFILE.json').write_text(json.dumps(dict(sorted(profile.items())),indent=2),encoding='utf8')
    print(json.dumps({'status':'PASS_REGENERATED_LITERAL_FILL','roles':len(fill['rows']),'fill_payload_sha256':hashlib.sha256(payload).hexdigest(),'profile':dict(sorted(profile.items()))},sort_keys=True),flush=True)
if __name__=='__main__':main()
