"""Independent literal owner audit for enlarged protected-terminal providers.

Window/family semantics follow the retained PR3734 owner audit (and its credited
PR3638 predicates). Widths are read directly from the final canonical words.
"""
from pathlib import Path
from collections import Counter,defaultdict
import json,gzip,hashlib,time,sys
ROOT=Path(__file__).resolve().parent;OUT=ROOT
OLD=ROOT/'inputs';CLOCK=OLD
def need(ok,why):
    if not ok:raise ValueError(why)
def read(path):
    raw=path.read_bytes();return json.loads(gzip.decompress(raw) if path.suffix=='.gz' else raw)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def dimension(rows):
    pivots={}
    for vector in rows:
        while vector:
            leading=vector.bit_length()
            if leading not in pivots:pivots[leading]=vector;break
            vector^=pivots[leading]
    return len(pivots)
def canonical_widths(word,origin):
    full=[1<<i for i in range(9)]
    return [dimension(word['helper_final'].get(str(r),full))-dimension(word['start_extra'].get(str(r),[])) for r in origin]
def main():
    begin=time.monotonic();T=360;N=9*T
    pins=read(ROOT/'SOURCES.json')['files']
    expected={'FILL.json.gz','PROFILE.json','MANIFEST.json','OWNER-LIFT-PROOF.txt','check_gain.py'}
    expected.update('inputs/'+name for name in ['old-FILL.json.gz','ORDER.json','DISPATCH.json','TEMPLATE.json.gz','PLAN.json.gz','old-PROFILE.json','old-NATIVE-base.json.gz','old-NATIVE-extra.json.gz','old-WORD-base.json.gz','old-WORD-extra.json.gz','WORD-base.json.gz','WORD-extra.json.gz'])
    need(set(pins)==expected,'complete source hash guard')
    for rel,digest in pins.items():need(sha(ROOT/rel)==digest,'source hash mismatch: '+rel)
    manifest=read(OUT/'MANIFEST.json');fill=read(OUT/'FILL.json.gz');oldfill=read(OLD/'old-FILL.json.gz')
    need(sha(OUT/'FILL.json.gz')==manifest['fill_sha256'],'fill source hash')
    need(hashlib.sha256(gzip.decompress((OUT/'FILL.json.gz').read_bytes())).hexdigest()==manifest['fill_payload_sha256'],'fill payload hash')
    order=read(OLD/'ORDER.json');dispatch=read(OLD/'DISPATCH.json')
    need(len(order)==N and set(map(tuple,order))=={(s,r) for s in range(9) for r in range(T)},'call coverage')
    clock={tuple(call):t for t,call in enumerate(order)}
    need(all(all(clock[s,r]<clock[s+1,r] for s in range(8)) for r in range(T)),'per-replica source chronology')
    provider={(s,r):tag for s,r,tag in dispatch}
    need(len(provider)==len(dispatch)==N and set(provider)==set(clock),'dispatch coverage')
    template=read(OLD/'TEMPLATE.json.gz');plan=read(OLD/'PLAN.json.gz')
    children={row['child']:row['rank'] for row in template['idle_children']};K=max(children)+1
    pair={row['stage']:row['pair'] for row in template['recipe'] if row['kind']=='local'}
    windows={}
    for owner in range(T):
        touched={};tick=None
        for event in template['execution_schedule']:
            kind=event['kind']
            if kind in ('bank_padding_cover_sweep','translate_every_live_family'):continue
            if kind=='local_cover_sweep':tick=clock[event['stage'],owner];touched[pair[event['stage']]]=tick
            elif kind=='bridge':
                need(tick is not None,'bridge source');touched[event['src']]=tick;touched[event['dst']]=tick
            elif kind=='idle_cover_sweep':
                need(tick is not None,'idle source')
                for child in event['child_indices']:
                    if child in children:
                        need((owner,child) not in windows,'repeated owner child')
                        windows[owner,child]=(children[child],touched.get(event['pair'],-1),tick)
            else:need(kind in ('terminal_Y_translations','terminal'),'unknown schedule event')
        need(all((owner,c) in windows for c in children),'missing owner window')
    family={c:2*plan['pairs'][str(c)]+plan['sides'][str(c)] for c in children}
    shared=defaultdict(list)
    for (o,c),(cap,lo,hi) in windows.items():shared[o,family[c]].append((lo,hi))
    for intervals in shared.values():
        need(all(max(a,c)>=min(b,d) for j,(a,b) in enumerate(intervals) for c,d in intervals[j+1:]),'overlapping source-family windows')
    widths={};oldwidths={};local={};newpins={}
    for tag,pid in [('base','gauge44'),('extra','extra')]:
        final_path=ROOT/f'inputs/WORD-{tag}.json.gz'
        word=read(final_path);oldword=read(OLD/f'old-WORD-{tag}.json.gz');native=read(OLD/f'old-NATIVE-{tag}.json.gz')
        origin=native['helper_origin']
        need(origin==sorted(set(range(224,704))-set(word['removed'])),'helper origin changed')
        need(hashlib.sha256(gzip.decompress(final_path.read_bytes())).hexdigest()==manifest['source_hashes']['new_word_'+tag],'final gzip/source payload bridge')
        widths[pid]=canonical_widths(word,origin);oldwidths[pid]=canonical_widths(oldword,origin)
        need(widths[pid]==manifest['new_widths'][pid],'builder/checker width disagreement')
        need(Counter(widths[pid])=={5:28,6:18,7:53,8:96,9:213},'literal inventory')
        local[pid]={int(r):n for r,n in word['report']['local_profile'].items()}
        need(local[pid]=={1:2482,2:712,3:218,4:109,5:20,8:8},'literal local profile')
        newpins[tag]=sha(final_path)
    rows=fill['rows'];oldrows=oldfill['rows']
    need(fill['columns']==oldfill['columns']==['time','helper','owner_replica','child','port','offset'] and fill['dispatch_overrides']==oldfill['dispatch_overrides']==[],'fill schema')
    need(len(rows)==len(oldrows)==N*408,'role coverage')
    need(all(a[:5]==b[:5] for a,b in zip(rows,oldrows)),'owner assignment moved')
    def walk(data,charges):
        used={};families=set()
        for j,row in enumerate(data):
            need(len(row)==6 and all(type(x)is int for x in row),'noninteger owner role')
            t,h,o,c,p,z=row;need(t==j//408 and h==j%408,'ordered helper identity')
            if h==0:families=set()
            need((o,c) in windows and 0<=p<112,'owner child/port range')
            cap,lo,hi=windows[o,c];need(lo<t<=hi,'owner chronology')
            f=o*896+family[c]*112+p
            need(f not in families,'same-call family distinctness');families.add(f)
            key=(o*K+c)*112+p;prefix=used.get(key,0)
            need(z==prefix,'compact chronological prefix')
            width=charges[provider[tuple(order[t])]][h]
            need(z+width<=cap,'owner capacity')
            used[key]=prefix+width
        idle=Counter()
        for (o,c),(cap,lo,hi) in windows.items():
            for p in range(112):
                left=cap-used.get((o*K+c)*112+p,0);need(left>=0,'negative bin idle')
                if left:idle[left]+=1
        return used,idle
    oldused,oldidle=walk(oldrows,oldwidths);newused,newidle=walk(rows,widths)
    need(oldidle==Counter({1:84240,3:23328,5:4320,48:137952}),'old idle source')
    need(newidle==Counter({3:23328,5:4320,48:137952}),'new idle exact')
    changed=sum(oldused.get(k,0)!=newused.get(k,0) for k in set(oldused)|set(newused))
    need(changed==84240,'changed bin count')
    need(all(newused.get(k,0)-oldused.get(k,0) in (0,1) for k in set(oldused)|set(newused)),'changed bin increment')
    profile=Counter(newidle)
    for call in order:profile.update(local[provider[tuple(call)]])
    need(dict(profile)=={int(r):n for r,n in read(OUT/'PROFILE.json').items()},'global paid histogram')
    need(sum(r*n for r,n in profile.items())==65088*T,'mass deficit')
    controls=[]
    saved=list(rows[1]);rows[1][2:5]=rows[0][2:5]
    try:
        try:walk(rows,widths)
        except ValueError as e:need('same-call family' in str(e),'wrong family rejection');controls.append('same_call_family')
        else:raise ValueError('family mutant accepted')
    finally:rows[1][:]=saved
    saved=rows[0][5];rows[0][5]+=1
    try:
        try:walk(rows,widths)
        except ValueError as e:need('compact chronological' in str(e),'wrong offset rejection');controls.append('offset_hole')
        else:raise ValueError('offset mutant accepted')
    finally:rows[0][5]=saved
    saved=list(rows[0]);bad=next((o,c) for (o,c),(cap,lo,hi) in windows.items() if lo>=0)
    rows[0][2:4]=list(bad)
    try:
        try:walk(rows,widths)
        except ValueError as e:need('owner chronology' in str(e),'wrong chronology rejection');controls.append('closed_owner_window')
        else:raise ValueError('chronology mutant accepted')
    finally:rows[0][:]=saved
    firsttag=provider[tuple(order[0])];saved=widths[firsttag][0];widths[firsttag][0]=windows[rows[0][2],rows[0][3]][0]+1
    try:
        try:walk(rows,widths)
        except ValueError as e:need('owner capacity' in str(e),'wrong overflow rejection');controls.append('overflow')
        else:raise ValueError('overflow mutant accepted')
    finally:widths[firsttag][0]=saved
    result={'status':'PASS_LITERAL_ENLARGED_OWNER_FILL','replicas':T,'calls':N,'roles':len(rows),'owner_windows':len(windows),'physical_bins':len(windows)*112,'changed_bins':changed,'idle_H1_consumed':84240,'idle':dict(sorted(newidle.items())),'profile':dict(sorted(profile.items())),'controls_rejected':controls,'fill_sha256':sha(OUT/'FILL.json.gz'),'final_word_sha256':newpins,'manifest_sha256':sha(OUT/'MANIFEST.json'),'checker_sha256':sha(Path(__file__)),'seconds':time.monotonic()-begin,'scope':'Exact full owner role, chronology, distinctness, capacity, compact-prefix and paid histogram audit. Native/caller/outer admission handled separately.'}

    if '--record' in sys.argv:(OUT/'CHECK.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps(result,sort_keys=True),flush=True)
if __name__=='__main__':main()
