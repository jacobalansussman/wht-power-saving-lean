import re,sys
def lists(path):
    t=open(path,encoding='utf-8').read()
    out={}
    for m in re.finditer(r'def\s+(\w+)\s*:\s*List[^:=]*:=\s*(\[[^\]]*\](?:\s*\]|))',t):
        pass
    # generic: name := [ ... ] possibly nested lists
    for m in re.finditer(r'def\s+(\w+)\s*:[^=]*?:=\s*\[',t):
        name=m.group(1); i=m.end()-1; depth=0; j=i
        while True:
            c=t[j]
            if c=='[': depth+=1
            elif c==']':
                depth-=1
                if depth==0: break
            j+=1
        out[name]=t[i:j+1]
    return out
neg=lists(sys.argv[1]); pos=lists(sys.argv[2]); sch=lists(sys.argv[3])
print('neg names',list(neg)); print('pos names',list(pos)); print('sched names',[k for k in sch][:20])
allpos=dict(pos); allpos.update({k:v for k,v in sch.items()})
import ast
for k,v in neg.items():
    key = k if k in allpos else ('schedRes' if k=='res' else None)
    if key is None or key not in allpos: print('no counterpart for',k); continue
    a=ast.literal_eval(v); b=ast.literal_eval(allpos[key])
    if len(a)!=len(b): print(k,'length differs',len(a),len(b)); continue
    def flat(x):
        if isinstance(x,list):
            for y in x: yield from flat(y)
        else: yield x
    fa=list(flat(a)); fb=list(flat(b))
    d=[(i,fb[i],fa[i]) for i in range(min(len(fa),len(fb))) if fa[i]!=fb[i]]
    print(k,'vs',key,'len',len(a),'flat',len(fa),len(fb),'diffs',len(d),d[:5])
