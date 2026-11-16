"""Export signed subgroup actions and match the 64 congruence classes."""
import ast
import csv
import json
import re
from pathlib import Path

root=Path(__file__).resolve().parent
out=root
out.mkdir(parents=True,exist_ok=True)
def canonical(s,r):
    n=len(s);best=None
    for base in range(n):
        order=[base];lab={base:0};code=[]
        for a in order:
            for perm in (s,r):
                b=perm[a]
                if b not in lab:lab[b]=len(order);order.append(b)
                code.append(lab[b])
        assert len(order)==n
        key=tuple(code)
        if best is None or key<best:best=key
    return best

maps={}
for line in (root/'maps.txt').read_text().splitlines()[1:]:
    f=line.split('|');i,mu,e3,g,c,reg,aut,raw,cl=map(int,f[:9])
    r=list(map(int,f[9].split(',')));s=list(map(int,f[10].split(',')))
    vertices=[];seen=set()
    for j in range(mu):
        if j in seen:continue
        v=[];k=j
        while k not in seen:seen.add(k);v.append(k);k=r[k]
        vertices.append(v)
    tree={a for v in vertices[1:] for a in (min(v),s[min(v)])}
    free=[j for j in range(mu) if j<s[j] and j not in tree]
    maps[i]=dict(id=i,mu=mu,e3=e3,genus=g,cusps=c,regular_cusps=reg,automorphisms=aut,
                 raw_lifts=raw,lift_classes=cl,R=r,S=s,free_edges=free,lift_bits=list(map(int,f[11].split(','))))
def lift(m,bits):
    s=m['S'];r=m['R'];n=m['mu'];voltage=[int(i>s[i]) for i in range(n)]
    for j,i in enumerate(m['free_edges']):
        if bits>>j&1:voltage[i]^=1;voltage[s[i]]^=1
    sl=[2*s[i]+(e^voltage[i]) for i in range(n) for e in range(2)]
    rl=[2*r[i]+e for i in range(n) for e in range(2)]
    return sl,rl
old={}
for line in (root/'congruence_actions.txt').read_text().splitlines()[1:]:
    i,rs,ss=line.split('|');key=canonical(list(map(int,ss.split(','))),list(map(int,rs.split(','))))
    assert key not in old;old[key]=int(i)
records=[];found=set()
for line in (root/'verified-lifts.txt').read_text().splitlines()[1:]:
    i,b,cong,n,g,reg,irr=line.split('|');i,b,n,g=map(int,(i,b,n,g));yes=cong=='true'
    rec=dict(id=f'{i}.{b}',map_id=i,lift_bits=b,SL2_index=n,genus=g,congruence=yes,
             regular_widths=ast.literal_eval(reg),irregular_widths=ast.literal_eval(irr))
    if yes:
        sl,rl=lift(maps[i],b);key=canonical(sl,rl);assert key in old
        rec['previous_congruence_id']=old[key];found.add(old[key])
    records.append(rec)
assert found==set(range(1,65))
assert len(records)==47244 and sum(not r['congruence'] for r in records)==47180
with (out/'subgroups.csv').open('w') as f:
    w=csv.writer(f);w.writerow(['id','SL2_index','genus','e3','congruence','regular_widths','irregular_widths','previous_congruence_id'])
    for r in records:w.writerow([r['id'],r['SL2_index'],r['genus'],maps[r['map_id']]['e3'],r['congruence'],r['regular_widths'],r['irregular_widths'],r.get('previous_congruence_id','')])
(out/'maps.json').write_text(json.dumps(list(maps.values()),separators=(',',':'))+'\n')
(out/'subgroups.json').write_text(json.dumps(records,separators=(',',':'))+'\n')
cal_arrays=re.findall(r'\[[^\]]*\]',(root/'calegari-action.txt').read_text())
cr,cs=[ast.literal_eval(x) for x in cal_arrays]
ck=canonical(cs,cr)
cal=[m for m in maps.values() if m['mu']==36 and m['e3']==0 and m['genus']==0 and canonical(m['S'],m['R'])==ck]
assert len(cal)==1
cal=cal[0];assert len(cal['lift_bits'])==1
print('Calegari geometric class',cal['id'],cal['lift_bits'])
(out/'calegari-identification.json').write_text(json.dumps(dict(map_id=cal['id'],lift_bits=cal['lift_bits'][0],
    subgroup_id=f"{cal['id']}.{cal['lift_bits'][0]}",method='Numerical monodromy repeated with 240 and 480 steps; exact permutation comparison',
    provenance='Frank Calegari, An explicit rank-four Asai-type motive of Hodge type (2,0)^2+(0,2)^2, unpublished supplied note, 2026',
    source_note='Unpublished supplied note; not redistributed'),indent=2)+'\n')
print("PASS: 47244 total classes, 47180 noncongruence; all 64 congruence classes matched")
