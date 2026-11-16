"""Regenerate the machine-readable congruence census from results.txt."""
from pathlib import Path
from collections import Counter
import re,ast,json,csv
p=Path(__file__).resolve().parent
rows=[]
for block in re.split(r'\n(?=\d+[A-Z]+[01]\|)',(p/'results.txt').read_text())[1:]:
 f=block.split('|');n=int(re.match(r'\d+',f[0])[0])
 rows.append(dict(id=len(rows)+1,cp=f[0],reflection=int(f[1]),projective_level=n,modulus=2*n,level=int(f[2]),mu=int(f[3]),genus=int(f[4]),e3=int(f[5]),regular_widths=ast.literal_eval(f[6]),irregular_widths=ast.literal_eval(f[7]),orbit_size=int(f[8]),generators=ast.literal_eval(f[9])))
assert len(rows)==64 and sum(r['orbit_size'] for r in rows)==558
assert Counter(r['genus'] for r in rows)=={0:32,1:32}
(p/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
with (p/'results.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print("PASS: 64 congruence classes, 558 actual subgroups; JSON and CSV regenerated")
