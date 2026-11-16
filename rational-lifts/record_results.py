"""Regenerate lift tables from the saved Magma output and LMFDB cache."""
import ast
import csv
import json
import re
from pathlib import Path

root = Path(__file__).resolve().parent
out = root
out.mkdir(parents=True, exist_ok=True)
geometric = json.loads((root.parent/'congruence-enumeration/results.json').read_text())
db = {r['label']:r for r in json.loads((root/'lmfdb-candidates.json').read_text())}
gens = {}
for m in re.finditer(r'Append\(~RationalLifts,<(\d+),(\d+),(\d+),(\d+),(\[.*?\])>\);', (root/'lift_representatives.m').read_text(), re.S):
    i,N,h,k = map(int,m.groups()[:4])
    gens[i,N,k] = ast.literal_eval(m[5])
matches = []
for line in (root/'matches.txt').read_text().splitlines()[1:]:
    i,N,h,k,L,label,coarse,conj = line.split('|')
    i,N,h,k,L = map(int,(i,N,h,k,L))
    matches.append(dict(id=i,modulus=N,H_level=h,lift_number=k,G_level=L,label=label,coarse_label=coarse,
                        computed_generators=[[a%L for a in v] for v in gens[i,N,k]],
                        conjugator=ast.literal_eval(conj),lmfdb_generators=db[label]['generators']))
unique = {}
for r in matches:
    if r['label'] in unique:
        assert unique[r['label']]['id'] == r['id']
    else:
        unique[r['label']] = r
def key(label):
    return [int(p) if p.isdigit() else p for p in re.split(r'[.\-]',label)]
for h in geometric:
    found = [r for r in unique.values() if r['id']==h['id']]
    h['exact_level_labels'] = sorted([r['label'] for r in found if r['G_level']==h['level']],key=key)
    h['additional_proof_modulus_labels'] = sorted([r['label'] for r in found if r['G_level']!=h['level']],key=key)
    h['admits_rational_lift'] = bool(found)
assert len(unique)==636
assert sum(len(h['exact_level_labels']) for h in geometric)==457
assert sum(h['admits_rational_lift'] for h in geometric)==61
assert [h['id'] for h in geometric if not h['admits_rational_lift']]==[2,6,7]
(out/'geometric-groups.json').write_text(json.dumps(geometric,indent=2)+'\n')
(out/'all-lifts.json').write_text(json.dumps(sorted(unique.values(),key=lambda r:(r['id'],key(r['label']))),indent=2)+'\n')
with (out/'labels-by-geometric-group.csv').open('w') as f:
    w=csv.writer(f); w.writerow(['id','CP','level_N','projective_level','genus','SL2_index','admits_Q_lift','exact_level_lift_count','example_fine_label','all_exact_level_fine_labels','additional_labels_mod_2Nprime'])
    for h in geometric:
        labels=h['exact_level_labels']
        w.writerow([h['id'],h['cp'],h['level'],h['projective_level'],h['genus'],2*h['mu'],h['admits_rational_lift'],len(labels),labels[0] if labels else '', '; '.join(labels),'; '.join(h['additional_proof_modulus_labels'])])
with (out/'all-lifts.csv').open('w') as f:
    w=csv.writer(f); w.writerow(['id','CP','H_level','G_level','fine_label','coarse_label','LMFDB_url'])
    for r in sorted(unique.values(),key=lambda r:(r['id'],key(r['label']))):
        w.writerow([r['id'],geometric[r['id']-1]['cp'],r['H_level'],r['G_level'],r['label'],r['coarse_label'],'https://beta.lmfdb.org/ModularCurve/Q/'+r['label']])

print("PASS: 61 geometric classes, 457 exact-level lifts, 636 proof-modulus labels")
