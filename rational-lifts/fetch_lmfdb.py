"""Cache the public LMFDB fine-group records needed for the finite lift test."""
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

root = Path(__file__).resolve().parent
(root / "sources").mkdir(exist_ok=True)
levels = sorted({int(line.split('|')[2]) for line in (root/'lifts.txt').read_text().splitlines()[1:]})
fields = 'label,level,index,genus,generators,contains_negative_one,coarse_label,CPlabel,sl2level,sl2label,determinant_label'
query = dict(_format='json', _fields=fields, level='py'+repr({'$in':levels}),
             genus="py{'$lte': 1}", index="py{'$lte': 72}", contains_negative_one='pyFalse')
records = []
offset = 0
while True:
    target = root/'sources'/f'lmfdb-page-{offset:05}.json'
    if target.exists():
        page = json.loads(target.read_text())
    else:
        url = 'https://beta.lmfdb.org/api/gps_gl2zhat_fine/?'+urlencode(query | {'_offset':offset})
        req = Request(url, headers={'Cookie':'human=1'})
        with urlopen(req, timeout=45) as response:
            data = response.read()
        page = json.loads(data)
        target.write_bytes(data)
    records.extend(page['data'])
    print(offset, len(page['data']), flush=True)
    if len(page['data']) < 100:
        break
    offset = page['offset']
assert len({r['label'] for r in records}) == len(records)
(root/'lmfdb-candidates.json').write_text(json.dumps(records, indent=2)+'\n')
with (root/'lmfdb_candidates.m').open('w') as f:
    f.write('// Downloaded from the LMFDB; see fetch_lmfdb.py and sources/*.json.\nLMFDB := [* *];\n')
    for r in records:
        f.write('Append(~LMFDB,<'+','.join([json.dumps(r['label']),str(r['level']),str(r['index']),str(r['genus']),json.dumps(r['coarse_label']),json.dumps(r['generators'])])+'>);\n')
print('Saved', len(records), 'records; levels', levels)
