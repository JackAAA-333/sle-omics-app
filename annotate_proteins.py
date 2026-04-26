import requests
import time
import pandas as pd
from pathlib import Path

OUT = Path('outputs_prot/protein_annotations.tsv')

def fetch_uniprot(acc):
    url = f'https://rest.uniprot.org/uniprotkb/{acc}.json'
    try:
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            return None
        j = r.json()
        pname = ''
        gname = ''
        # recommended name
        pdsc = j.get('proteinDescription', {})
        rec = pdsc.get('recommendedName') if pdsc else None
        if rec and 'fullName' in rec:
            pname = rec['fullName'].get('value','')
        # gene
        genes = j.get('genes', [])
        if genes:
            g = genes[0].get('geneName') or genes[0].get('primaryGeneName') or {}
            if isinstance(g, dict):
                gname = g.get('value','')
            else:
                gname = str(g)
        return {'protein_name':pname, 'gene_name':gname}
    except Exception:
        return None

def translate_zh(text):
    if not text:
        return ''
    url = 'https://libretranslate.de/translate'
    payload = {'q': text, 'source': 'en', 'target': 'zh', 'format': 'text'}
    try:
        r = requests.post(url, data=payload, timeout=10)
        if r.status_code == 200:
            j = r.json()
            return j.get('translatedText','')
    except Exception:
        return ''
    return ''

def gather_ids():
    ids = set()
    p1 = Path('outputs_prot/xgb_shap_prot.tsv')
    p2 = Path('outputs_prot/differential_prot.tsv')
    for p in (p1,p2):
        if not p.exists():
            continue
        try:
            df = pd.read_csv(p, sep='\t', header=0, index_col=None)
        except Exception:
            continue
        # look for first column that looks like accession (starts with P or Q and 5-6 chars)
        if df.shape[1] >= 1:
            col0 = df.columns[0]
            for v in df[col0].astype(str).values:
                ids.add(v.strip())
        # also check index
        try:
            idx = pd.read_csv(p, sep='\t', index_col=0)
            for v in idx.index.astype(str):
                ids.add(v.strip())
        except Exception:
            pass
    return sorted([i for i in ids if i])

def main():
    ids = gather_ids()
    rows = []
    for i,acc in enumerate(ids):
        info = fetch_uniprot(acc)
        time.sleep(0.4)
        if info is None:
            pname = ''
            gname = ''
        else:
            pname = info.get('protein_name','')
            gname = info.get('gene_name','')
        # try translate
        zh = ''
        try:
            zh = translate_zh(pname)
            time.sleep(0.3)
        except Exception:
            zh = ''
        rows.append({'protein_id':acc, 'protein_name':pname, 'gene_name':gname, 'protein_name_zh':zh})
        if (i+1) % 20 == 0:
            print(f'Processed {i+1}/{len(ids)}')

    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, sep='\t', index=False)
    print('Protein annotations written to', OUT)

if __name__=='__main__':
    main()
