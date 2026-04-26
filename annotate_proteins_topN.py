import requests
import time
import pandas as pd
from pathlib import Path

OUT = Path('outputs_prot/protein_annotations_top50.tsv')

def fetch_uniprot(acc):
    url = f'https://rest.uniprot.org/uniprotkb/{acc}.json'
    try:
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            return None
        j = r.json()
        pname = ''
        gname = ''
        pdsc = j.get('proteinDescription', {})
        rec = pdsc.get('recommendedName') if pdsc else None
        if rec and 'fullName' in rec:
            pname = rec['fullName'].get('value','')
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

def main(top_n=50):
    p = Path('outputs_prot/xgb_shap_prot.tsv')
    if not p.exists():
        raise SystemExit('xgb_shap_prot.tsv not found')
    df = pd.read_csv(p, sep='\t', header=0, index_col=0)
    ids = list(df.index.astype(str))[:top_n]
    rows = []
    for i, acc in enumerate(ids):
        info = fetch_uniprot(acc)
        time.sleep(0.3)
        if info is None:
            pname = ''
            gname = ''
        else:
            pname = info.get('protein_name','')
            gname = info.get('gene_name','')
        zh = ''
        try:
            zh = translate_zh(pname)
        except Exception:
            zh = ''
        rows.append({'protein_id':acc, 'protein_name':pname, 'gene_name':gname, 'protein_name_zh':zh})
        print(f'Annotated {i+1}/{len(ids)}')
    OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT, sep='\t', index=False)
    print('Wrote', OUT)

if __name__=='__main__':
    main()
