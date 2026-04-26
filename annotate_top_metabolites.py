import pandas as pd
import requests
import time

IN_TOP = 'outputs/top100_candidates.tsv'
OUT = 'outputs_advanced/metabolite_annotations_top100.tsv'

def kegg_find_by_name(name):
    if not isinstance(name,str) or not name.strip():
        return ''
    q = requests.utils.requote_uri(name)
    url = f'https://rest.kegg.jp/find/compound/{q}'
    r = requests.get(url, timeout=10)
    if r.status_code!=200 or not r.text.strip():
        return ''
    line = r.text.splitlines()[0]
    pid = line.split('\t')[0]
    return pid.split(':')[-1]

def main():
    df = pd.read_csv(IN_TOP, sep='\t', index_col=0)
    rows = []
    for feat in df.index[:100]:
        name = df.loc[feat].get('Metabolite Name','') if 'Metabolite Name' in df.columns else ''
        k = ''
        if pd.notna(name) and name:
            try:
                k = kegg_find_by_name(name)
            except Exception:
                k = ''
            time.sleep(0.3)
        rows.append({'feature':feat,'name':name,'kegg':k})
    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep='\t', index=False)
    print('Top100 metabolite annotation complete:', OUT)

if __name__=='__main__':
    main()
