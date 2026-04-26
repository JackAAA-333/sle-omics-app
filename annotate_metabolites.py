import pandas as pd
import requests
import time

IN_XLS = 'metab.xlsx'
OUT = 'outputs_advanced/metabolite_annotations.tsv'

def kegg_conv_from_hmdb(hmdb_id):
    url = f'https://rest.kegg.jp/conv/compound/hmdb:{hmdb_id}'
    r = requests.get(url, timeout=10)
    if r.status_code!=200 or not r.text.strip():
        return None
    # response like: hmdb:HMDB00001\tcpd:C00001
    line = r.text.splitlines()[0]
    parts = line.split('\t')
    if len(parts)>=2:
        return parts[1].split(':')[-1]
    return None

def kegg_find_by_name(name):
    if not isinstance(name,str) or not name.strip():
        return None
    q = requests.utils.requote_uri(name)
    url = f'https://rest.kegg.jp/find/compound/{q}'
    r = requests.get(url, timeout=10)
    if r.status_code!=200 or not r.text.strip():
        return None
    line = r.text.splitlines()[0]
    pid = line.split('\t')[0]
    return pid.split(':')[-1]

def main():
    df = pd.read_excel(IN_XLS, sheet_name=0)
    # attempt to find standard columns
    cols = df.columns.tolist()
    # common names
    met_col = next((c for c in cols if 'Met ID' in c or c.lower().startswith('met id') or c.lower().startswith('met_id')), cols[0])
    name_col = next((c for c in cols if 'Name' in c or 'Metabolite' in c), None)
    hmdb_col = next((c for c in cols if 'HMDB' in c.upper()), None)
    kegg_col = next((c for c in cols if 'KEGG' in c.upper()), None)

    out_rows = []
    for _, row in df.iterrows():
        met = row.get(met_col)
        name = row.get(name_col) if name_col else ''
        hmdb = row.get(hmdb_col) if hmdb_col else None
        kegg = row.get(kegg_col) if kegg_col else None
        kegg = str(kegg).strip() if pd.notna(kegg) else ''
        if not kegg and pd.notna(hmdb):
            hmdb_id = str(hmdb).strip()
            try:
                k = kegg_conv_from_hmdb(hmdb_id)
                if k:
                    kegg = k
            except Exception:
                kegg = ''
            time.sleep(0.2)
        if not kegg and pd.notna(name):
            try:
                k = kegg_find_by_name(str(name))
                if k:
                    kegg = k
            except Exception:
                kegg = ''
            time.sleep(0.2)
        out_rows.append({'Met ID':met, 'Name':name, 'HMDB':hmdb if pd.notna(hmdb) else '', 'KEGG':kegg})

    out = pd.DataFrame(out_rows)
    out.to_csv(OUT, sep='\t', index=False)
    print('Annotation complete:', OUT)

if __name__=='__main__':
    main()
