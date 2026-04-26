import json
import requests
import time
import xml.etree.ElementTree as ET

IN = 'outputs_advanced/pubmed_summary.json'
OUT = 'outputs_advanced/pubmed_expanded.tsv'

def fetch_pmid(pmid):
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    params = {'db':'pubmed','id':pmid,'retmode':'xml'}
    r = requests.get(url, params=params, timeout=20)
    r.raise_for_status()
    return r.text

def parse_article(xml_text):
    root = ET.fromstring(xml_text)
    art = root.find('.//Article')
    title = art.findtext('ArticleTitle') if art is not None else ''
    abstract_nodes = root.findall('.//AbstractText')
    abstract = ' '.join([ET.tostring(n, encoding='unicode', method='text') for n in abstract_nodes])
    return title, abstract

def main():
    data = json.load(open(IN))
    rows = []
    for feat, info in data.items():
        pmids = info.get('pmids', []) or []
        if not pmids:
            rows.append([feat, info.get('name', ''), '', ''])
            continue
        for pmid in pmids:
            try:
                xml = fetch_pmid(pmid)
                title, abstract = parse_article(xml)
            except Exception as e:
                title = ''
                abstract = ''
            rows.append([feat, info.get('name',''), pmid, title + ' ' + abstract])
            time.sleep(0.34)
    import csv
    with open(OUT,'w',newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['feature','name','pmid','title_abstract'])
        w.writerows(rows)
    print('PubMed expansion complete:', OUT)

if __name__=='__main__':
    main()
