import pandas as pd
from pathlib import Path

IN = Path('outputs_prot/protein_annotations_top50.tsv')
OUT_CSV = Path('outputs_prot/protein_annotations_top50_cn.csv')
OUT_XLSX = Path('outputs_prot/protein_annotations_top50_supplement.xlsx')

def main():
    if not IN.exists():
        raise SystemExit(f'Input not found: {IN}')
    df = pd.read_csv(IN, sep='\t', dtype=str)
    df = df.fillna('')
    df_cn = df.rename(columns={
        'protein_id': '蛋白ID',
        'protein_name': '蛋白英文名',
        'gene_name': '基因名',
        'protein_name_zh': '蛋白中文名'
    })
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df_cn.to_csv(OUT_CSV, index=False)
    try:
        df_cn.to_excel(OUT_XLSX, index=False)
    except Exception as e:
        print('写入 Excel 失败:', e)
    print('Wrote:', OUT_CSV, OUT_XLSX)

if __name__ == '__main__':
    main()
