import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt

PANEL = 'outputs_advanced/panel_selected_features.tsv'
METAB = 'outputs/filtered_matrix.tsv'
PROT = 'outputs_prot/filtered_prot_matrix.tsv'
META = 'outputs/sample_metadata.csv'
OUT_CSV = 'outputs_advanced/final_panel.csv'
OUT_ROC = 'outputs_advanced/final_panel_roc.png'

def load_matrix(path):
    df = pd.read_csv(path, sep='\t', index_col=0)
    return df

def main():
    panel = pd.read_csv(PANEL, sep='\t', header=None)
    features = panel[0].tolist()
    metab = load_matrix(METAB)
    prot = load_matrix(PROT)
    # collect feature series into a dataframe with samples as rows
    used = []
    data = {}
    for f in features:
        if f in metab.index:
            data[f] = metab.loc[f]
            used.append(f)
        elif f in prot.index:
            data[f] = prot.loc[f]
            used.append(f)
    if not data:
        raise SystemExit('No panel features found in matrices')
    df_feat = pd.DataFrame(data)
    # samples x features
    X = df_feat.fillna(0)
    meta = pd.read_csv(META, index_col=0)
    # build y aligned to X rows (samples); use 'group' column
    if 'group' in meta.columns:
        y = meta.loc[X.index,'group'].map({'SLE':1,'HC':0})
    elif 'label' in meta.columns:
        y = meta.loc[X.index,'label'].map({'SLE':1,'HC':0})
    else:
        raise SystemExit('sample metadata missing group/label column')

    # remove QC or unknown groups
    mask = y.notna()
    X2 = X.loc[mask.index[mask]].fillna(0)
    y2 = y.loc[mask.index[mask]]
    clf = LogisticRegression(penalty='l2', solver='liblinear', max_iter=2000)
    clf.fit(X2, y2)
    probs = clf.predict_proba(X.fillna(0))[:,1]
    probs2 = clf.predict_proba(X2)[:,1]
    fpr, tpr, _ = roc_curve(y2, probs2)
    roc_auc = auc(fpr, tpr)

    # save panel CSV with coefficients and probs
    coef = pd.Series(clf.coef_[0], index=X.columns)
    out = pd.DataFrame({'feature':X.columns, 'coef':coef.values})
    out.to_csv(OUT_CSV, index=False)

    # save probs per sample
    pred = pd.DataFrame({'sample':X2.index, 'prob':probs2, 'label':y2.values})
    pred.to_csv('outputs_advanced/final_panel_predictions.tsv', sep='\t', index=False)

    plt.figure(figsize=(5,5))
    plt.plot(fpr, tpr, label=f'AUC={roc_auc:.3f}')
    plt.plot([0,1],[0,1],'--',color='gray')
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.title('Final panel ROC')
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT_ROC, dpi=150)
    print('Final panel exported:', OUT_CSV, OUT_ROC)

if __name__=='__main__':
    main()
