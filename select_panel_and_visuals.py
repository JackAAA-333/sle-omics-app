import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve, auc

FREQ = 'outputs_advanced/nested_cv_strict_feature_freq.tsv'
METAB = 'outputs/filtered_matrix.tsv'
PROT = 'outputs_prot/filtered_prot_matrix.tsv'
META = 'outputs/sample_metadata.csv'
OUT_PANEL = 'outputs_advanced/final_panel_by_freq.tsv'
OUT_COEF = 'outputs_advanced/panel_coef_heatmap.png'
OUT_ROCBOX = 'outputs_advanced/panel_cv_auc_boxplot.png'
OUT_ROC = 'outputs_advanced/panel_cv_roc_mean.png'

def load_matrix(path):
    return pd.read_csv(path, sep='\t', index_col=0)

def main(top_n=10):
    freq = pd.read_csv(FREQ, sep='\t')
    selected = freq.head(top_n)['feature'].tolist()
    metab = load_matrix(METAB)
    prot = load_matrix(PROT)
    meta = pd.read_csv(META, index_col=0)
    # assemble X samples x features
    data = {}
    for f in selected:
        if f in metab.index:
            data[f] = metab.loc[f]
        elif f in prot.index:
            data[f] = prot.loc[f]
    X = pd.DataFrame(data).fillna(0)
    # align samples and remove QC
    if 'group' in meta.columns:
        y_full = meta['group'].map({'SLE':1,'HC':0})
    else:
        y_full = meta['label'].map({'SLE':1,'HC':0})
    mask = y_full.notna()
    X = X.loc[mask.index[mask]]
    y = y_full.loc[mask.index[mask]].values

    # save panel
    pd.DataFrame({'feature':selected}).to_csv(OUT_PANEL, sep='\t', index=False)

    # cross-validated logistic to get per-fold AUCs and coefs
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    aucs = []
    rocs = []
    coefs = []
    for train_idx, test_idx in skf.split(X, y):
        Xtr, Xte = X.iloc[train_idx], X.iloc[test_idx]
        ytr, yte = y[train_idx], y[test_idx]
        clf = LogisticRegression(solver='liblinear', max_iter=2000)
        clf.fit(Xtr, ytr)
        probs = clf.predict_proba(Xte)[:,1]
        aucs.append(roc_auc_score(yte, probs))
        fpr, tpr, _ = roc_curve(yte, probs)
        rocs.append((fpr, tpr))
        coefs.append(pd.Series(clf.coef_[0], index=X.columns))

    # save coefs heatmap
    coef_df = pd.DataFrame(coefs)
    coef_df.index = [f'fold{i+1}' for i in range(len(coef_df))]
    plt.figure(figsize=(6, max(3, len(selected)*0.4)))
    sns.heatmap(coef_df.T, cmap='vlag', center=0)
    plt.title('Panel coefficients across CV folds')
    plt.tight_layout()
    plt.savefig(OUT_COEF, dpi=150)
    plt.close()

    # ROC boxplot
    plt.figure(figsize=(4,3))
    sns.boxplot(data=aucs)
    plt.ylabel('AUC')
    plt.title('Cross-validated AUCs')
    plt.tight_layout()
    plt.savefig(OUT_ROCBOX, dpi=150)
    plt.close()

    # mean ROC curve
    # interpolate to mean
    mean_fpr = np.linspace(0,1,100)
    tprs = []
    for fpr,tpr in rocs:
        tprs.append(np.interp(mean_fpr, fpr, tpr))
    mean_tpr = np.mean(tprs, axis=0)
    mean_auc = np.mean(aucs)
    plt.figure(figsize=(4.5,4.5))
    plt.plot(mean_fpr, mean_tpr, label=f'Mean AUC={mean_auc:.3f}')
    plt.plot([0,1],[0,1],'--', color='gray')
    plt.xlabel('FPR'); plt.ylabel('TPR'); plt.title('Mean ROC')
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT_ROC, dpi=150)
    plt.close()

    print('Panel selection and visuals saved:', OUT_PANEL)

if __name__=='__main__':
    main()
