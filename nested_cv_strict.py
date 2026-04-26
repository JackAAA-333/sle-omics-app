import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import roc_auc_score
import json

METAB = 'outputs/filtered_matrix.tsv'
PROT = 'outputs_prot/filtered_prot_matrix.tsv'
META = 'outputs/sample_metadata.csv'
OUT_SUM = 'outputs_advanced/nested_cv_strict.json'
OUT_FREQ = 'outputs_advanced/nested_cv_strict_feature_freq.tsv'

def load_matrix(path):
    return pd.read_csv(path, sep='\t', index_col=0)

def main():
    metab = load_matrix(METAB)
    prot = load_matrix(PROT)
    meta = pd.read_csv(META, index_col=0)
    # combine by intersecting samples
    samples = [s for s in metab.columns if s in meta.index and s in prot.columns]
    # if no overlap with prot, use metab only
    if not samples:
        samples = [s for s in metab.columns if s in meta.index]
        X = metab[samples].T
    else:
        X = pd.concat([metab[samples], prot[samples]], axis=0).T
    # map group column; filter out QC
    if 'group' in meta.columns:
        y = meta.loc[X.index,'group'].map({'SLE':1,'HC':0})
    elif 'label' in meta.columns:
        y = meta.loc[X.index,'label'].map({'SLE':1,'HC':0})
    else:
        raise SystemExit('sample metadata missing group/label column')
    mask = y.notna()
    X = X.loc[mask.index[mask]]
    y = y.loc[mask.index[mask]].values

    outer_reps = 5
    outer_splits = 5
    aucs = []
    feat_counts = {}
    for rep in range(outer_reps):
        skf = StratifiedKFold(n_splits=outer_splits, shuffle=True, random_state=42+rep)
        for train_idx, test_idx in skf.split(X, y):
            X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]
            # inner CV for hyperparam with L1 selection
            inner = StratifiedKFold(n_splits=5, shuffle=True, random_state=1)
            lr = LogisticRegressionCV(Cs=10, cv=inner, penalty='l1', solver='liblinear', scoring='roc_auc', max_iter=2000)
            lr.fit(X_train.fillna(0), y_train)
            probs = lr.predict_proba(X_test.fillna(0))[:,1]
            aucs.append(roc_auc_score(y_test, probs))
            coefs = lr.coef_[0]
            for f,c in zip(X.columns, coefs):
                if abs(c)>1e-6:
                    feat_counts[f] = feat_counts.get(f,0)+1

    summary = {'n_runs': len(aucs), 'auc_mean': float(np.mean(aucs)), 'auc_std': float(np.std(aucs))}
    with open(OUT_SUM,'w') as f:
        json.dump(summary, f, indent=2)
    # save frequencies
    freq = pd.DataFrame(list(feat_counts.items()), columns=['feature','count']).sort_values('count', ascending=False)
    freq.to_csv(OUT_FREQ, sep='\t', index=False)
    print('Nested strict CV complete:', OUT_SUM, OUT_FREQ)

if __name__=='__main__':
    main()
