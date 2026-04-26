import json
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import xgboost as xgb

OUT = "outputs_advanced"

def load_data():
    mat = pd.read_csv('outputs/filtered_matrix.tsv', sep='\t', index_col=0)
    meta = pd.read_csv('outputs/sample_metadata.csv', index_col=0)
    # transpose to samples x features
    X = mat.T
    X.index.name = 'sample'
    # align metadata
    meta = meta.loc[X.index]
    mask = meta['group'] != 'QC'
    X = X.loc[mask]
    y = (meta.loc[mask,'group'] == 'SLE').astype(int).values
    return X, y

def nested_cv(X, y):
    outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    inner = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)

    lr_params = {'C':[0.01,0.1,1,10]}
    xgb_params = {'n_estimators':[50,100], 'max_depth':[3,5]}

    lr_auc = []
    xgb_auc = []
    selected_features = []

    for train_idx, test_idx in outer.split(X, y):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # logistic L1 via liblinear
        lr = LogisticRegression(penalty='l1', solver='liblinear', max_iter=2000)
        gs = GridSearchCV(lr, {'C':lr_params['C']}, cv=inner, scoring='roc_auc')
        gs.fit(X_train, y_train)
        best = gs.best_estimator_
        p = best.predict_proba(X_test)[:,1]
        lr_auc.append(roc_auc_score(y_test, p))
        coef = pd.Series(best.coef_[0], index=X.columns)
        sel = list(coef[coef!=0].index)
        selected_features.extend(sel)

        # xgboost
        clf = xgb.XGBClassifier(use_label_encoder=False, eval_metric='logloss', verbosity=0)
        gs2 = GridSearchCV(clf, xgb_params, cv=inner, scoring='roc_auc')
        gs2.fit(X_train, y_train)
        p2 = gs2.best_estimator_.predict_proba(X_test)[:,1]
        xgb_auc.append(roc_auc_score(y_test, p2))

    return lr_auc, xgb_auc, selected_features

def summarize(lr_auc, xgb_auc, sel):
    out1 = {
        'lr_cv_auc_mean': float(np.mean(lr_auc)),
        'lr_cv_auc_std': float(np.std(lr_auc)),
        'xgb_cv_auc_mean': float(np.mean(xgb_auc)),
        'xgb_cv_auc_std': float(np.std(xgb_auc)),
    }
    # stability: count features selected across folds
    from collections import Counter
    ctr = Counter(sel)
    n_outer = 5
    stable = [f for f,c in ctr.items() if c >= 3]
    pd.Series(ctr).sort_values(ascending=False).to_csv(f"{OUT}/nested_feature_counts.tsv", sep='\t')
    pd.Series(stable).to_csv(f"{OUT}/panel_selected_features.tsv", index=False, header=False)
    with open(f"{OUT}/nested_cv_results.json","w") as f:
        json.dump(out1, f, indent=2)
    return out1, stable

def main():
    X, y = load_data()
    lr_auc, xgb_auc, sel = nested_cv(X, y)
    summary, stable = summarize(lr_auc, xgb_auc, sel)
    print('Nested CV done. Summary:', summary)

if __name__ == '__main__':
    main()
