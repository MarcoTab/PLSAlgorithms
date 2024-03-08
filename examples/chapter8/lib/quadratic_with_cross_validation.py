from pls_matrices import pls_matrices
import numpy as np
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
import functools as ft

def quadratic_discrimination(Xes, y, K, d, para, loocv=False):

    p_pls_1 = np.empty_like(y).flatten()
    p_pls_2 = np.empty_like(y).flatten()

    n = Xes.shape[0]

    # Generate random samples
    random_samples = np.random.choice(np.arange(1, n+1), size=n, replace=False)

    bins = None
    folds = None
    if not loocv:
        # Perform binning
        bins = np.linspace(1, n, K+1)
        folds = np.digitize(random_samples, bins, right=False)-1
    else:
        folds = random_samples-1

    for kk in range(K):

        test_idxs = (folds == kk)
        train_idxs = (folds != kk)

        holdout = Xes[train_idxs, :]
        X_test = Xes[test_idxs, :]
        
        y_holdout = y[train_idxs]

        S_X = np.cov(holdout, rowvar=False, ddof=1)
        uY = np.array([int(x) for x in np.unique(y_holdout)])

        S_k = np.array([x for x in map(lambda k: np.cov(holdout[(y_holdout==k).flatten()], rowvar=False, ddof=1), uY)])

        cov_mean = np.array(ft.reduce(lambda a, b: a+b, [x for x in map(lambda k: S_k[k] * np.mean(y_holdout == k), uY)]))

        cov_ss1 = np.array(ft.reduce(lambda a, b: a+b, [x for x in map(lambda k: (S_X - S_k[k]) @ (S_X - S_k[k]).T * np.mean(y_holdout == k), uY)]))

        cov_ss2 = np.array(ft.reduce(lambda a, b: a+b, [x for x in map(lambda k: (cov_mean - S_k[k]) @ (cov_mean - S_k[k]).T * np.mean(y_holdout == k), uY)]))

        dif_d_delta2 = (S_X - cov_mean) @ (S_X - cov_mean).T
        A_al_para = para * dif_d_delta2 + (1-para)*cov_ss2
        A_al1 = cov_ss1
        AUX_para = pls_matrices(A_al_para, cov_mean, d)
        AUX = pls_matrices(A_al1, cov_mean, d)

        proj_data_2 = holdout @ AUX_para['Gamma']
        proj_test_data_2 = X_test @ AUX_para['Gamma']

        proj_data_1 = holdout @ AUX['Gamma']
        proj_test_data_1 = X_test @ AUX['Gamma']

        quad_pls_1 = QuadraticDiscriminantAnalysis()
        quad_pls_1.fit(proj_data_1, y_holdout)

        p_pls_1[test_idxs] = quad_pls_1.predict(proj_test_data_1)

        quad_pls_2 = QuadraticDiscriminantAnalysis()
        quad_pls_2.fit(proj_data_2, y_holdout)
        p_pls_2[test_idxs] = quad_pls_2.predict(proj_test_data_2)

    return {
        "predict_quad_1": p_pls_1,
        "predict_quad_2": p_pls_2
    }


# Trains and predicts on a dataset using QDA. Leave One Out.
def simple_qda(Xes, y):
    n = Xes.shape[0]
    random_samples = np.random.choice(np.arange(1, n+1), size=n, replace=False) - 1

    preds = np.empty_like(y).flatten()

    for i in random_samples:
        test_x = Xes[random_samples == i, :]
        train_x = Xes[random_samples != i, :]

        train_y = y[random_samples != i]


        mod = QuadraticDiscriminantAnalysis()
        mod.fit(train_x, train_y)
        preds[random_samples == i] = mod.predict(test_x)

    return preds