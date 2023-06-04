import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.linear_model import LinearRegression
from sliced import SlicedInverseRegression

#  Xes is the data n times p, y is n by 1, K the number of folds
def lineal_junto_small_R_more_classes(Xes, y, K, d, loocv=False):
    p_pls = np.empty_like(y.reshape(-1,1))
    p_iso = np.empty_like(y.reshape(-1,1))
    p_pfc = np.empty_like(y.reshape(-1,1))
    # d should be smaller than the number of samples

    n = Xes.shape[0]

    # Generate random samples
    random_samples = np.random.choice(np.arange(1, n+1), size=n, replace=False)

    max_encoding = int(np.max(y).item())

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

        auxxx = np.zeros((y[train_idxs].size, max_encoding+1))
        auxxx[np.arange(y[train_idxs].size), [int(x) for x in y[train_idxs]]] = 1

        mod = PLSRegression(n_components=d)
        mod.fit(Xes[train_idxs, :], auxxx)
        AUX = mod.x_loadings_

        proj_data = Xes[train_idxs, :] @ AUX
        proj_test_data = Xes[test_idxs, :] @ AUX

        # PLS
        linear_pls = LinearDiscriminantAnalysis()
        linear_pls.fit(X=proj_data, y=y[train_idxs])

        p_pls[test_idxs] = (linear_pls.predict(proj_test_data)).reshape(-1,1)

        # ISO
        reg = LinearRegression()
        # Inverse regression
        reg.fit(auxxx, Xes[train_idxs, :])
        residuals = Xes[train_idxs, :] - reg.predict(auxxx) 

        _, evecs = np.linalg.eigh(np.cov(residuals, rowvar=False))

        evecs = evecs[:, ::-1]
        evecs = evecs[:, :d]
        
        linear_iso = LinearDiscriminantAnalysis()
        linear_iso.fit(X=Xes[train_idxs, :] @ evecs, y=y[train_idxs])

        p_iso[test_idxs] = (linear_iso.predict(Xes[test_idxs,:] @ evecs)).reshape(-1,1)

        # PFC
        pfc_train_x = None
        pfc_test_x = None
        s0 = None
        if (d==1 or d==2):
            pfc_train_x = proj_data
            pfc_test_x = proj_test_data
        else:
            s0 = SlicedInverseRegression()
            s0.fit(X=proj_data, y=y[train_idxs].flatten())
            vecs = s0.directions_[np.argsort(s0.eigenvalues_)[::-1][:2]].T

            pfc_train_x = proj_data @ vecs
            pfc_test_x = proj_test_data @ vecs
        
        linear_pfc = LinearDiscriminantAnalysis()
        linear_pfc.fit(pfc_train_x, y[train_idxs])

        p_pfc[test_idxs] = (linear_pfc.predict(pfc_test_x)).reshape(-1, 1)
    
    return {
        'predict_linear_pls': p_pls, 
        'predict_linear_iso': p_iso, 
        'predict_linear_pfc': p_pfc
    }

