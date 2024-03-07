import warnings
warnings.filterwarnings("ignore")

import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.cross_decomposition import PLSRegression
import pandas as pd
from tqdm import tqdm
from tabulate import tabulate

np.random.seed(5)

# PPLS Partial OLS and ENV Partial table 6.2

# Columns 2,3,4

df = pd.read_csv("data/chapter6/Process_data_Skagerberg1992 - Sheet1.csv")

# For reduction
X1 = df.iloc[:, :20].to_numpy()

# No reduction
X2 = df.iloc[:, 20:22].to_numpy()

Y = np.log(df.iloc[:, 22:28].to_numpy())

Y = Y @ (np.diag(np.diag(np.cov(Y, ddof=1,rowvar=False)) ** (-1/2)))

# df = df.assign(Y=Y)

aux = df.to_numpy()

Dmax = 18

n = len(df)

pls_predict = np.zeros((n, Dmax, 6))
ols_predict = np.zeros((n, Dmax, 6))
# TODO no envelope module in python yet
# env_predict = np.zeros((n, Dmax, 6))

# PPLS partial OLS and ENV partial

for d in tqdm(range(Dmax), desc="Training pls, ols, and env models"):
    for i in range(n):
        # Training
        aux1 = np.delete(aux, i, 0)
        aux1_Y = np.delete(Y, i, 0)

        # Testing
        aux2 = aux[i]
        aux2_Y = Y[i]

        # Matrix X to reduce
        AA = aux1[:, :20]

        # print((aux1).shape, (aux1[df.columns.get_loc("S"), :]).shape)
        
        # Regression of the continuous X that we do not reduce
        aux1_Tw_S = np.stack((aux1[:, df.columns.get_loc("Tw")], aux1[:, df.columns.get_loc("S")])).T
        A = LinearRegression()
        A.fit(X=aux1_Tw_S, y=AA)

        # Save residuals
        aux1_SS = AA - A.predict(aux1_Tw_S)

        # Predict on new data
        aux2_Tw_S = np.stack((aux2[df.columns.get_loc("Tw")], aux2[df.columns.get_loc("S")])).reshape(-1,1).T
        A_1 = A.predict(aux2_Tw_S)

        # Center response in training
        aux1_ynew = aux1_Y - aux1_Y.mean(axis=0)

        # Do PLS on the training residuals
        m = PLSRegression(n_components=d+1, scale=False)
        m.fit(X=aux1_SS, Y=aux1_ynew)

        # Perform regression on non-reduced variables
        j = LinearRegression()
        j.fit(X=aux1_Tw_S, y=aux1_ynew)

        # Predict on the new datapoint from non-reduced model
        j_2 = j.predict(aux2_Tw_S)

        # Final prediction on new data

        yhat = aux1_Y.mean(axis=0) + j_2 + (aux2[0:20] - A_1) @ m.coef_.T

        pls_predict[i, d] = yhat
        
        aux1_a_predictors = np.delete(aux1, [df.columns.get_loc("Mw_inverse"), df.columns.get_loc("LCB"), df.columns.get_loc("SCB"), df.columns.get_loc("VNL"), df.columns.get_loc("VND")], 1)
        a = LinearRegression()
        a.fit(X=aux1_a_predictors, y=aux1_Y)
        
        aux2_a_predictors = np.delete(aux2, [df.columns.get_loc("Mw_inverse"), df.columns.get_loc("LCB"), df.columns.get_loc("SCB"), df.columns.get_loc("VNL"), df.columns.get_loc("VND")], 0).reshape(-1, 1).T
        ols_predict[i,d] = a.predict(aux2_a_predictors)

        # TODO No envelopes package
        # env_predict[i,d] = ...

p_ols = np.zeros((Dmax,))
p_pls = np.zeros((Dmax,))
# TODO No envelopes package
# p_env = np.zeros((Dmax,))

for d in range(Dmax):
    p_ols[d] = np.mean(np.sqrt(((Y - ols_predict[:,d])**2).mean(axis=1)))
    p_pls[d] = np.mean(np.sqrt(((Y - pls_predict[:,d])**2).mean(axis=1)))
    # TODO No envelopes package
    # p_env[d] = np.mean(np.sqrt(((Y - env_predict[:,d])**2).mean(axis=1)))

table6_2 = {
    "Method": ['OLS', 'PLS (P.Predictor)', 'Envelope (P.Predictor)', 'PLS (P.Response)', 'Envelope (P.Response)'],
    "No. components": ["not implemented"]*5,
    "predRMSE": ["not implemented"]*5
}

table6_2["No. components"][0] = "p1 = 20"
table6_2["predRMSE"][0] = p_ols[0]
table6_2["No. components"][1] = np.argmin(p_pls)+1
table6_2["predRMSE"][1] = np.min(p_pls)
# TODO No envelopes package
# table6_2["No. components"][2] = np.argmin(p_env)+!
# table6_2["predRMSE"][2] = np.min(p_env)


# Columns 5 and 6

# TODO No envelopes package
# fit_penv = np.zeros((X1.shape[0], Y.shape[1], X1.shape[1]))
# error_penv = np.zeros((X1.shape[0], Y.shape[1], X1.shape[1]))

# for d in range(Y.shape[1]):
#     for i in range(X1.shape[0]):
#         # Do envelope stuff


fit_prpls = np.zeros((X1.shape[0], Y.shape[1], X1.shape[1]))
error_prpls = np.zeros((X1.shape[0], Y.shape[1], X1.shape[1]))

# PPLS for d=1 we do it separately

for i in tqdm(range(X1.shape[0]), desc="Fitting PPLS for d=1"):

    Ythis = np.delete(Y, i, 0)
    X1this = np.delete(X1, i, 0)
    X2this = np.delete(X2, i, 0)

    lm = LinearRegression()
    lm.fit(X=X2this, y=Ythis)
    resY2 = Ythis - lm.predict(X2this)

    lm.fit(X=X2this, y=X1this)
    res12 = X1this - lm.predict(X2this)

    m = PLSRegression(n_components=1, scale=False)
    m.fit(X=resY2, Y=res12)

    A = Ythis @ (-m.x_rotations_)

    ss_X = np.hstack((X1this, X2this))
    ss = LinearRegression()
    ss.fit(X=ss_X, y=A)

    auxx = Ythis - X1this @ ss.coef_.T[0:20] @ (-m.x_rotations_.T)

    auxx1 = LinearRegression()
    auxx1.fit(X=X2this, y=auxx)
    auxx1 = auxx1.coef_

    fit_prpls[i, :, 0] = Ythis.mean(axis=0) + (X1[i] - X1this.mean(axis=0)) @ ss.coef_.T[0:20] @ (-m.x_rotations_.T) + (X2[i] - X2this.mean(axis=0)) @ auxx1.T

    error_prpls[i,:,0] = (Y[i] - fit_prpls[i,:,0])**2


# PPLS for the rest of the values for d
for d in tqdm(range(1, Y.shape[1]), desc="Fitting PPLS for 1<d<=ncol(Y)"):
    for i in range(X1.shape[0]):
        Ythis = np.delete(Y, i, 0)
        X1this = np.delete(X1, i, 0)
        X2this = np.delete(X2, i, 0)

        lm = LinearRegression()
        lm.fit(X=X2this, y=Ythis)
        resY2 = Ythis - lm.predict(X2this)

        lm.fit(X=X2this, y=X1this)
        res12 = X1this - lm.predict(X2this)

        m = PLSRegression(n_components=d+1, scale=False)
        m.fit(X=resY2, Y=res12)

        A = Ythis @ (-m.x_rotations_)

        ss_X = np.hstack((X1this, X2this))
        ss = LinearRegression()
        ss.fit(X=ss_X, y=A)

        auxx = Ythis - X1this @ ss.coef_.T[0:20] @ (-m.x_rotations_.T)

        auxx1 = LinearRegression()
        auxx1.fit(X=X2this, y=auxx)
        auxx1 = auxx1.coef_

        fit_prpls[i, :, d] = Ythis.mean(axis=0) + (X1[i] - X1this.mean(axis=0)) @ ss.coef_.T[0:20] @ (-m.x_rotations_.T) + (X2[i] - X2this.mean(axis=0)) @ auxx1.T

        error_prpls[i,:,d] = (Y[i] - fit_prpls[i,:,d])**2

# print(error_prpls)

p_pls = np.zeros((Y.shape[1],))
# TODO No envelopes package
# p_env = np.zeros((Y.shape[1],))

for d in range(Y.shape[1]):
    p_pls[d] = np.mean(np.sqrt(error_prpls[:,:,d].mean(axis=1)))
    # TODO No envelopes package
    # p_env[d] = np.mean(np.sqrt(error_penv[:,:,d].mean(axis=1)))

table6_2["No. components"][3] = np.argmin(p_pls)+1
table6_2["predRMSE"][3] = np.min(p_pls)
# TODO No envelopes package
# table6_2["No. components"][4] = np.argmin(p_env)+1
# table6_2["predRMSE"][4] = np.min(p_env)

print()
print("Table 6.2 (transposed):")
print(tabulate(table6_2, headers=table6_2.keys()))
print()