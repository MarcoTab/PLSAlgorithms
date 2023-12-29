import warnings
warnings.filterwarnings("ignore")

import numpy as np
from sklearn.linear_model import LinearRegression, LassoCV
from sklearn.cross_decomposition import PLSRegression
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

np.random.seed(5)

df = pd.read_csv("data/chapter6/DepuredNewDataGovernance_withLagGDP_22.csv", header=0)

df = df.iloc[:, 1:]

aux_2003 = df[df.year >= 2003].to_numpy()

aux_2003 = np.delete(aux_2003, 1, 0)

n = aux_2003.shape[0]

X1 = aux_2003[:, [3,4,5,6,7,8,10,11,12,14,15,16,18,19,20,22,23,24,26,27,28,30,31,32]]

X2 = aux_2003[:, [33,34,35,36,37,38,40] + [i for i in range(42, 52+1)] + [54]]

y = aux_2003[:, 0]

lasso_prediction = [0] * n
pls_prediction = [0] * n

for i in tqdm(range(n), desc="Fitting Lasso and PPLS models..."):
    # Training data
    X11_train = np.delete(X1, i, 0)
    X22_train = np.delete(X2, i, 0)
    y_train = np.delete(y, i, 0)

    # Testing data
    X11_test = X1[i,:].reshape(-1, 1).T
    X22_test = X2[i,:].reshape(-1, 1).T
    y_test = y[i]

    #### PPLS
    A = LinearRegression()
    A.fit(X=X22_train, y=X11_train)

    res_train = X11_train - A.predict(X22_train)

    A1 = A.predict(X22_test)

    y_train_center = y_train - y_train.mean()

    # Perform PLS on residuals in training
    m = PLSRegression(n_components=8, scale=False)
    m.fit(X=res_train, Y=y_train_center)

    j = LinearRegression()
    j.fit(X=X22_train, y=y_train_center)

    j_2 = j.predict(X22_test)

    yhat = np.mean(y_train) + j_2.item() + ((X11_test - A1) @ (m.coef_).T).item()

    pls_prediction[i] = yhat

    #### LASSO
    covariates = np.concatenate((X11_train, X22_train), axis=1)
    z = y_train

    cv_model = LassoCV(cv=10)
    # LassoCV Automatically chooses the best model
    cv_model.fit(covariates, z)

    covariates_test = np.concatenate((X11_test, X22_test), axis=1)

    # Use the Lasso regression model to predict response value
    yhat = cv_model.predict(covariates_test)

    lasso_prediction[i] = yhat


plt.scatter(x=pls_prediction, y=y, alpha=0.7, c="black")
plt.xlabel("Leave-one-out fitted values from PPLS-7")
plt.ylabel("Y")
plt.show()

plt.scatter(x=lasso_prediction, y=y, alpha=0.7, c="black")
plt.xlabel("Leave-one-out fitted values from LASSO")
plt.ylabel("Y")
plt.show()
