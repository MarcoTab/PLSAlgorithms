import warnings
warnings.filterwarnings("ignore")

import numpy as np
from sklearn.metrics import confusion_matrix, mean_squared_error
from sklearn.linear_model import LinearRegression, LassoCV
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression
from tabulate import tabulate
import matplotlib.pyplot as plt
import pandas as pd
from progress.bar import ChargingBar

np.random.seed(5)

df = pd.read_csv("data/chapter6/DepuredNewDataGovernance_withLagGDP_22.csv", header=0)

df = df.iloc[:, 1:]

aux_2003 = df[df.year >= 2003].to_numpy()
aux_2008 = df[df.year >= 2008].to_numpy()
aux_2010 = df[df.year >= 2010].to_numpy()
aux_2010_2014 = df[df.year >= 2010]
aux_2010_2014 = aux_2010_2014[aux_2010_2014.year <= 2014].to_numpy()
aux_2013 = df[df.year >= 2013].to_numpy()
aux_2014 = df[df.year >= 2014].to_numpy()
aux_2015 = df[df.year >= 2015].to_numpy()

table_lines = [aux_2003, aux_2008, aux_2010, aux_2010_2014, aux_2013, aux_2014, aux_2015]

table61 = {
    'n': [161, 109, 86, 53, 51, 41, 32],
    'Period': ["2003-2018", "2008-2018", "2010-2018", "2010-2014", "2013-2018", "2014-2018", "2015-2018"],
    'Countries': [12, 12, 12, 12, 10, 9, 9],
    'PPLS.q1': [np.nan]*7,
    'PPCA.q1': [np.nan]*7,
    'PPLS.1': [np.nan]*7,
    'PPCA.1': [np.nan]*7,
    'LASSO': [np.nan]*7,
    'FULL': [np.nan]*7,
    # 'PPENV.q1': [np.nan]*7,
    # 'PPENV.1': [np.nan]*7
    # TODO Implement an envelope package for Python 😬
}

index = 0

print("Generating Table 6.1")

for this_aux in table_lines:

    this_aux = np.delete(this_aux, 1, 0)

    n = this_aux.shape[0]

    X1 = this_aux[:, [3,4,5,6,7,8,10,11,12,14,15,16,18,19,20,22,23,24,26,27,28,30,31,32]]

    X2 = this_aux[:, [33,34,35,36,37,38,40,42,43,44,45,46,47,48,49,50,51,52,54]]

    y = this_aux[:, 0]

    mm = X1.shape[1]-2

    pls_predict = np.empty((n, mm))
    pls_predict.fill(float('inf'))
    MSE_PLS = np.empty((n, mm))
    MSE_PLS.fill(float('inf'))

    pca_predict = np.empty((n, mm))
    pca_predict.fill(float('inf'))
    MSE_PCA = np.empty((n, mm))
    MSE_PCA.fill(float('inf'))

    ols_predict = np.empty((n, mm))
    ols_predict.fill(float('inf'))
    MSE_OLS = np.empty((n, mm))
    MSE_OLS.fill(float('inf'))

    lasso_predict = np.empty((n, mm))
    lasso_predict.fill(float('inf'))
    MSE_LASSO = np.empty((n, mm))
    MSE_LASSO.fill(float('inf'))

    # TODO No envelope package
    # env_predict = np.empty((n, mm))
    # MSE_EPLS = np.empty((n, mm))

    progbar = ChargingBar(f"Row {index+1}", max=mm*n, suffix="%(percent)d%%: %(index)d / %(max)d | ~%(eta)ds")

    for d in range(mm):
        for i in range(n):
            progbar.next()
            X11_train = np.delete(X1, i, 0)
            X22_train = np.delete(X2, i, 0)
            y_train = np.delete(y, i)
            X11_test = X1[i].reshape(-1, 1).T
            X22_test = X2[i].reshape(-1, 1).T
            y_test = y[i]

            ###### PPLS
            # Regression of the continuous X which we do not reduce
            A = LinearRegression()
            A.fit(X=X22_train, y=X11_train)

            # Save residuals
            res_train = X11_train - A.predict(X22_train)

            # Predict on new data
            A_1 = A.predict(X22_test)
            
            # Center training response variable
            y_train_center = y_train - np.mean(y_train)

            # Perform PLS on residuals in training
            m = PLSRegression(n_components=d+1, scale=False)
            m.fit(X=res_train, Y=y_train_center)

            # Perform regression on non-reduced variables
            j = LinearRegression()
            j.fit(X=X22_train, y=y_train_center)

            # Predict non-reduced on new data
            j_2 = j.predict(X22_test)

            # Final prediction on new data

            yhat = np.mean(y_train) + j_2.item() + ((X11_test - A_1) @ (m.coef_).T).item()

            pls_predict[i,d] = yhat
            MSE_PLS[i,d] = yhat - y_test


            #### LASSO
            covariates = np.concatenate((X11_train, X22_train), axis=1)
            # Find and fit best model
            cv_model = LassoCV(cv=10)    
            cv_model.fit(covariates, y_train)

            covariates_test = np.concatenate((X11_test, X22_test), axis=1)

            # Use the Lasso regression model to predict response value
            yhat = cv_model.predict(covariates_test)
            lasso_predict[i,d] = yhat
            MSE_LASSO[i, d] = yhat-y_test

            ####### PPCA
            # Perform PCR on training residuals
            m = PCA(n_components=d+1)
            m.fit(X=res_train)
            
            pcatransformed = m.transform(res_train)

            mod = LinearRegression()
            mod.fit(pcatransformed, y_train_center)

            # Perform regression on non-reduced variables
            j = LinearRegression()
            j.fit(X=X22_train, y=y_train_center)

            # Predict non-reduced on new data
            j_2 = j.predict(X22_test)

            # Final prediction on new data
            yhat = np.mean(y_train) + j_2.item() + ((X11_test - A_1) @ (m.components_.T @ mod.coef_)).item()

            

            pca_predict[i,d] = yhat
            MSE_PCA[i,d] = yhat-y_test

            if (X11_train.shape[0] > (X1.shape[1] + X2.shape[1])):
                train_predictors = np.concatenate((X11_train, X22_train), axis=1)

                fit=LinearRegression()
                fit.fit(X=train_predictors, y=y_train)

                # Predicting
                test_predictors = np.concatenate((X11_test, X22_test), axis=1)
                yhat = fit.predict(test_predictors)

                ols_predict[i,d] = yhat
                MSE_OLS[i,d] = yhat - y_test

            # TODO Envelope package
            # if (X11_train.shape[0] > this_aux.shape[1]):

    progbar.finish()

    print()
    mses_pls = np.mean(MSE_PLS**2, axis=0)
    d_PLS = np.where(mses_pls == np.min(mses_pls))
    table61['PPLS.q1'][index] = mses_pls[d_PLS[0][0]]
    table61['PPLS.1'][index] = mses_pls[0]

    mse_lasso = np.mean(np.mean(MSE_LASSO**2, axis=0))
    table61['LASSO'][index] = mse_lasso

    mses_pca = np.mean(MSE_PCA**2, axis=0)
    d_PCA = np.where(mses_pca == np.min(mses_pca))
    table61['PPCA.q1'][index] = mses_pca[d_PCA[0][0]]
    table61['PPCA.1'][index] = mses_pca[0]

    mse_ols = np.mean(MSE_OLS**2)
    table61['FULL'][index] = mse_ols

    index+=1


print(tabulate(table61, headers='keys'))

