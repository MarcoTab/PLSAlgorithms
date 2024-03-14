import warnings
warnings.filterwarnings("ignore")

import numpy as np
# from sklearn.linear_model import LinearRegression
# from sklearn.cross_decomposition import PLSRegression
import pandas as pd
import pprint
# from tqdm import tqdm
# from tabulate import tabulate
import semopy

np.random.seed(5)

# Needs to be run with N=100 and N=1000
N = 10000

# Loadings are scaled os that errors have unit variances
Lrawx = np.array([.8, -.7, .6])
Lrawy = np.array([.8, -.7, .6])

Lambda = np.concatenate((Lrawx, np.zeros((6,)), Lrawy)).reshape(2,6).T

auxxx = Lrawx/np.sqrt(Lrawx.T @ Lrawx)
auxxy = Lrawy/np.sqrt(Lrawy.T @ Lrawy)


MODEL = """
X =~ x1 + x2 + x3
Y =~ y1 + y2 + y3
Y ~ X
"""

# Vary the population beta

BETAS = np.arange(200/1500, (1000+1)/1500, 20/1500)

nrep = 10

nose = np.empty((BETAS.shape[0], 7, nrep))
nose[:] = np.NaN

sigmacond = np.eye(6)


for h in range(nrep):
    for j in range(BETAS.shape[0]):
        Beta = BETAS[j]

        eta = np.random.multivariate_normal(mean=np.zeros((2,)), cov=np.array([[1, Beta], [Beta, 1]]), size=(N,))

        print(eta.mean(axis=0))
        print(np.cov(eta, rowvar=False, ddof=1))

        data = eta @ Lambda.T + np.random.multivariate_normal(mean=np.zeros((6,)), cov=sigmacond, size=(N,))

        dataframe = pd.DataFrame(
            data=data,
            columns=["x1", "x2", "x3", "y1", "y2", "y3"]
        )

        dataframe.to_csv("example.csv", index=False)

        # CB | SEM Estimates

        ml = semopy.Model(MODEL)
        ml.fit(dataframe, obj="FIML")
        # pprint.pprint(ml.inspect(mode="mx", std_est=True))
        # TODO Unfinished script
        exit()