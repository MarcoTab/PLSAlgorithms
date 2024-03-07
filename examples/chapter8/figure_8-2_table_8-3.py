import warnings
warnings.filterwarnings("ignore")

import numpy as np
from sklearn.metrics import confusion_matrix, mean_squared_error
from tabulate import tabulate
import matplotlib.pyplot as plt
import pandas as pd

from quadratic_with_cross_validation import quadratic_discrimination
from pls_matrices import pls_matrices
import functools as ft

np.random.seed(5)

MIR_Fruit_purees = pd.read_csv("data/chapter8/MIR_Fruit_purees.csv", header=None)
# print(MIR_Fruit_purees)
Yaux = MIR_Fruit_purees.iloc[0, 1:].to_numpy()
Y = np.empty_like(Yaux)
Y[Yaux == "NON-Strawberry"] = 0
Y[Yaux == "Strawberry"] = 1

Y = Y.reshape(-1,1)
Y = Y.astype('int')

X = MIR_Fruit_purees.iloc[1:, 1:].to_numpy().T
X = np.array([[float(i) for i in x] for x in X])

Dmax = 20
K = 10
h = np.linspace(0, 1, 6)

quad_1 = np.zeros((Dmax, len(h)))
quad_2 = np.zeros((Dmax, len(h)))

for d in range(1, Dmax+1):
    for s in range(len(h)):

        print(f"{d}, {s+1} / {Dmax}, {len(h)}")
        alp = h[s]

        aux = quadratic_discrimination(X, Y, K, d, alp)

        li_1 = aux['predict_quad_1']
        li_2 = aux['predict_quad_2']

        quad_1[d-1, s] = mean_squared_error(Y, li_1)
        quad_2[d-1, s] = mean_squared_error(Y, li_2)

M1 = np.where(quad_1 == np.min(quad_1))
if len(M1[0]) <= 1:
    M1 = M1[0][0]+1, M1[1][0]
else:
    M1 = M1[0]
    M1 = M1[0]+1, M1[1]

M2 = np.where(quad_2 == np.min(quad_2))
if len(M2[0]) <= 1:
    M2 = M2[0][0]+1, M2[1][0]
else:
    M2 = M2[0]
    M2 = M2[0]+1, M2[1]

# Run to get data for table 8.3
# AN-Q1
print("Calculating LOOCV for AN-Q1 (may take a long time)")
d = M1[0]
alp = h[M1[1]]
aux = quadratic_discrimination(X, Y, X.shape[0], d, alp, loocv=True)['predict_quad_1']
cm=confusion_matrix(y_true=Y, y_pred=aux)
anq1 = cm.diagonal().sum() / cm.sum()

# AN-Q2
print("Calculating LOOCV for AN-Q2 (may take a long time)")
d = M2[0]
alp = h[M2[1]]
aux = quadratic_discrimination(X, Y, X.shape[0], d, alp, loocv=True)['predict_quad_2']
cm=confusion_matrix(y_true=Y, y_pred=aux)
anq2 = cm.diagonal().sum() / cm.sum()

table = {
    "Dataset": [
        "Fruit"
    ],
    "AN-Q1": [
        round(100*(anq1), 2)
    ],
    "AN-Q2": [
        round(100*(anq2), 2)
    ]
}

print("Table 8.3")
print(tabulate(table, headers=table.keys()))

## Figure 8.2
holdout=X
y_holdout=Y

S_X = np.cov(holdout, rowvar=False, ddof=1)
uY = np.array([int(x) for x in np.unique(y_holdout)])

S_k = np.array([x for x in map(lambda k: np.cov(holdout[(y_holdout==k).flatten()], rowvar=False, ddof=1), uY)])

cov_mean = np.array(ft.reduce(lambda a, b: a+b, [x for x in map(lambda k: S_k[k] * np.mean(y_holdout == k), uY)]))

cov_ss1 = np.array(ft.reduce(lambda a, b: a+b, [x for x in map(lambda k: (S_X - S_k[k]) @ (S_X - S_k[k]).T * np.mean(y_holdout == k), uY)]))

cov_ss2 = np.array(ft.reduce(lambda a, b: a+b, [x for x in map(lambda k: (cov_mean - S_k[k]) @ (cov_mean - S_k[k]).T * np.mean(y_holdout == k), uY)]))

d = M1[0]
alp = h[M1[1]]

dif_d_delta2 = (S_X - cov_mean) @ (S_X - cov_mean).T
A_al1 = cov_ss1
AUX = pls_matrices(A_al1, cov_mean, d)

d = M2[0]
alp = h[M2[1]]

A_al_para = alp * dif_d_delta2 + (1-alp)*cov_ss2
AUX_para = pls_matrices(A_al_para, cov_mean, d)

proj_data_2 = holdout @ AUX_para['Gamma']
proj_data_1 = holdout @ AUX['Gamma']

print(proj_data_1[:, 0].shape)

plt.title("Figure 8.2.a: AN-Q2 Projections (Fruit)")
plt.xlabel("First projection using AN-Q2")
plt.ylabel("Second projection using AN-Q2")
plt.scatter(x=-proj_data_2[(Y==0).flatten()][:, 0], y=-proj_data_2[(Y==0).flatten()][:, 1], label="NON-Strawberry", marker="o", color="red", alpha=0.7)
plt.scatter(x=-proj_data_2[(Y==1).flatten()][:, 0], y=-proj_data_2[(Y==1).flatten()][:, 1], label="Strawberry", marker="^", color="cyan", alpha=0.7)
plt.legend(loc='best')
plt.show()

plt.title("Figure 8.2.b: AN-Q1 Projections (Fruit)")
plt.xlabel("First projection using AN-Q1")
plt.ylabel("Second projection using AN-Q1")
plt.scatter(x=proj_data_1[(Y==0).flatten()][:, 0], y=proj_data_1[(Y==0).flatten()][:, 1], label="NON-Strawbery", marker="o", color="red", alpha=0.7)
plt.scatter(x=proj_data_1[(Y==1).flatten()][:, 0], y=proj_data_1[(Y==1).flatten()][:, 1], label="Strawberry", marker="^", color="cyan", alpha=0.7)
plt.legend(loc='best')
plt.show()
