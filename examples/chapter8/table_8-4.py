import warnings
warnings.filterwarnings("ignore")

import numpy as np
from sklearn.metrics import confusion_matrix, mean_squared_error
from tabulate import tabulate

from quadratic_with_cross_validation import quadratic_discrimination
from pls_matrices import pls_matrices
import functools as ft

np.random.seed(5)

FITR_Spectra_instant_coffee = np.genfromtxt("data/chapter8/coffee.csv", delimiter=",")

Y = (FITR_Spectra_instant_coffee[1, 1:57].T - 1).reshape(-1, 1)
X = (FITR_Spectra_instant_coffee[3:289, 1:57].T)

Dmax = 12
K = 10
h = np.linspace(0, 1, 6)

quad_1 = np.zeros((Dmax, len(h)))
quad_2 = np.zeros((Dmax, len(h)))

print("Performing crossvalidation on coffee data.")
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

# Table 8.4
d = M1[0]
alp = M1[1]

print("Calculating LOOCV for AN-Q1 (may take a long time)")
aux = quadratic_discrimination(X, Y, X.shape[0], d, alp, loocv=True)['predict_quad_1']
cm=confusion_matrix(y_true=Y, y_pred=aux)
canq1 = cm.diagonal().sum() / cm.sum()

d = M2[0]
alp = M2[1]

print("Calculating LOOCV for AN-Q2 (may take a long time)")
aux = quadratic_discrimination(X, Y, X.shape[0], d, alp, loocv=True)['predict_quad_2']
cm=confusion_matrix(y_true=Y, y_pred=aux)
canq2 = cm.diagonal().sum() / cm.sum()


oil = np.genfromtxt("data/chapter8/FTIR_Spectra_olive_oils.csv", delimiter=",")

Y = (oil[1, 1:121].T - 1).reshape(-1, 1)
X = (oil[3:571, 1:121].T)

Dmax = 14
K = 10
h = np.linspace(0, 1, 6)

quad_1 = np.zeros((Dmax, len(h)))
quad_2 = np.zeros((Dmax, len(h)))

print("Performing crossvalidation on olive oil data.")
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

d = M1[0]
alp = M1[1]

print("Calculating LOOCV for AN-Q1 (may take a long time)")
aux = quadratic_discrimination(X, Y, X.shape[0], d, alp, loocv=True)['predict_quad_1']
cm=confusion_matrix(y_true=Y, y_pred=aux)
oanq1 = cm.diagonal().sum() / cm.sum()

d = M2[0]
alp = M2[1]

print("Calculating LOOCV for AN-Q2 (may take a long time)")
aux = quadratic_discrimination(X, Y, X.shape[0], d, alp, loocv=True)['predict_quad_2']
cm=confusion_matrix(y_true=Y, y_pred=aux)
oanq2 = cm.diagonal().sum() / cm.sum()


# Table 8.4 (FULL)
table = {
    "Dataset": [
        "Coffee", "Olive Oil"
    ],
    "AN-Q1": [
        round(100*(canq1), 2), round(100*(oanq1), 2)
    ],
    "AN-Q2": [
        round(100*(canq2), 2), round(100*(oanq2), 2)
    ]
}

print("Table 8.4")
print(tabulate(table, headers=table.keys()))