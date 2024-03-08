# Chapter 7 Coffee Data Examples

import warnings
warnings.filterwarnings("ignore")

import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import confusion_matrix
from sklearn.linear_model import LinearRegression
from tabulate import tabulate
from sliced import SlicedInverseRegression
import matplotlib.pyplot as plt

from lib.lineal_with_R import lineal_junto_small_R
# import pls_matrices


## TABLE 7.1
np.random.seed(180)

FITR_Spectra_instant_coffee = np.genfromtxt("data/chapter7/coffee.csv", delimiter=",")

Y = (FITR_Spectra_instant_coffee[1, 1:57].T - 1).reshape(-1, 1)
X = (FITR_Spectra_instant_coffee[3:289, 1:57].T)


# Choose d with k=10 folds
Dmax = 15
K = 10

linear_error_pfc = [-1]*Dmax
linear_error_pls = [-1]*Dmax
linear_error_iso = [-1]*Dmax

for d in range(1, Dmax+1):

    aux = lineal_junto_small_R(X, Y, K, d)

    li = aux['predict_linear_pls']
    li_pfc = aux['predict_linear_pfc']
    li_iso = aux['predict_linear_iso']

    linear_error_pls[d-1] = np.mean((li - Y)**2)
    linear_error_pfc[d-1] = np.mean((li_pfc - Y)**2)
    linear_error_iso[d-1] = np.mean((li_iso - Y)**2)

d_pls = np.argmin(linear_error_pls)+1
d_pfc = np.argmin(linear_error_pfc)+1
d_iso = np.argmin(linear_error_iso)+1

aux_pls = lineal_junto_small_R(X, Y, X.shape[0], d_pls, loocv=True)['predict_linear_pls']
aux_pfc = lineal_junto_small_R(X, Y, X.shape[0], d_pls, loocv=True)['predict_linear_pfc']
aux_iso = lineal_junto_small_R(X, Y, X.shape[0], d_pls, loocv=True)['predict_linear_iso']

cm = confusion_matrix(y_true=Y, y_pred=aux_pls-1)
plserr = 100*(cm.diagonal().sum() / cm.sum())

cm = confusion_matrix(y_true=Y, y_pred=aux_iso-1)
isoerr = 100*(cm.diagonal().sum() / cm.sum())

cm = confusion_matrix(y_true=Y, y_pred=aux_pfc-1)
pls_pfcerr = 100*(cm.diagonal().sum() / cm.sum())

table = {
    "Dataset": ["Coffee"],
    "PLS": [f"{100-plserr:.2f}"],
    "ISO": [f"{100-isoerr:.2f}"],
    "PLS+PFC": [f"{100-pls_pfcerr:.2f}"]
}
# print(table)

print("Table 7.1")
print(tabulate(table, headers=table.keys()))
print()

# Figure 7.2

Dmax = 15
K = X.shape[0]
AA = []

for d in range(1,Dmax+1):
    aux_pls = lineal_junto_small_R(X, Y, X.shape[0], d, loocv=True)['predict_linear_pls']

    cm = confusion_matrix(y_true=Y, y_pred=aux_pls-1)
    plserr = 100-100*(cm.diagonal().sum() / cm.sum())

    AA.append(plserr)

XX = np.arange(1, Dmax+1)

plt.plot(XX, AA, color="black")
plt.scatter(XX, AA, color="black")
plt.xlabel("Number of components")
plt.ylabel("Accuracy")
plt.title("Figure 7.2")
plt.show()


## Figure 7.3 a
# PLS
d = d_pls
mod = PLSRegression(n_components=d, scale=False)
mod.fit(X, Y)
AUX = mod.x_loadings_

proj_data = X @ AUX

plt.title("Figure 7.3.a: PLS Projections")
plt.xlabel("First PLS Projection")
plt.ylabel("Second PLS Projection")
plt.scatter(x=proj_data[(Y==0).flatten()][:, 0], y=proj_data[(Y==0).flatten()][:, 1], label="robusta", marker="^", color="blue")
plt.scatter(x=proj_data[(Y!=0).flatten()][:, 0], y=proj_data[(Y!=0).flatten()][:, 1], label="arabica", color="green", marker="o")
plt.legend(loc='best')
plt.show()


## Figure 7.3 b
# PFC
d=d_pfc
mod = PLSRegression(n_components=d, scale=False)
mod.fit(X, Y)
AUX = mod.x_loadings_

proj_data = X @ AUX

s0 = SlicedInverseRegression()
s0.fit(X=proj_data, y=Y.flatten())
vecs = s0.directions_[np.argmax(s0.eigenvalues_)].reshape(-1,1)

proj_data_pfc = proj_data @ vecs

plt.title("Figure 7.3.b: PLS+PFC Projections")
plt.xlabel("PFC Projection")
plt.ylabel("Number of cases")

plt.hist(x=proj_data_pfc[Y==0], alpha=0.5, label='robusta', color='blue', edgecolor="black")
plt.hist(x=proj_data_pfc[Y!=0], alpha=1.0, label='arabica', color='green', edgecolor="black")
plt.legend(loc='best')
plt.show()


# Figure 7.3.c
d=Dmax
# ISO
reg = LinearRegression()
# Inverse regression
reg.fit(Y, X)
residuals = X - reg.predict(Y) 

_, evecs = np.linalg.eigh(np.cov(residuals, rowvar=False))

evecs = evecs[:, ::-1]
evecs = evecs[:, :d]

proj_data_iso = X @ evecs

plt.title("Figure 7.3.c: Isotropic (ISO) Projections")
plt.xlabel("First ISO Projection")
plt.ylabel("Second ISO Projection")
plt.scatter(x=proj_data_iso[(Y==0).flatten()][:, 0], y=proj_data_iso[(Y==0).flatten()][:, 1], label="robusta", marker="^", color="blue")
plt.scatter(x=proj_data_iso[(Y!=0).flatten()][:, 0], y=proj_data_iso[(Y!=0).flatten()][:, 1], label="arabica", color="green", marker="o")
plt.legend(loc='best')
plt.show()