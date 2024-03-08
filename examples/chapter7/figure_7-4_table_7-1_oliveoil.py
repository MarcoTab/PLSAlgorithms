import warnings
warnings.filterwarnings("ignore")

import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import confusion_matrix
from sklearn.linear_model import LinearRegression
from tabulate import tabulate
from sliced import SlicedInverseRegression
import matplotlib.pyplot as plt

from lib.lineal_with_R_more_classes import lineal_junto_small_R_more_classes

## TABLE 7.1
np.random.seed(180)

FITR_Spectra_olive_oils = np.genfromtxt("data/chapter7/FTIR_Spectra_olive_oils.csv", delimiter=",")

Y = (FITR_Spectra_olive_oils[1, 1:121].T - 1).reshape(-1, 1)
X = (FITR_Spectra_olive_oils[3:571, 1:121].T)

# Choose d with k=5 folds
Dmax = 33
K = 5

linear_error_pfc = [-1]*Dmax
linear_error_pls = [-1]*Dmax
linear_error_iso = [-1]*Dmax

for d in range(1, Dmax+1):

    aux = lineal_junto_small_R_more_classes(X, Y, K, d)

    li = aux['predict_linear_pls']
    li_pfc = aux['predict_linear_pfc']
    li_iso = aux['predict_linear_iso']

    linear_error_pls[d-1] = np.mean((li - Y)**2)
    linear_error_pfc[d-1] = np.mean((li_pfc - Y)**2)
    linear_error_iso[d-1] = np.mean((li_iso - Y)**2)

d_pls = np.argmin(linear_error_pls)+1
d_pfc = np.argmin(linear_error_pfc)+1
d_iso = np.argmin(linear_error_iso)+1

aux_pls = lineal_junto_small_R_more_classes(X, Y, X.shape[0], d_pls, loocv=True)['predict_linear_pls']
aux_pfc = lineal_junto_small_R_more_classes(X, Y, X.shape[0], d_pls, loocv=True)['predict_linear_pfc']
aux_iso = lineal_junto_small_R_more_classes(X, Y, X.shape[0], d_pls, loocv=True)['predict_linear_iso']

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

print("Table 7.1")
print(tabulate(table, headers=table.keys()))
print()

## Figure 7.4 a
# PLS
d = d_pls
mod = PLSRegression(n_components=d, scale=False)
mod.fit(X, Y)
AUX = mod.x_loadings_

proj_data = X @ AUX

plt.title("Figure 7.4.a: PLS Projections")
plt.xlabel("First PLS Projection")
plt.ylabel("Second PLS Projection")
plt.scatter(x=proj_data[(Y==0).flatten()][:, 0], y=proj_data[(Y==0).flatten()][:, 1], label="Greece", marker="o")
plt.scatter(x=proj_data[(Y==1).flatten()][:, 0], y=proj_data[(Y==1).flatten()][:, 1], label="Italy", marker="^")
plt.scatter(x=proj_data[(Y==2).flatten()][:, 0], y=proj_data[(Y==2).flatten()][:, 1], label="Portugal", marker="s")
plt.scatter(x=proj_data[(Y==3).flatten()][:, 0], y=proj_data[(Y==3).flatten()][:, 1], label="Spain", marker="+")
plt.legend(loc='best')
plt.show()


## Figure 7.4 b
# PFC
d=d_pfc
mod = PLSRegression(n_components=d, scale=False)
mod.fit(X, Y)
AUX = mod.x_loadings_

proj_data = X @ AUX

s0 = SlicedInverseRegression()
s0.fit(X=proj_data, y=Y.flatten())
vecs = s0.directions_[np.argsort(s0.eigenvalues_)[::-1][:2]].T

proj_data_pfc = proj_data @ vecs

plt.title("Figure 7.4.b: PLS+PFC Projections")
plt.xlabel("First PFC Projection")
plt.ylabel("Second PFC Projection")
plt.scatter(x=proj_data_pfc[(Y==0).flatten()][:, 0], y=proj_data_pfc[(Y==0).flatten()][:, 1], label="Greece", marker="o")
plt.scatter(x=proj_data_pfc[(Y==1).flatten()][:, 0], y=proj_data_pfc[(Y==1).flatten()][:, 1], label="Italy", marker="^")
plt.scatter(x=proj_data_pfc[(Y==2).flatten()][:, 0], y=proj_data_pfc[(Y==2).flatten()][:, 1], label="Portugal", marker="s")
plt.scatter(x=proj_data_pfc[(Y==3).flatten()][:, 0], y=proj_data_pfc[(Y==3).flatten()][:, 1], label="Spain", marker="+")
plt.legend(loc='best')
plt.show()

# Figure 7.4.c
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

plt.title("Figure 7.4.c: Isotropic (ISO) Projections")
plt.xlabel("First ISO Projection")
plt.ylabel("Second ISO Projection")
plt.scatter(x=proj_data_iso[(Y==0).flatten()][:, 0], y=proj_data_iso[(Y==0).flatten()][:, 1], label="Greece", marker="o")
plt.scatter(x=proj_data_iso[(Y==1).flatten()][:, 0], y=proj_data_iso[(Y==1).flatten()][:, 1], label="Italy", marker="^")
plt.scatter(x=proj_data_iso[(Y==2).flatten()][:, 0], y=proj_data_iso[(Y==2).flatten()][:, 1], label="Portugal", marker="s")
plt.scatter(x=proj_data_iso[(Y==3).flatten()][:, 0], y=proj_data_iso[(Y==3).flatten()][:, 1], label="Spain", marker="+")
plt.legend(loc='best')
plt.show()