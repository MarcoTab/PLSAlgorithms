# Chapter 4
# Figure 4.1, Tables 4.1, 4.2

import warnings
warnings.filterwarnings("ignore")

import numpy as np
from sklearn.model_selection import LeaveOneOut
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import pandas as pd
from tabulate import tabulate
import scipy as sp

# Read in mussel data
df = pd.read_csv('data/chapter4/Mussel/musselslogMinus3.txt', sep=" ", header=0)
Y = df["logM"].to_numpy().reshape(-1, 1)
X = df.drop("logM", axis=1).to_numpy()

# Instantiate LOO CV
loocv = LeaveOneOut()

errors = []

total_ndims = 4

# For every number of dimensions we want to test
print("Fitting for mussel muscle mass", end="")
for i in range(1,total_ndims+1):
    print("." if (i) % 4 != 0 else "\b\b\b   \b\b\b", end="", flush=True)
    rmses = []
    # For every split in the dataset
    for (train_idxs, test_idx) in loocv.split(X):
        # Train a model on this dataset
        this_pls = PLSRegression(n_components=i, scale=False)

        this_pls.fit(X[train_idxs], Y[train_idxs])

        rmses.append(mean_squared_error(Y[test_idx], this_pls.predict(X[test_idx]), squared=False))
    
    errors.append(np.mean(rmses))

print()
print(f"Best n_components: {np.argmin(errors)+1}")
print()

# Train a model with the best number of components we saw.
best_pls = PLSRegression(n_components=np.argmin(errors)+1, scale=False)
best_pls.fit(X, Y)
preds = best_pls.predict(X)

line= np.linspace(1.5, 4, 100)

# Plot predictions against actual outcomes, with a line to indicate where a perfect predictor would be.
plt.scatter(Y, preds, facecolors="black")
plt.plot(line, line, '-r', color="blue")
plt.xlabel("Observed response, Y")
plt.ylabel("Predicted response")
plt.title("Figure 4.1")
plt.show()

# compute SD using the formulas from 
# Envelopes and partial least squares regression
# R. D. Cook, I. S. Helland, Z. Su
# Pages 851-877
# Number of pages 27
# Journal of the Royal Statistical Society. Series B: Statistical Methodology
# Volume 75
# Issue number 5
# Date Published - Nov 2013
shat = np.cov(np.concatenate((Y, X), axis=1), rowvar=False)[0, 1:].reshape(-1,1)
sts = (shat.T @ shat)[0,0]
sigma_0hat, r = np.linalg.qr(shat, mode="complete")
sigma_0hat = sigma_0hat[:,1:X.shape[1]]
SIhat = np.cov(X, rowvar=False)
Omegahat = ((shat.T @ SIhat @ shat) / sts)[0,0]
Omega0hat = (sigma_0hat.T @ SIhat @ sigma_0hat)
SIhat2 = (Omegahat/sts) * (shat@shat.T) + (sigma_0hat @ Omega0hat @ sigma_0hat.T)
sigmaY2hat = np.var(Y, ddof=1)
m = Y.shape[0] - 1

# Keep track of stds
sds=[]

# Compute Beta 1 SD
XN = np.array([[1,0,0,0]]).reshape(-1,1)
Vhat = ( (Omegahat**(-2)/m) * XN.T @ (SIhat2 * sigmaY2hat - shat@shat.T) @ XN )[0,0]
sds.append(round(np.sqrt(Vhat),4))

# Compute Beta 2 SD
XN = np.array([[0,1,0,0]]).reshape(-1,1)
Vhat = ( (Omegahat**(-2)/m) * XN.T @ (SIhat2 * sigmaY2hat - shat@shat.T) @ XN )[0,0]
sds.append(round(np.sqrt(Vhat),4))

# Compute Beta 3 SD
XN = np.array([[0,0,1,0]]).reshape(-1,1)
Vhat = ( (Omegahat**(-2)/m) * XN.T @ (SIhat2 * sigmaY2hat - shat@shat.T) @ XN )[0,0]
sds.append(round(np.sqrt(Vhat),4))

# Compute Beta 4 SD
XN = np.array([[0,0,0,1]]).reshape(-1,1)
Vhat = ( (Omegahat**(-2)/m) * XN.T @ (SIhat2 * sigmaY2hat - shat@shat.T) @ XN )[0,0]
sds.append(round(np.sqrt(Vhat),4))

# Get betas from trained model
coefs = [round(c[0], 3) for c in best_pls.coef_]

# For pretty print
SUB = str.maketrans("1234", "₁₂₃₄")
table = {
    "Estimate": coefs,
    "S.D.": sds
}
print("Table 4.1: PLS")
print(tabulate(table, showindex=["β" + str(i+1).translate(SUB) for i in range(4)], headers="keys"))
print()

SIhat2 = np.round(SIhat2, 3)

# For pretty print
table = {
    "log H": SIhat2[0],
    "log L": SIhat2[1],
    "log S": SIhat2[2],
    "log W": SIhat2[3]
}
print("Table 4.2: PLS")
print(tabulate(table, showindex=["log H", "log L", "log S", "log W"], headers="keys"))

