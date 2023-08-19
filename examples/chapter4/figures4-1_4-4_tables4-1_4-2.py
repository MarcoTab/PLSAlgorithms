# Chapter 4
# Figures 4.1, 4.4, Tables 4.1, 4.2
#### READ BEFORE RUNNING ALL THE WAY THROUGH
# There is a very long simulation that runs as part of this script. 
# To avoid having to run this simulation again and again, please uncomment the code below (you might have to scroll a
# long ways, or search "uncomment") to save the data to your working directory.
# You must have the joblib package installed to do so.

# Some numbers will be different between the book and this script, since the random number generation for Python and R are different.


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
df = pd.read_csv('data/chapter4/musselslogMinus3.txt', sep=" ", header=0)
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




print("Simulating", end="", flush=True)

def estimate_beta(X, y, pls=True):
    # Compute the beta coefficient for the regression
    # y = E(y) + beta^T (X-E(X)) + error
    # y in R, X in R^p, Sigma_x=Cov(X), sigma_xy=cov(X,y)
    # assuming that the PLS model with d=1 is true, i.e.
    # beta = sigma_xy (sigma_xy^T Sigma_x sigma_xy)^{-1} (sigma_xy^T sigma_xy)

    n = X.shape[0]
    degf = n
    p = X.shape[1]
    mX = np.mean(X, axis=0)
    my = np.mean(y)
    Xc = X - mX
    yc = y - my
    sXX = Xc.T @ Xc / degf
    sXy = Xc.T @ yc / degf
    syy = np.sum(yc ** 2) / degf
    
    if pls:
        aux1 = sXy.T @ sXy
        aux2 = sXy.T @ (sXX @ sXy)
        b_pls = np.linalg.inv(aux2) @ aux1 * sXy
        b_pls = np.asarray(b_pls).flatten()

    return b_pls

#############################
#############################
#simulation to get Figure 4.4
#############################
#############################
np.random.seed(180)
reps = 500
p_vec = np.array([2**i for i in range(3,12)])
pm = p_vec[-1]
pow = 1/2
s = pow
a = 1
ovec = np.array([1] + [0]*pm)

###### Start commenting here to avoid running simulation again
if (pow > 0):
    ovec = np.array([1]*round(p_vec[0]**pow) + [0]*(p_vec[0] - round(p_vec[0]**pow)))

    for k in range(1, len(p_vec)):
        b = len(ovec)
        p = p_vec[k]
        cnz = np.sum(ovec)
        tnz = round(p**pow) 
        nnz = tnz - cnz
        newvec = np.array([1]*nnz + [0]*(b-nnz))
        ovec = np.concatenate((ovec, newvec))

tau = 1/2
pe_pls1 = np.empty((len(p_vec), reps))
pe_pls2 = np.empty((len(p_vec), reps))
nm = pm/2
n_vec = p_vec/2

for i in range(reps):
    sigma_xym = (ovec*np.random.normal(0, 1, pm)).reshape(-1, 1)
    sigma_0, _ = np.linalg.qr(sigma_xym, mode="complete")
    sigma_0 = sigma_0[:, 1:pm]
    for k in range(len(p_vec)):
        print(f". {i}/{reps}", end="", flush=True)
        p = p_vec[k]
        n = int(p/2)
        m = n-1
        Cp1 = p
        Cp2 = p
        sigma_xymp = sigma_xym[:p]
        sTs = (sigma_xymp.T @ sigma_xymp).item()
        Omega = sTs
        sigmaY2 = tau ** 2 + sTs * Omega**(-1)
        sigma_0, _ = np.linalg.qr(sigma_xymp, mode="complete")
        sigma_0 = sigma_0[:, 1:pm]
        Omega0 = np.eye(p-1)
        Sigma = (Omega/sTs) * (sigma_xymp @ sigma_xymp.T) + sigma_0 @ (Omega0 @ sigma_0.T)
        out_eigenvals, out_eigenvecs = np.linalg.eigh(Sigma)
        Sigma_sqrt = out_eigenvecs @ np.diag(np.sqrt(out_eigenvals)) @ out_eigenvecs.T
        X = np.random.normal(0, 1, (n,p)) @ Sigma_sqrt
        true_beta = Omega**(-1) * sigma_xymp
        y = np.array([np.random.normal((X @ true_beta)[i], tau, 1) for i in range(n)])
        Xnewc = sigma_xymp
        b_pls = estimate_beta(X=X, y=y)
        bias = (tau**2 - sTs/Omega) * (Cp1/(m*sTs)) - sigmaY2 * (Cp1 ** 2) / (m**2*Omega*sTs) - sigmaY2 * Cp2 / (m*Omega*sTs)
        V = ((Omega**(-2)/n)*Xnewc.T @ (Sigma*sigmaY2-sigma_xymp @ sigma_xymp.T) @ Xnewc).item()
        pe_pls1[k, i] = V**(-1/2) * (b_pls.T - true_beta.T) @ Xnewc
        pe_pls2[k, i] = V**(-1/2) * (b_pls.T - true_beta.T * (1+bias)) @ Xnewc 

        deletelen = len(f" {i}/{reps}")
        print(("\b"*deletelen) + (" "*deletelen) + ("\b"*deletelen), end="", flush=True)

    print(("\b"*len(p_vec)) + (" "*len(p_vec)) + ("\b"*len(p_vec)), end="", flush=True)
print()
###### Stop commenting here to avoid rerunning simulation


###### Uncomment to avoid rerunning the simulation
# import joblib
# joblib.dump(pe_pls1, "pe_pls1.pkl")
# joblib.dump(pe_pls2, "pe_pls2.pkl")

# Uncomment only when the simulation is done running
# pe_pls1 = joblib.load("pe_pls1.pkl")
# pe_pls2 = joblib.load("pe_pls2.pkl")

# Figure 4.4
nc = 25
lim = [-7,7]

mu = 0
sigma = 1
x = np.linspace(mu - 3*sigma, mu + 3*sigma)

fig, axs = plt.subplots(3,3)

axs[0,0].hist(pe_pls1[0], color="green", bins=nc, range=lim, edgecolor="black", alpha=0.5, density=True)
axs[0,0].hist(pe_pls2[0, pe_pls2[0] < 6], color="blue", bins=nc, range=lim, edgecolor="black", alpha=0.5, density=True)
axs[0,0].title.set_text('n=4, p=8')
axs[0,0].plot(x, sp.stats.norm.pdf(x, mu, sigma), color="red")

axs[0,1].hist(pe_pls1[1], color="green", bins=nc, range=lim, edgecolor="black", alpha=0.5, density=True)
axs[0,1].hist(pe_pls2[1, pe_pls2[1] < 6], color="blue", bins=nc, range=lim, edgecolor="black", alpha=0.5, density=True)
axs[0,1].title.set_text('n=8, p=16')
axs[0,1].plot(x, sp.stats.norm.pdf(x, mu, sigma), color="red")

axs[0,2].hist(pe_pls1[2], color="green", bins=nc, range=lim, edgecolor="black", alpha=0.5, density=True)
axs[0,2].hist(pe_pls2[2, pe_pls2[2] < 6], color="blue", bins=nc, range=lim, edgecolor="black", alpha=0.5, density=True)
axs[0,2].title.set_text('n=16, p=32')
axs[0,2].plot(x, sp.stats.norm.pdf(x, mu, sigma), color="red")

axs[1,0].hist(pe_pls1[3], color="green", bins=nc, range=lim, edgecolor="black", alpha=0.5, density=True)
axs[1,0].hist(pe_pls2[3, pe_pls2[3] < 6], color="blue", bins=nc, range=lim, edgecolor="black", alpha=0.5, density=True)
axs[1,0].title.set_text('n=32, p=64')
axs[1,0].plot(x, sp.stats.norm.pdf(x, mu, sigma), color="red")

axs[1,1].hist(pe_pls1[4], color="green", bins=nc, range=lim, edgecolor="black", alpha=0.5, density=True)
axs[1,1].hist(pe_pls2[4, pe_pls2[4] < 6], color="blue", bins=nc, range=lim, edgecolor="black", alpha=0.5, density=True)
axs[1,1].title.set_text('n=64, p=128')
axs[1,1].plot(x, sp.stats.norm.pdf(x, mu, sigma), color="red")

axs[1,2].hist(pe_pls1[5], color="green", bins=nc, range=lim, edgecolor="black", alpha=0.5, density=True)
axs[1,2].hist(pe_pls2[5, pe_pls2[5] < 6], color="blue", bins=nc, range=lim, edgecolor="black", alpha=0.5, density=True)
axs[1,2].title.set_text('n=128, p=256')
axs[1,2].plot(x, sp.stats.norm.pdf(x, mu, sigma), color="red")

axs[2,0].hist(pe_pls1[6], color="green", bins=nc, range=lim, edgecolor="black", alpha=0.5, density=True)
axs[2,0].hist(pe_pls2[6, pe_pls2[6] < 6], color="blue", bins=nc, range=lim, edgecolor="black", alpha=0.5, density=True)
axs[2,0].title.set_text('n=256, p=512')
axs[2,0].plot(x, sp.stats.norm.pdf(x, mu, sigma), color="red")

axs[2,1].hist(pe_pls1[7], color="green", bins=nc, range=lim, edgecolor="black", alpha=0.5, density=True)
axs[2,1].hist(pe_pls2[7, pe_pls2[7] < 6], color="blue", bins=nc, range=lim, edgecolor="black", alpha=0.5, density=True)
axs[2,1].title.set_text('n=512, p=1024')
axs[2,1].plot(x, sp.stats.norm.pdf(x, mu, sigma), color="red")

axs[2,2].hist(pe_pls1[8], color="green", bins=nc, range=lim, edgecolor="black", alpha=0.5, density=True)
axs[2,2].hist(pe_pls2[8, pe_pls2[8] < 6], color="blue", bins=nc, range=lim, edgecolor="black", alpha=0.5, density=True)
axs[2,2].title.set_text('n=1024, p=2048')
axs[2,2].plot(x, sp.stats.norm.pdf(x, mu, sigma), color="red")

plt.show()

A = np.zeros((9,2))
B = np.zeros((9,2))

fig, axs = plt.subplots(1,2)
for j in range(9):
    A[j] = np.array([np.median(pe_pls1[j]), np.std(pe_pls1[j])])
    B[j] = np.array([np.median(pe_pls2[j]), np.std(pe_pls2[j,pe_pls2[j]<15])])

axs[0].plot(np.log(p_vec), A[:,0], color="red", linewidth=2)
axs[0].plot(np.log(p_vec), B[:,0], color="blue", linewidth=2)
axs[0].plot(np.log(p_vec), A[:,0], '.', marker="o", markerfacecolor="none", color="red", linewidth=2)
axs[0].plot(np.log(p_vec), B[:,0], '.', marker="o", markerfacecolor="none", color="blue", linewidth=2)
axs[0].axline(xy1=(0,0), slope=0, color="black")
axs[0].set_ylim(-3,10)
axs[0].set_xlim(1.5,8)
axs[0].set_ylabel("Mean")
axs[0].set_xlabel("log p")

axs[1].plot(np.log(p_vec), A[:,1], color="red", linewidth=2)
axs[1].plot(np.log(p_vec), B[:,1], color="blue", linewidth=2)
axs[1].plot(np.log(p_vec), A[:,1], '.', marker="o", markerfacecolor="none", color="red", linewidth=2)
axs[1].plot(np.log(p_vec), B[:,1], '.', marker="o", markerfacecolor="none", color="blue", linewidth=2)
axs[1].axline(xy1=(0,1), slope=0, color="black")
axs[1].set_ylim(0.5,4.25)
axs[1].set_xlim(1.5,8)
axs[1].set_ylabel("Standard deviation")
axs[1].set_xlabel("log p")

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

