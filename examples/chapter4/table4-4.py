# Chapter 4
# Table 4.4
#### READ BEFORE RUNNING ALL THE WAY THROUGH
# There is a very long simulation that runs as part of this script. 
# To avoid having to run this simulation again and again, please uncomment the code below (you might have to scroll a
# long ways, or search "uncomment") to save the data to your working directory.
# You must have the joblib package installed to do so.

# Some numbers will be different between the book and this script, since the random number generation for Python and R are different.

import warnings
warnings.filterwarnings("ignore")
from IPython.display import display

import numpy as np
import pandas as pd
import scipy as sp
from sklearn.cross_decomposition import PLSRegression


###### Simulation to get Table 4.4

p_vec = np.array([2**i for i in range(8,11)])
pmax = p_vec[-1]
reps = 1000
n = int(np.floor(np.sqrt(pmax)))
pow = 3/4
s = pow
ovec = np.array([1] + [0]*(pmax-1)).reshape(-1, 1)
p = p_vec[0]    

if (pow > 0):
    ovec = np.array([1]*(round(p_vec[0]**pow)) + [0]*(p_vec[0] - round(p_vec[0]**pow))).reshape(-1, 1)
    for pk in p_vec[1:]:
        b = len(ovec)
        p = pk
        cnz = sum(ovec)
        tnz = round(p**pow)
        nnz = (tnz-cnz).item()
        newvec = np.array([1]*nnz + [0]*(b-nnz)).reshape(-1, 1)
        ovec = np.concatenate((ovec, newvec))
        
sigmaxy = ovec*(np.random.normal(0, 1, p).reshape(-1, 1))
sigma_0 = np.linalg.qr(sigmaxy, mode="complete")[0][:,1:]

sTs = (sigmaxy.T @ sigmaxy).item()
Omega = sTs
tau = 1
sigmaY2 = tau**2 + sTs*(Omega**(-1))
Omega0 = np.eye(p-1)
Sigma = (sigmaxy @ sigmaxy.T / sTs) * Omega + sigma_0 @ Omega0 @ sigma_0.T
out_eigenvals, out_eigenvecs = np.linalg.eigh(Sigma)

out_eigenvals = out_eigenvals[::-1]
out_eigenvecs = out_eigenvecs[:,::-1]

Sigma_sqrt = out_eigenvecs @ np.diag(np.sqrt(out_eigenvals)) @ out_eigenvecs.T

Xbig = np.empty((n, p, reps))
Ybig = np.empty((n, reps))
true_beta = (Omega**(-1))*sigmaxy

for j in range(reps):
    X = np.random.normal(0, 1, (n,p)) @ Sigma_sqrt
    y = np.array([np.random.normal((X @ true_beta)[i], tau, 1) for i in range(n)])
    Xbig[:,:,j] = X
    Ybig[:,j] = y.flatten()

# Prep
nmax = Xbig.shape[0]
pmax = Xbig.shape[1]
k = Xbig.shape[2]
nvec = np.floor(np.sqrt(p_vec))
XN = np.random.multivariate_normal([0]*pmax, cov=np.eye(pmax))
cover1 = np.empty((len(p_vec), k))
cover2 = np.empty((len(p_vec), k))
alfa = 0.05
z = sp.stats.norm.ppf(1-alfa/2, loc=0, scale=1)

# Sim1
print("First Sim (Table 4.4)", end="", flush=True)
for i in range(len(p_vec)):
    n = int(nvec[i])
    p = p_vec[i]

    m = n-1
    G = (XN[:p]).reshape(-1, 1)
    Omega = (sigmaxy[:p].T @ sigmaxy[:p]).item()
    beta = Omega**(-1) * sigmaxy[:p]
    true_val = (beta.T @ G).item()
    Cp1 = p-1
    Cp2 = p-1
    sts = (sigmaxy[:p].T @ sigmaxy[:p]).item()
    bias = -((sts/Omega - tau**2)*(Cp1/(m*sts)) + sigmaY2*(Cp1**2)/((m**2)*Omega*sts) + sigmaY2*Cp2/(m*Omega*sts))
    div = k//5
    for j in range(k):
        print("." if j%div == 0 else "", end="", flush=True)
        X = Xbig[:n, :p, j]
        y = (Ybig[:n, j]).reshape(-1, 1)
        mX = np.mean(X, axis=1).reshape(-1,1)
        X = (X - mX)
        y = y-np.mean(y)

        mod = PLSRegression(n_components=1, scale=False)
        mod.fit(X, y)
        betahat = mod.coef_
        pred = (betahat.T @ G).item()
        shat = np.cov(np.concatenate((y, X), axis=1), rowvar=False)[0, 1:].reshape(-1,1)
        SIhat = np.cov(X, rowvar=False)
        stshat = (shat.T @ shat).item()

        sigmaY2hat = np.var(y, ddof=1)
        sigma_0hat = np.linalg.qr(shat, mode="complete")[0][:,1:p]
        Omegahat = (shat.T @ SIhat @ shat / stshat).item()
        Omega0hat = sigma_0hat.T @ SIhat @ sigma_0hat

        SIhat2 = (Omegahat/stshat)*(shat @ shat.T) + sigma_0hat @ Omega0hat @ sigma_0hat.T

        Vhat = (((Omegahat**(-2))/n) * G.T @ (SIhat2*sigmaY2hat - shat@shat.T)@G).item()
        
        IC1 = pred + np.array([-1, 1])*z*np.sqrt(Vhat)

        cover1[i,j] = IC1[0] <= true_val and true_val <= IC1[1]
        cover2[i,j] = IC1[0] <= true_val*(1+bias) and true_val*(1+bias) <= IC1[1]

    print("\b"*(5)+" "*(5)+"\b"*(5), end="", flush=True)
print()

############ Second Sim
p_vec = np.array([2**i for i in range(8,11)])
pmax = p_vec[-1]
reps = 1000
n = int(np.floor(np.sqrt(pmax)))
pow = 3/4
s = pow
ovec = np.array([1] + [0]*(pmax-1)).reshape(-1, 1)
p = p_vec[0]    

if (pow > 0):
    ovec = np.array([1]*(round(p_vec[0]**pow)) + [0]*(p_vec[0] - round(p_vec[0]**pow))).reshape(-1, 1)
    for pk in p_vec[1:]:
        b = len(ovec)
        p = pk
        cnz = sum(ovec)
        tnz = round(p**pow)
        nnz = (tnz-cnz).item()
        newvec = np.array([1]*nnz + [0]*(b-nnz)).reshape(-1, 1)
        ovec = np.concatenate((ovec, newvec))
        
sigmaxy = ovec*(np.random.normal(0, 1, p).reshape(-1, 1))

sigma_0 = np.linalg.qr(sigmaxy, mode="complete")[0][:,1:]

sTs = (sigmaxy.T @ sigmaxy).item()
Omega = sTs
tau = 1/2
sigmaY2 = tau**2 + sTs*(Omega**(-1))
Omega0 = np.eye(p-1)
Sigma = (sigmaxy @ sigmaxy.T / sTs) * Omega + sigma_0 @ Omega0 @ sigma_0.T
out_eigenvals, out_eigenvecs = np.linalg.eigh(Sigma)

out_eigenvals = out_eigenvals[::-1]
out_eigenvecs = out_eigenvecs[:,::-1]

Sigma_sqrt = out_eigenvecs @ np.diag(np.sqrt(out_eigenvals)) @ out_eigenvecs.T
Xbig = np.empty((n, p, reps))
Ybig = np.empty((n, reps))
true_beta = (Omega**(-1))*sigmaxy

for j in range(reps):
    X = np.random.normal(0, 1, (n,p)) @ Sigma_sqrt
    y = np.array([np.random.normal((X @ true_beta)[i], tau, 1) for i in range(n)])
    Xbig[:,:,j] = X
    Ybig[:,j] = y.flatten()

# Prep
nmax = Xbig.shape[0]
pmax = Xbig.shape[1]
k = Xbig.shape[2]
nvec = np.floor(np.sqrt(p_vec))
cover3 = np.empty((len(p_vec), k))
cover4 = np.empty((len(p_vec), k))
cover5 = np.empty((len(p_vec), k))
alfa = 0.05
z = sp.stats.norm.ppf(1-alfa/2, loc=0, scale=1)

# Sim2
print("Second Sim (Table 4.4)", end="", flush=True)
for i in range(len(p_vec)):
    n = int(nvec[i])
    p = p_vec[i]
    # print(f"{n=}, {p=}")

    m = n-1
    Omega = (sigmaxy[:p].T @ sigmaxy[:p]).item()
    G = (sigmaxy[:p]).reshape(-1, 1)
    beta = Omega**(-1) * sigmaxy[:p]
    true_val = (beta.T @ G).item()
    Cp1 = p-1
    Cp2 = p-1
    sts = (sigmaxy[:p].T @ sigmaxy[:p]).item()
    bias = -((sts/Omega - tau**2)*(Cp1/(m*sts)) + sigmaY2*(Cp1**2)/((m**2)*Omega*sts) + sigmaY2*Cp2/(m*Omega*sts))
    div = k//5
    for j in range(k):
        print("." if j%div == 0 else "", end="", flush=True)
        X = Xbig[:n, :p, j]
        y = (Ybig[:n, j]).reshape(-1, 1)
        mX = np.mean(X, axis=1).reshape(-1,1)
        X = (X - mX)
        y = y-np.mean(y)

        mod = PLSRegression(n_components=1, scale=False)
        mod.fit(X, y)
        betahat = mod.coef_
        pred = (betahat.T @ G).item()
        shat = np.cov(np.concatenate((y, X), axis=1), rowvar=False)[0, 1:].reshape(-1,1)
        SIhat = np.cov(X, rowvar=False)
        stshat = (shat.T @ shat).item()

        sigmaY2hat = np.var(y, ddof=1)
        sigma_0hat = np.linalg.qr(shat, mode="complete")[0][:,1:p]
        Omegahat = (shat.T @ SIhat @ shat / stshat).item()
        Omega0hat = sigma_0hat.T @ SIhat @ sigma_0hat
        Cp1hat = np.trace(Omega0hat)
        Cp2hat = np.trace(Omega0hat @ Omega0hat)
        tau2hat = np.var(y - (mod.predict(X)), ddof=1)

        biashat = -((stshat / Omegahat - tau2hat)*(Cp1hat/(m*stshat)) + sigmaY2hat*(Cp1hat**2)/((m**2)*Omegahat*stshat) + sigmaY2hat*Cp2hat/(m*Omegahat*stshat))

        V = (((Omega**(-2))/n) * G.T @ (Sigma[:p, :p] * sigmaY2 - sigmaxy[:p]@sigmaxy[:p].T) @ G).item()
        IC1 = pred + np.array([-1, 1])*z*np.sqrt(V)
        IC2 = IC1/(1+biashat)

        cover3[i,j] = IC1[0] <= true_val and true_val <= IC1[1]
        cover4[i,j] = IC1[0] <= true_val*(1+bias) and true_val*(1+bias) <= IC1[1]
        cover5[i,j] = IC2[0] <= true_val and true_val <= IC2[1]

    # print(f"5th:{np.mean(cover3, axis=1)[i]}")
    # print(f"6th:{np.mean(cover4, axis=1)[i]}")
    # print(f"7th:{np.mean(cover5, axis=1)[i]}")
    # print()

    print("\b"*(5)+" "*(5)+"\b"*(5), end="", flush=True)
print()


colnames = pd.DataFrame([
    ["", 'n'],
    ["", 'p'],
    ["√n|b| → 0", 'βᵀG'],
    ["√n|b| → 0", '(1+b)βᵀG'],
    ['√n|b| ↛ 0', 'βᵀG'],
    ['√n|b| ↛ 0', '(1+b)βᵀG'],
    ['', 'βᵀG']
])

rows = np.array([nvec, p_vec, np.mean(cover1, axis=1), np.mean(cover2, axis=1), np.mean(cover3, axis=1), np.mean(cover4, axis=1), np.mean(cover5, axis=1)]).T

cols = pd.MultiIndex.from_frame(colnames)
df = pd.DataFrame(rows, columns=cols)
df.style.set_properties(**{'text-align': 'center'})
pd.set_option('display.max_columns', None)
display(df)