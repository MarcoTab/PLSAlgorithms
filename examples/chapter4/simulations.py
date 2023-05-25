# Chapter 4
# Figure 4.4 and Table 4.4
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