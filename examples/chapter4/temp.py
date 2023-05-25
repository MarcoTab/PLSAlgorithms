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
from IPython.display import display

np.random.seed(180)

ss = 10
JJ = np.array([int(i) for i in np.linspace(150, 1150, 11)])
nreps = 40
PLSS = np.zeros((nreps, len(JJ), 4))

qq = [0,0,0,0]

yeast = dict()

for i in range(len(JJ)):
    p = JJ[i]

    ## qq is how many predictors do not make contribution
    ##qq=0 all make contribution means Case 1 or d=p
    ##qq=40 all except 40 make contribution (almost Case 1
    ##qq=floor(p-p^(2/3)) means p^(2/3) makes contribution Case 2
    ##qq=p-2 only two make contribution 
    qq = [0, 40, p-int(np.floor(p**(2/3))), p-2]

    for jj in range(4):
        n = int(np.floor(p/3))
        q = qq[jj]

        for hh in range(nreps):
           
            if (q != p):
                H1 = np.random.multivariate_normal(mean=[0]*n, cov=25*np.eye(n), size=1)[0,:]
                H2 = np.random.multivariate_normal(mean=[0]*n, cov=25*np.eye(n), size=1)[0,:]
                H3 = None

                if (q != 0):
                    H3 = np.random.multivariate_normal(mean=[0]*n, cov=25*np.eye(n), size=1)[0,:]

                ns = [int(0), int(np.floor((p-q)/ 2)), int(p-q), int(p)]

                X = np.zeros((n, p))

                H1aux = np.stack([H1 for _ in range(ns[1])],axis=1)
                H2aux = np.stack([H2 for _ in range(ns[2]-ns[1])], axis=1)
                H3aux = None

                if (q!=0):
                    H3aux = np.stack([H3 for _ in range(ns[3]-ns[2])], axis=1)

                if (q!=0):
                    X = np.concatenate((
                        H1aux + np.random.normal(0,1,ns[1]*n).reshape(n, ns[1]),
                        H2aux + np.random.normal(0,1,(ns[2]-ns[1])*n).reshape(n, ns[2]-ns[1]),
                        H3aux + np.random.normal(0,1,(ns[3]-ns[2])*n).reshape(n, ns[3]-ns[2])
                    ), axis=1)


                if (q==0):
                    X = np.concatenate((
                        H1aux + np.random.normal(0,1,ns[1]*n).reshape(n, ns[1]),
                        H2aux + np.random.normal(0,1,(ns[2]-ns[1])*n).reshape(n, ns[2]-ns[1])
                    ), axis=1)

                means = 3*H1 - 4*H2

                y = np.array([np.random.normal(means[i], ss, 1) for i in range(len(means))])

                # print(y.shape)
                # exit()

                yeast['X'] = X
                yeast['y'] = y

            if (q == p):
                X = np.zeros((n, p))
                
                H3 = np.random.multivariate_normal(mean=[0]*n, cov=25*np.eye(n), size=1)[0,:]

                H3aux = np.stack([H3 for _ in range(q)], axis=1).reshape(n, q)

                X = H3aux + np.random.normal(0, 1, q*n).reshape(n,q)
                
                y = np.random.normal(0, ss, n).reshape(-1, 1)
                
                yeast['X'] = X
                yeast['y'] = y
            
            j = 1
            k = 2

            a1 = 75/(1+25*np.floor((p-q)/2))
            a2 = 100/(1+25*np.ceil((p-q)/2))

            beta_true = np.array([a1]*int(np.floor((p-q)/2)) + [-a2]*int(np.ceil((p-q)/2)) + [0]*q).reshape(-1,1)
            
            A1 = PLSRegression(n_components=2, scale=False)
            A1.fit(yeast['X'], yeast['y'])

            PLSS[hh,i,jj] = ((beta_true - A1.coef_).T @ np.cov(yeast['X'], ddof=1, rowvar=False) @ (beta_true - A1.coef_)).item()


H = np.zeros((len(JJ), len(qq)))

for i in range(nreps):
    H += PLSS[i,:,:]

H /= nreps

## FIGURE 4.7

plt.plot(JJ, H[:,0], color="blue", linewidth=2.5)
# plt.scatter(JJ, H[:,0], facecolor="none")
plt.plot(JJ, H[:,1], color="red", linewidth=2.5)

plt.plot(JJ, H[:,2], color="orange", linewidth=2.5)

plt.plot(JJ, H[:,3], color="green", linewidth=2.5)

plt.ylim(0, 25)
plt.xlim(100, 1200)
plt.title('Figure 4.7')
plt.rcParams['text.usetex'] = True
plt.xlabel('Number of predictors, $p$')
plt.ylabel(r"$||\beta - \hat{\beta}||^2_{S_X}$")

plt.show()

                
