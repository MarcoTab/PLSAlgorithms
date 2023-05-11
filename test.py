from PLSLib.nipals import nipals
import numpy as np

X_train = np.genfromtxt('mvc1_data/Parameters in corn/Xcal.txt', delimiter="  ")
X_train = X_train.T
Y_train = np.atleast_2d(np.genfromtxt('mvc1_data/Parameters in corn/y1cal.txt', delimiter="  "))
Y_train = Y_train.T

X_test = np.genfromtxt('mvc1_data/Parameters in corn/Xtest.txt', delimiter="  ")
X_test = X_test.T
Y_test = np.atleast_2d(np.genfromtxt('mvc1_data/Parameters in corn/y1test.txt', delimiter="  "))
Y_test = Y_test.T

xbar = np.atleast_2d(np.mean(X_train, axis=0))
ybar = np.atleast_2d(np.mean(Y_train, axis=0))

# print(f"{X_train.shape=} {Y_train.shape=} {xbar.shape=} {ybar.shape=} {X_test.shape=} {Y_test.shape=}")

my_nipals = nipals()

my_nipals.fit(X_train, Y_train, 9)
W, beta = my_nipals.get_transformers()

RMSE = (Y_test - (beta.T @ (X_test - xbar).T + ybar)) ** 2

print(RMSE)