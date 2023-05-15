# Chapter 1

import warnings
warnings.filterwarnings("ignore")

import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_val_score
import matplotlib.pyplot as plt

# Read in corn data
X_train = np.genfromtxt('mvc1_data/Parameters in corn/Xcal.txt', delimiter="  ")
X_train = X_train.T
Y_train = np.genfromtxt('mvc1_data/Parameters in corn/y1cal.txt', delimiter="  ")

X_test = np.genfromtxt('mvc1_data/Parameters in corn/Xtest.txt', delimiter="  ")
X_test = X_test.T
Y_test = np.genfromtxt('mvc1_data/Parameters in corn/y1test.txt', delimiter="  ")

scores = []
stds = []
print("Fitting for moisture in corn", end="")
for i in range(1, 51):
    print("." if (i) % 4 != 0 else "\b\b\b   \b\b\b", end="", flush=True)
    # Instantiate a PLSRegression object. This uses NIPALS as the algorithm for computing PLS.
    # n_components decides how many components to keep, scale decides whether to scale the data by the s.d., and copy decides whether to center and scale the data in-place or not.
    nipals_pls = PLSRegression(n_components=i, scale=False, copy=True)

    # Use 5 fold cross validation to score this instance of the model, score with RMSE.
    this_score = cross_val_score(nipals_pls, X_train, Y_train, scoring='neg_mean_squared_error', cv=10)
    scores.append(np.mean(this_score))
    stds.append(np.std(this_score))

print()

# Get the best performing n from the bunch. Remember that cross_val_score _SCORES_ the model, it doesn't measure error in the traditional sense.
best_n = np.argmax(scores)+1

# Again, instantiate the object first.
best_nipals = PLSRegression(n_components=best_n, scale=False, copy=True)
# Fit the model to the training data
best_nipals.fit(X_train, Y_train)
# Predict on the testing data
preds = best_nipals.predict(X_test)

# Compute RMSE
rmse = mean_squared_error(Y_test, preds, squared=False)

print(f"Moisture in corn: Best n = {best_n}, with RMSE of {rmse}")

plt.scatter(Y_test, preds, facecolors='none', edgecolors='black')
# plt.title("NIPALS Observed vs Predicted Response")
plt.xlabel("Observed response, Y")
plt.ylabel("Predicted response")
plt.show()


# Read in meat data
X_train = np.genfromtxt('mvc1_data/Parameters in meat/Xcal.txt', delimiter="  ")
X_train = X_train.T
Y_train = np.genfromtxt('mvc1_data/Parameters in meat/y3cal.txt', delimiter="  ")

X_test = np.genfromtxt('mvc1_data/Parameters in meat/Xtest.txt', delimiter="  ")
X_test = X_test.T
Y_test = np.genfromtxt('mvc1_data/Parameters in meat/y3test.txt', delimiter="  ")

scores = []
stds = []

print("Fitting for protein in meat", end="")
for i in range(1, 51):
    print("." if i % 4 != 0 else "\b\b\b   \b\b\b", end="", flush=True)
    # Instantiate a PLSRegression object. This uses NIPALS as the algorithm for computing PLS.
    # n_components decides how many components to keep, scale decides whether to scale the data by the s.d., and copy decides whether to center and scale the data in-place or not.
    nipals_pls = PLSRegression(n_components=i, scale=False, copy=True)

    # Use 5 fold cross validation to score this instance of the model, score with RMSE.
    this_score = cross_val_score(nipals_pls, X_train, Y_train, scoring='neg_mean_squared_error', cv=10)
    scores.append(np.mean(this_score))
    stds.append(np.std(this_score))

    # print(f"Scored {scores[-1]:.10f} with n_components = {i}")
print()

# Get the best performing n from the bunch. Remember that cross_val_score _SCORES_ the model, it doesn't measure error in the traditional sense.
best_n = np.argmax(scores)+1

# Again, instantiate the object first.
best_nipals = PLSRegression(n_components=best_n, scale=False, copy=True)
# Fit the model to the training data
best_nipals.fit(X_train, Y_train)
# Predict on the testing data
preds = best_nipals.predict(X_test)

# Compute RMSE
rmse = mean_squared_error(Y_test, preds, squared=False)

print(f"Protein in meat: Best n = {best_n}, with RMSE of {rmse}")

plt.scatter(Y_test, preds, facecolors='none', edgecolors='black')
# plt.title("NIPALS Observed vs Predicted Response")
plt.xlabel("Observed response, Y")
plt.ylabel("Predicted response")
plt.show()


# Read in Tetracycline data
X_train = np.genfromtxt('mvc1_data/Tetracycline in serum/Xcal.txt', delimiter="  ")
X_train = X_train.T
Y_train = np.genfromtxt('mvc1_data/Tetracycline in serum/ycal.txt', delimiter="  ")

X_test = np.genfromtxt('mvc1_data/Tetracycline in serum/Xtest.txt', delimiter="  ")
X_test = X_test.T
Y_test = np.genfromtxt('mvc1_data/Tetracycline in serum/ytest.txt', delimiter="  ")

scores = []
stds = []

print("Fitting for tetracycline in serum", end="")
for i in range(1, 51):
    print("." if i % 4 != 0 else "\b\b\b   \b\b\b", end="", flush=True)
    # Instantiate a PLSRegression object. This uses NIPALS as the algorithm for computing PLS.
    # n_components decides how many components to keep, scale decides whether to scale the data by the s.d., and copy decides whether to center and scale the data in-place or not.
    nipals_pls = PLSRegression(n_components=i, scale=False, copy=True)

    # Use 5 fold cross validation to score this instance of the model, score with RMSE.
    this_score = cross_val_score(nipals_pls, X_train, Y_train, scoring='neg_mean_squared_error', cv=10)
    scores.append(np.mean(this_score))
    stds.append(np.std(this_score))

    # print(f"Scored {scores[-1]:.10f} with n_components = {i}")
print()

# Get the best performing n from the bunch. Remember that cross_val_score _SCORES_ the model, it doesn't measure error in the traditional sense.
best_n = np.argmax(scores)+1

# Again, instantiate the object first.
best_nipals = PLSRegression(n_components=best_n, scale=False, copy=True)
# Fit the model to the training data
best_nipals.fit(X_train, Y_train)
# Predict on the testing data
preds = best_nipals.predict(X_test)

# Compute RMSE
rmse = mean_squared_error(Y_test, preds, squared=False)

print(f"Tetracycline in serum: Best n = {best_n}, with RMSE of {rmse}")

plt.scatter(Y_test, preds, facecolors='none', edgecolors='black')
# plt.title("NIPALS Observed vs Predicted Response")
plt.xlabel("Observed response, Y")
plt.ylabel("Predicted response")
plt.show()

