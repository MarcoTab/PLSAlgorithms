import numpy as np

def is_psd(A):
    return np.all(np.linalg.eigh(A)[0] > -1e-10)

# Given a positive definite matrix and a vector, finds the orthogonal projection wrt the dot product generated with that matrix
def orth_proj(M, v):
    # if not is_psd(M):
    #     raise ValueError("M must be positive definite")
    
    m = M.shape[0]
    return np.eye(m) - v @ np.linalg.inv(v.T @ M @ v) @ v.T @ M

# Given positive semi-definite matrices A and S, calculates the d first PLS projections.
# A is Sigma_xy sigma_y^-1 Sigma_xy
# S is Sigma_x
# Returns Gamma, Omega, beta.
def pls_matrices(A, S, d):
    # if not is_psd(A):
    #     raise ValueError("A must be positive semi-definite")
    # if not is_psd(S):
    #     raise ValueError("S must be positive semi-definite")
    
    V = np.zeros((A.shape[0], d))

    _, evecs = np.linalg.eigh(A)
    evecs = evecs[:,::-1]

    V[:,0] = evecs[:,0]

    if (d > 1):
        for i in range(1, d):
            aux = orth_proj(S, V[:, 0:i])
            _, evecs = np.linalg.eigh(aux.T @ A @ aux)
            evecs = evecs[:,::-1]

            V[:, i] = evecs[:, 0]
    

    return {
        'Gamma':V,
        'Omega':(V.T @ A @ V),
        'Atilde': V @ V.T @ A @ V @ V.T,
        'beta': V @ np.linalg.inv(V.T @ S @ V) @ V.T @ A
    }