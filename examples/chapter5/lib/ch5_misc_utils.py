import numpy as np

### Function to calculate Q of A under inner product matrix M
def cal_Q_fun(A, M):
    P_A = A @ np.linalg.inv(A.T @ M @ A) @ A.T @ M

    k = A.shape[0]

    return np.eye(k) - P_A


### Function to calculate two blocks PLS in a single step
def twoblock_pls_onestep_fun(mymatrix):
    try:
        svdu, svdsv, svdv = np.linalg.svd(mymatrix, full_matrices=False)
        svdu = svdu[:, 0]
        svdv = svdv[:, 0]
    except Exception as e:
        svdu, svdsv, svdv = np.linalg.svd(mymatrix.T, full_matrices=False)
        svdu = svdu[:, 0]
        svdv = svdv[:, 0]

    return {
        "u": svdu,
        "v": svdv
    }


### Function to calculate two blocks PLS
def twoblock_pls_bymatrix_fun(d, SigmaX, SigmaY, SigmaXY):
    mymatrix = SigmaXY
    u_m = np.zeros((SigmaX.shape[1],d))
    v_m = np.zeros((SigmaY.shape[1],d))

    for i in range(d):
        twoblock_pls_onestep_result = twoblock_pls_onestep_fun(mymatrix)
        u_current = twoblock_pls_onestep_result["u"]
        v_current = twoblock_pls_onestep_result["v"]

        u_m[:, i] = u_current
        v_m[:, i] = v_current

        mymatrix = cal_Q_fun(u_m[:,:i], SigmaX).T @ SigmaXY @ cal_Q_fun(v_m[:, :i], SigmaY)

    _, svs, _ = np.linalg.svd(mymatrix, full_matrices=False)

    return {
        "u_m": u_m,
        "v_m": v_m,
        "singular_value_v": svs
    }


### Function to calculate the dimension of two block PLS in population
def cal_d_twoblock_fun(SigmaX, SigmaY, SigmaXY, mytol=1e-3):
    maxdim = min(SigmaX.shape[0], SigmaY.shape[0])

    d = 0

    _, svs, _ = np.linalg.svd(SigmaXY, full_matrices=False)
    singular_value_sum = svs.sum()

    while d <= maxdim and singular_value_sum > mytol:
        d += 1

        twoblock_result = twoblock_pls_bymatrix_fun(d,SigmaX=SigmaX, SigmaY=SigmaY, SigmaXY=SigmaXY)

        _, svs, _ = np.linalg.svd(twoblock_result["singular_value_v"], full_matrices=False)
        singular_value_sum = svs.sum()
    
    return d


