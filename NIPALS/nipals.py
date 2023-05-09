class nipals:
    """
        Orthogonal weights: :math:`W^T_q W_q = I_q`

        Envelope connection: :math:`\mathrm{span}(W_q) = \mathcal{E}_{\Sigma_X}(\mathcal{B})`, the :math:`\Sigma_X`-envelope of :math:`\mathcal{B} \mathrel{\\vcenter{:}}= \mathrm{span}(\\beta)`.

        Score matrix :math:`S_d`: These are traditional computational intermediaries,
        although they are not needed in the computation of :math:`\hat{\\beta}_{\mathrm{npls}}`.

        Algorithm :math:`\mathbb{N}`: This is an instance of Algorithm :math:`\mathbb{N}` discussed in §1.5.3.

        PLS1 v. PLS2: Algorithm is applicable for PLS1 or PLS2 fits; See §3.8.
    """
    def __init__(self):
        self.q = None
        self.W = None
        self.beta = None

    def fit(self, X, Y, q, version='sample'):
        """    
            Fit this model to the training data `X`, `Y` using `q` dimensions. 

            :param X: Predictor of shape (`n_samples`, `p_features`)
            :type X: array-like
            :param Y: Response of shape (`n_samples`, `r_features`)
            :type Y: array-like
            :param q: Value between `1` and `p_features`. The number of projections used.
            :type q: int
            :param version: either 'sample' or 'population', defaults to 'sample'
            :type version: str
            :return: Nothing.
        """ 
        if (version != 'sample' or version != 'population'):
            raise ValueError(f"version must be 'sample' or 'population', not {version}.")
        
        if (len(X.shape) != 2 or len(Y.shape) != 2):
            raise ValueError(f"X and Y must be 2 dimensional array-like. Current dimensions: {len(X.shape)} and {len(Y.shape)}.")
        
        if (X.shape[0] != Y.shape[0]):
            raise ValueError(f"X and Y must have the same first dimension. Current shapes: {X.shape} and {Y.shape}.")
        
        q = int(q)

        if (q < 1 or q > X.shape[1]):
            raise ValueError(f"q should be 1 and X's second dimension. Currently {q} and {X.shape[1]}.")
        
        
        ...

    def transform(self, X):
        """
            Transform data using the NIPALS algorithm. Must run :meth:`NIPALS.nipals.nipals.fit` before running this function.

            :param X: Predictor of shape (`n_samples`, `p_features`)
            :type X: array-like
            :return: The :math:`W` and :math:`\\beta` transformed data, respectively.
            :rtype: tuple(array-like, array-like)
        """

        if (self.q is None or self.W is None or self.beta is None):
            raise Exception("You must run `fit` before running this function.")
        
        return self.W.T @ X, self.beta.T @ X
        ...

