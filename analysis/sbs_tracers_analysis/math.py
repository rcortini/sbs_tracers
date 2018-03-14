from scipy.stats import entropy

def KL_divergence(P, Q) :
    """
    Calculates the Kullback-Leibler divergence between probability distributions
    P and Q. The distributions can be non-normalized, but they must have
    positive sum, otherwise the function returns NaN.
    """
    if P.sum()>0 and Q.sum()>0 :
        return entropy(P, Q)
    else :
        return np.nan
