import numpy as np
from scipy import optimize

def __T1Bp_help(tao, S0, T1):
    return S0 * (1 - np.exp(-tao/T1))
## For a function S(tao) = S0 (1 - exp(-tao/T1)). Start by using non-linear least square, then iterate.
def T1Bp(tao, S):
    T1_prev = 0
    popt, var = optimize.curve_fit(__T1Bp_help, tao, S)
    S0 = popt[0]
    T1 = popt[1]
    print('Initial guess (S0 first): ' + str([S0, T1]))
    print('Initial std (S0 first): ' + str(np.sqrt(np.diag(var))))

    ### Can be done, but I cannot figure out the variance matrix and how it relates to shit
    # tol = 1e-4
    # print("Initial uncertainty: " + str(np.sqrt(np.diag(var))))
    # while np.abs(T1-T1_prev) > tol:
    #     diff = S0 - S
    #     z = np.log(diff)
    #     [p, var] = np.polyfit(tao,z,1, full=False, cov=True)
    #     #print("...uncertainty: " + str((np.sqrt(np.diag(var)))))
    #     S0 = np.exp(p[1])
    #     T1_prev = T1
    #     T1 = -1/p[0]
    # print('Polyfit guess (T1 first): ' + str([T1, S0]))
    # print('Polyfit std (T1 first): ' + str(np.sqrt(np.diag(var))))
    # popt, var = optimize.curve_fit(__T1Bp_help, tao, S, p0 = [S0, T1])
    # S0 = popt[0]
    # T1 = popt[1]
    # print('Final guess (S0 first): ' + str([S0, T1]))
    # print('Final std (S0 first): ' + str(np.sqrt(np.diag(var))))

    curve_of_best_fit = __T1Bp_help(tao, S0, T1)
    return [curve_of_best_fit, np.array([S0, T1]), var]


def __T1BE_help(tao, S0, T1):
    return S0 * np.exp(-tao/T1)
## For a function S(tao) = S0 exp(-tao/T1)
def T1BE(tao, S):
    logS = np.log(S)
    [p, var] = np.polyfit(tao, logS, 1, full=False, cov=True)
    S0 = np.exp(p[1])
    T1 = -1/p[0]
    print('Polyfitted coefficients (T1 first): ' + str([T1,S0]))
    print('Polyfitted variance (T1 first): ' + str(np.diag(var)))
    print('Polyfitted std (T1 first): ' + str(np.sqrt(np.diag(var))))

    ### I do this to get a variance matrix I understand
    popt, var = optimize.curve_fit(__T1BE_help, tao, S, p0 = [S0, T1])
    S0 = popt[0]
    T1 = popt[1]

    curve_of_best_fit = __T1BE_help(tao, S0, T1)
    return [curve_of_best_fit, np.array([S0, T1]), var]

def __T2_help(two_tao, S0, T2):
    return S0 * np.exp(-two_tao/T2)
## For a function S(tao) = S0 exp(-2tao/T2). Please note 2tao is a variable, not 2 times tao
def T2__(two_tao, S):
    logS = np.log(S)
    [p, var] = np.polyfit(two_tao, logS, 1, full=False, cov=True)
    S0 = np.exp(p[1])
    T2 = -1/p[0]
    print('Polyfitted coefficients (T2 first): ' + str([T2, S0]))
    print('Polyfitted variance (T2 first): ' + str(np.diag(var)))
    print('Polyfitted std (T2 first): ' + str(np.sqrt(np.diag(var))))

    ### I do this to get a variance matrix I understand
    popt, var = optimize.curve_fit(__T2_help, two_tao, S, p0 = [S0, T2])
    S0 = popt[0]
    T2 = popt[1]

    curve_of_best_fit = __T2_help(two_tao, S0, T2)
    return [curve_of_best_fit, np.array([S0, T2]), var]