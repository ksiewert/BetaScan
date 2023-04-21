import sys
import numpy as np
import argparse
import math


def find_win_indx(prev_start_i, prev_end_i, snp_i, data_list, win_size):
    """Takes in the previous indices of the start_ing and end of the window,
    then returns the appropriate start_ing and ending index for the next SNP

    Parameters:
        prev_start_i: start_ing index in the array of SNP for the previous core SNP's window, inclusive
        prev_end_i: ending index in the array for the previous SNP's window, inclusive
        snp_i, the index in the array for the current SNP under consideration
        data_list: the numpy array of all SNP locations & frequencies
    """

    loc_snp = data_list[snp_i, 0]  # the coordinates of the core SNP
    win_start = loc_snp - win_size/2

    # array index of start of window, inclusive
    firstI = prev_start_i + np.searchsorted(data_list[prev_start_i:, 0], win_start, side='left')
    winEnd = loc_snp + win_size/2

    # array index of end of window, exclusive
    endI = prev_end_i - 1 + np.searchsorted(data_list[prev_end_i:, 0], winEnd, side='right')
    return (firstI, endI)


def calc_beta_folded(snp_freq_list, core_freq, num_ind, p):
    """Calculates the value of the folded beta statistic

    Parameters:
        SNPFreq: freq of SNP under consideration, ranges from 1 to sample size
        core_freq: freq of coresite, ranges from 0 to 1
        p: the p parameter specificying sharpness of peak
        num_ind: the number of haploid individuals used to calculate frequency of core site
    """

    if snp_freq_list.size == 0:
        return 0
    a1 = np.sum(1./np.arange(1, num_ind))
    thetaW = len(snp_freq_list[:, 0])/a1
    thetaBNum = np.sum(calc_d(snp_freq_list[:, 0]/snp_freq_list[:, 1], core_freq, p))

    i = np.arange(1, num_ind)
    thetaBDenom = np.sum((1./i)*calc_d(i/float(num_ind), core_freq, p))

    thetaB = thetaBNum/thetaBDenom
    return thetaB - thetaW


def calc_beta_unfolded(snp_freq_list, core_freq, num_ind, p):
    """Calculates the unfolded version of Beta from Siewert and Voight
        For use when the ancestral and derived alleles can be confidently called

    Parameters:
        snp_freq_list: a list of frequencies, one for each SNP in the window,
            first column ranges from 1 to number of individuals, second columns is # individuals
        core_freq: the frequency of the core SNP, must range from 0 to 1, exclusive
        num_ind: number of individuals used to calculate the core site frequency
        p: value of parameter p
    """
    if snp_freq_list.size == 0:
        return 0
    a1 = np.sum(1./np.arange(1, num_ind))
    thetaW = len(snp_freq_list[:, 0])/a1
    thetaBNum = sum(calc_d(snp_freq_list[:, 0]/snp_freq_list[:, 1], core_freq, p) * snp_freq_list[:, 0])
    thetaBDenom = np.sum(calc_d(np.arange(1, num_ind)/float(num_ind), core_freq, p))
    thetaB = thetaBNum/thetaBDenom
    return thetaB - thetaW


def calc_thetabeta_unfolded(snp_freq_list, core_freq, num_ind, p):
    """Calculates theta_Beta usign the unfolded SFS

    Parameters:
        SNPFreq: freq of SNP under consideration, ranges from 1 to sample size
        core_freq: freq of coresite, ranges from 0 to 1
        p: the p parameter specificying sharpness of peak
        num_ind: the number of haploid individuals used to calculate frequency of core site
    """

    if snp_freq_list.size == 0:
        return 0
    thetaBNum = np.sum(calc_d(snp_freq_list[:, 0]/snp_freq_list[:, 1], core_freq, p) * snp_freq_list[:, 0])

    thetaBDenom = np.sum(calc_d(np.arange(1, num_ind)/float(num_ind), core_freq, p))

    thetaB = thetaBNum/thetaBDenom
    return thetaB


def calc_thetabeta_folded(snp_freq_list, core_freq, num_ind, p):
    """Calculates theta_Beta using the folded SFS

    Parameters:
        SNPFreq: freq of SNP under consideration, ranges from 1 to sample size
        core_freq: freq of coresite, ranges from 0 to 1
        p: the p parameter specificying sharpness of peak
        num_ind: the number of haploid individuals used to calculate frequency of core site
    """

    if snp_freq_list.size == 0:
        return 0
    thetaBNum = np.sum(calc_d(snp_freq_list[:, 0]/snp_freq_list[:, 1], core_freq, p))

    thetaBDenom = np.sum((1./np.arange(1, num_ind))*calc_d(np.arange(1, num_ind) / float(num_ind), core_freq, p))

    thetaB = thetaBNum/thetaBDenom
    return thetaB


def calc_thetaw_unfolded(snp_freq_list, num_ind):
    """Calculates watterson's theta

    Parameters:
        snp_freq_list: a list of frequencies, one for each SNP in the window,
            first column ranges from 1 to number of individuals, second columns is # individuals
        num_ind: number of individuals used to calculate the core site frequency
    """
    if snp_freq_list.size == 0:
        return 0

    a1 = np.sum(1./np.arange(1, num_ind))

    thetaW = len(snp_freq_list[:, 0])/a1
    return thetaW


def calc_theta_d(snp_freq_list, c, n):
    """
    Calculates theta_D

    Parameters:
        c: Speciation time in coalescent units
        n: Sample Size
    """
    if snp_freq_list.size == 0:
        return 0

    S = np.where(snp_freq_list[:, 0] == snp_freq_list[:, 1])[0].shape[0]
    return S/(c+1./n)


def calc_beta_2(snp_freq_list, c, n, core_freq, p):
    SNPs = snp_freq_list[np.where(snp_freq_list[:, 0] != snp_freq_list[:, 1])]
    return calc_thetabeta_unfolded(SNPs, core_freq, n, p) - calc_theta_d(snp_freq_list, c, n)


def calc_var_theta_d(c, n, theta):
    """Calculates the variance of Theta_S

    Parameters:
        c: Speciation time in coalescent units
        n: Sample Size
        theta: genome-wide estimate of the mutation rate
    """
    i = np.arange(2, n+1)
    x = np.sum(1./(i**2.*(i-1)**2.))
    return (1./(c+1./n))**2.*(theta**2.+c*theta+theta/n+theta**2.*x)


def calc_t_b2(snp_freq_list, core_freq, c, n, p, theta, var_dic):
    '''

    Parameters:
        core_freq: freq of SNP under consideration, ranges from 1 to sample size
        n: sample size of core SNP
        p: the p parameter specifying sharpness of peak
        theta: genome-wide estimate of the mutation rate
    '''
    notSubsList_noCore = snp_freq_list[np.where(snp_freq_list[:, 0] != snp_freq_list[:, 1])]
    thetaB = calc_thetabeta_unfolded(notSubsList_noCore, core_freq/n, n, p)
    thetasubs = calc_theta_d(snp_freq_list, c, n)
    if not (n, core_freq, theta) in var_dic:
        VarD = calc_var_theta_d(c, n, theta)
        VarB = calc_var_theta(n, theta, core_freq, p, False)
        denom = math.sqrt(VarD+VarB)
        var_dic[(n, core_freq, theta)] = denom
    else:
        denom = var_dic[(n, core_freq, theta)]
    return (thetaB-thetasubs)/denom


def calc_d(freq, x, p):
    """Calculates the value of d, the similarity measure

    Parameters:
        freq: freq of SNP under consideration, ranges from 0 to 1
        x: freq of coresite, ranges from 0 to 1
        p: the p parameter specifying sharpness of peak
    """
    xf = min(x, 1.-x)
    f = np.minimum(freq, 1.-freq)
    maxdiff = np.maximum(xf, .5-xf)
    corr = ((maxdiff-np.abs(xf-f))/maxdiff)**p
    return corr 


def calc_t_unfolded(snp_freq_list, core_freq, snp_n, p, theta, var_dic):
    """
    Using equation 8 from Achaz 2009

    Parameters:
        core_freq: freq of SNP under consideration, ranges from 1 to sample size
        snp_n: sample size of core SNP
        p: the p parameter specifying sharpness of peak
        theta: genome-wide estimate of the mutation rate
    """

    x = float(core_freq)/snp_n

    num = np.sum(snp_freq_list[:, 0]/snp_freq_list[:, 1]*snp_n*omegai(snp_freq_list[:, 0]/snp_freq_list[:, 1],
                                                                      snp_n, x, p))
    if not (snp_n, core_freq, theta) in var_dic:
        denom = math.sqrt(an(snp_n, x, p) * theta + Bn(snp_n, x, p) * theta**2.)
        var_dic[(snp_n, core_freq, theta)] = denom
    else:
        denom = var_dic[(snp_n, core_freq, theta)]
    return num/denom


def calc_var_theta(n, theta, core_freq, p, wattersons):
    """
    Calculates variance of a given estimator of theta, eq 7 from Achaz.

    Parameters:
        core_freq: freq of SNP under consideration, ranges from 1 to sample size
        n: sample size of core SNP
        p: the p parameter specifying sharpness of peak
        theta: genome-wide estimate of the mutation rate
        wattersons: whether to calculate wattersons theta instead of
    """
    wVector = None
    if wattersons:
        wVector = 1./np.arange(1, n)
    else:
        wVector = calc_d(np.arange(1, n)/float(n), float(core_freq)/n, p)
    t1 = np.sum(wVector)**(-2.)
    t2 = theta*np.sum(wVector**2. * np.arange(1, n))

    i = np.arange(1, n)
    s1 = np.sum(wVector**2*i**2*sigma(n, np.column_stack([i, i])))

    coords = np.asarray([(j, i) for i in range(1, n) for j in range(i+1, n)])
    iind = np.asarray([i-1 for i in range(1, n) for j in range(i+1, n)])
    jind = np.asarray([j-1 for i in range(1, n) for j in range(i+1, n)])

    s2 = np.sum(coords[:, 0] * coords[:, 1] * wVector[iind] * wVector[jind] * sigma(n, coords))

    t3 = theta**2.*(s1+2.*s2)
    return t1*(t2+t3)


def calc_var_theta_fold(n, theta, core_freq, p):
    """
        Parameters:
        core_freq: freq of SNP under consideration, ranges from 1 to sample size
        n: sample size of core SNP
        p: the p parameter specifying sharpness of peak
        theta: genome-wide estimate of the mutation rate
    """

    wVector = calc_d(np.arange(1, int(n/2)+1)/float(n), float(core_freq)/n, p)
    r = np.arange(1, int(n/2)+1)
    t1 = sum(wVector*(1./r+1./(n-r)) * 1./(1+(r == n-r)))**-2.
    t2 = sum([wVector[i-1]**2.*(phi(n, i)*theta+rho_p_ii(n, i)*theta**2.) for i in range(1, int(n/2)+1)])

    coords = np.asarray([(j, i) for i in range(1, int(n/2)+1) for j in range(1, i)])
    t3 = np.sum(wVector[coords[:, 0]-1]*wVector[coords[:, 1]-1] * rho_p_ij(n, coords[:, 0], coords[:, 1]) * theta**2.)

    return t1*(t2+2.*t3)


def calc_cov_folded(n, theta, core_freq, p):
    """
    Parameters:
        core_freq: freq of SNP under consideration, ranges from 1 to sample size
        n: sample size of core SNP
        p: the p parameter specifying sharpness of peak
        theta: genome-wide estimate of the mutation rate
    """
    r = np.arange(1, int(n/2)+1)
    wVector = calc_d(r/float(n), float(core_freq)/n, p)
    t1 = 1./sum(wVector*(1./r+1./(n-r))*1./(1.+(r == n-r)))
    t2 = 1./sum((1./r+1./(n-r))*1./(1+(r == n-r)))
    coords = np.asarray([(i, j) for i in range(1, int(n/2)+1) for j in range(1, int(n/2)+1)])
    t3 = np.sum(wVector[coords[:, 0]-1]*rho_p_ij(n, coords[:, 0], coords[:, 1]) * theta**2.)
    return t1*t2*t3


def calc_var_folded_beta(n, theta, core_freq, p):
    """
    Parameters:
        n: sample size of core SNP
        theta: genome-wide estimate of the mutation rate
        core_freq: freq of SNP under consideration, ranges from 1 to sample size
        p: the p parameter specifying sharpness of peak
    """
    return calc_var_theta_fold(n, theta, core_freq, p) + calc_var_theta(n, theta, core_freq, p, True) - \
        2. * calc_cov_folded(n, theta, core_freq, p)


def omegai(i, snp_n, x, p):
    """Calculates 9a

    Parameters:
        i:freq of SNP under consideration, ranges between 0 and 1
        snp_n: number of chromosomes used to calculate frequency of core SNP
        x: freq of coresite, ranges from 0 to 1
        p: the p parameter specifying sharpness of peak
    """
    n1num = calc_d(i, x, p)
    n1denom = np.sum(calc_d(np.arange(1., snp_n)/snp_n, x, p))
    n1 = n1num/n1denom
    n2 = (1./(i*snp_n)) / (np.sum(1./np.arange(1., snp_n)))
    return n1 - n2


def phi(n, i):
    """
    Calculates equation 12a of Achaz

    Parameters:
        n:sample size
        i: frequency of SNP, in number of individuals
    """
    return n/((1.+(i == n-i)) * i * (n-i))


def rho_p_ii(n, i):
    """
    Calculates equation 12b of Achaz

    Parameters:
        n:sample size
        i: frequency of SNP, in number of individuals
    """
    return (sigma(n, np.column_stack([i, i]))+sigma(n, np.column_stack([n-i, n-i]))+2.
            * sigma(n, np.column_stack([i, n-i]))) / (1.+(i == (n-i)))**2.


def rho_p_ij(n, i, j):
    """
    Calculates equation 12c of Achaz

    Parameters:
        n:sample size
        i: frequency of SNP, in number of individuals
        j: second frequency
    """
    return (sigma(n, np.column_stack([i, j]))+sigma(n, np.column_stack([i, n-j])) + sigma(n, np.column_stack([n-i, j]))
            + sigma(n, np.column_stack([n-i, n-j]))) / ((1.+(i == n-i)) * (1. + (j == n-j)))


def an(snp_n, x, p):
    """
    Calculates alpha_n from Achaz 2009, eq 9b

        snp_n: Sample size
        x: frequency, ranges from 0 to 1
        p: value of p parameter
    """
    i = np.arange(1, snp_n)
    return np.sum(i * omegai(i/float(snp_n), snp_n, x, p)**2.)


def Bn(snp_n, x, p):
    '''
    Returns Beta_N from Achaz 2009, eq 9c

    Parameters:
        snp_n: Sample size
        x: frequency, ranges from 0 to 1
        p: value of p parameter
    '''

    i = np.arange(1, snp_n)
    n1 = np.sum(i**2.*omegai(i/float(snp_n), snp_n, x, p)**2.*sigma(snp_n, np.column_stack([i, i])))

    coords = np.asarray([(j, i) for i in range(1, snp_n) for j in range(1, i)])
    s2 = np.sum(coords[:, 0]*coords[:, 1]*omegai(coords[:, 0]/float(snp_n), snp_n, x, p)
                * omegai(coords[:, 1]/float(snp_n), snp_n, x, p)*sigma(snp_n, coords))

    n2 = 2.*s2
    return n1+n2


def calc_t_fold(snp_freq_list, core_freq, snp_n, p, theta, var_dic):
    """

    Parameters:
        core_freq: freq of SNP under consideration, ranges from 1 to sample size
        snp_n: sample size of core SNP
        p: the p parameter specifying sharpness of peak
        theta: genome-wide estimate of the mutation rate
    """

    x = float(core_freq)/snp_n
    num = calc_beta_folded(snp_freq_list, x, snp_n, p)
    if not (snp_n, core_freq, theta) in var_dic:
        denom = math.sqrt(calc_var_folded_beta(snp_n, theta, core_freq, p))
        var_dic[(snp_n, core_freq, theta)] = denom
    else:
        denom = var_dic[(snp_n, core_freq, theta)]
    return num/denom


def sigma(n, ij):
    """
    Returns sigma from eq 2 or 3 in Fu 1995

    Parameters:
        n: sample size
        ij: 2-d array of integers with 2 cols and no rows
    """
    np.seterr(all='raise')
    res = np.zeros(ij.shape[0])
    # i must be greater than j
    ij[:, 0], ij[:, 1] = ij.max(axis=1), ij.min(axis=1)  
    ci = np.logical_and(ij[:, 0] == ij[:, 1], ij[:, 0] == n/2)

    # Using eq 2
    if np.any(ci) > 0:
        res[ci] = 2.*((fu_an_vec([n]) - fu_an_vec(ij[ci, 0]))/(float(n)-ij[ci, 0]))-(1./(ij[ci, 0]**2.))

    ci = np.logical_and(ij[:, 0] == ij[:, 1], ij[:, 0] < n/2)
    if np.any(ci) > 0:
        res[ci] = fu_Bn(n, ij[ci, 0]+1)

    ci = np.logical_and(ij[:, 0] == ij[:, 1], ij[:, 0] > n/2)
    if np.any(ci) > 0:
        res[ci] = fu_Bn(n, ij[ci, 0])-1./(ij[ci, 0]**2.)

    # using eq 3
    ci = np.logical_and(ij[:, 0] > ij[:, 1], ij[:, 0]+ij[:, 1] == n)
    if np.any(ci) > 0:
        res[ci] = (fu_an_vec([n])-fu_an_vec(ij[ci, 0]))/(n-ij[ci, 0]) + \
         (fu_an_vec([n])-fu_an_vec(ij[ci, 1]))/(n-ij[ci, 1])
        - (fu_Bn(n, ij[ci, 0]) + fu_Bn(n, ij[ci, 1]+1))/2. - 1./(ij[ci, 0]*ij[ci, 1])

    ci = np.logical_and(ij[:, 0] > ij[:, 1], ij[:, 0]+ij[:, 1] < n)
    if np.any(ci) > 0:
        res[ci] = (fu_Bn(n, ij[ci, 0]+1)-fu_Bn(n, ij[ci, 0]))/2.

    ci = np.logical_and(ij[:, 0] > ij[:, 1], ij[:, 0]+ij[:, 1] > n)
    if np.any(ci) > 0:
        res[ci] = (fu_Bn(n, ij[ci, 1])-fu_Bn(n, ij[ci, 1]+1))/2.-(1./(ij[ci, 0] * ij[ci, 1]))

    return res


def fu_an_vec(n):
    """Calculates a_n from Fu 1995, eq 4"""
    a = np.insert(np.cumsum(1./np.arange(1, np.amax(n))), 0, 0)
    return a[np.asarray(n)-1]  # minus one for sum being only to n-1


def fu_Bn(n, i):
    """Calculates Beta_n(i) from Fu 1995, eq 5"""

    r = 2.0 * n/((n-i+1.)*(n-i)) * (fu_an_vec([n+1])-fu_an_vec(i)) - (2./(n-i))
    return r


def find_local_theta(theta_map, start_i, coordinate):
    """
    Given a numpy array of mutation rates finds the theta corresponding to the window that coordinate is in.
    Starts searching at the prior window index to save time
    """
    for i in range(start_i, theta_map.shape[0]):
        if coordinate < theta_map[i, 1] and coordinate >= theta_map[i, 0]:
            return (theta_map[i, 2], i)
    print(sys.exit("Error: Coordinate " + str(coordinate)+" is found in the SNP input file, but is not in any \
        of the windows in the theta_map file."))


def main():

    # Loads the input parameters given by the user
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="Name of input file with all SNPs", type=str, required=True)
    parser.add_argument("-o", help="Output file", type=str, default="/dev/stdout")
    parser.add_argument("-w", help="Maximum Window Size (in bp) to calculate Beta in for a single test SNP", type=int,
                        default=1000)
    parser.add_argument("-onewin", help="Calculate Beta on window which uses all SNPs in input file instead of using \
                         distance-based window", default=False, action="store_true")
    parser.add_argument("-p", help="Power to raise difference measure by", type=int, default=2)
    parser.add_argument("-fold", help="Use folded SFS version", action="store_true")
    parser.add_argument("-B2", help="Use the Beta2 statistic. Substiution data with an outgroup must be provided.",
                        action="store_true")
    parser.add_argument("-m", help="Minimum folded core SNP frequency, exclusive. Must be between 0 and 0.5.",
                        type=float, default=0)
    parser.add_argument("-std", help="Instead of returning Beta value, return normalized Beta Statistic", default=False,
                        action="store_true")
    parser.add_argument("-theta", help="Estimated genome wide theta value per basepair. Used for calculation of \
                        variance. It's equal to 2*l*N_e*u, where u is the locus neutral mutation rate, Ne is the \
                         effective population size and l is the ploidy", type=float)
    parser.add_argument("-theta_map", help="Filename of map of mutation rates. This file should contain estimated\
                         mutation rates in windows across the genomic area you are applying Beta on.", type=str)
    parser.add_argument("-thetaPerSNP", help="Filename of map of mutation rates. This file should contain estimated\
                         mutation rates around each SNP. This file should be two columns: position and estimated theta\
                         rate.", type=str)
    parser.add_argument("-DivTime", help="Divergence time, in coalescent units, between the two species. Only needed\
                         if using B^(2). This can be estimated using the BALLET software, or you can use prior \
                        estimates for your species of interest. In practice, this value affects power very little, \
                         but will affect the standardized statistic.  To convert from generations (g) to coalescent \
                         units (c), the formula is g=c*Ne*2 where Ne is the effective population size.", type=float)

    args = parser.parse_args()
    output = open(args.o, 'w')

    # Check for valid file format and parameters
    try:
        SNPs = np.loadtxt(args.i, dtype=float)
    except IOError:
        print(sys.exit("Error: Input file cannot be found"))
    except:
        print(sys.exit("Error: Input file in wrong format"))
    if args.m < 0 or args.m > .5:
        print(sys.exit("Error: Parameter m must be between 0 and 0.5."))
    if args.p <= 0:
        print(sys.exit("Error: Parameter p must be positive."))
    if len(SNPs.shape) <= 1:
        print(sys.exit("Error: Because the core SNP is excluded from calculations, there must be at least two SNPs in\
                        the input file."))
    if args.std and args.theta is None and args.theta_map is None and args.thetaPerSNP is None:
        print(sys.exit("Error: In order to normalize Beta statistics, a theta value must be provided using the -theta\
                        or -theta_map flags."))
    if args.onewin and (args.theta_map is not None or args.thetaPerSNP is not None):
        print(sys.exit("Error: onewin and theta_map options are not compatible. onewin clculates the mutation rate in\
                        the given window of arbitrary size"))
    if args.w < 2:
        print(sys.exit("Error: Window size must be 2 bp or above. However, you probably want to use a window size much\
                        larger than 2."))
    if args.std and args.theta_map is None and args.theta <= 0 and args.thetaPerSNP is None:
        print(sys.exit("Error: You must provide an estimate of theta (population-scaled mutation rate) and it must be a\
                        positive value."))
    if args.p > 50:
        print(sys.exit("Error: P is too large. Reduce value to prevent python numerical errors. See manual for more \
                       information."))
    if args.fold and args.B2:
        print(sys.exit("Error: You cannot use both B1* (folded Beta) and B2. B1* is for when you have no outgroup, \
                       and B2 is for when you can call substiutions with an outgroup. See manual for guidance about \
                       which to use."))
    if args.DivTime is not None and args.DivTime > 1000:
        print(sys.exit("Error: Your divergence time seems very high. Divergence time should be in coalescent units,\
                        not generations or years."))
    if args.B2 and not np.any(SNPs[:, 1] == SNPs[:, 2]):
        print(sys.exit("Error: You chose to calculate Beta2, but your input file contains no substiutions. If you do \
                        not have substiution data, please use Beta1 or Beta1*."))
    if args.B2 and args.DivTime is None:
        print(sys.exit("You must provide a divergence time using the -DivTime flag to use B2"))
    if args.theta_map is not None and args.thetaPerSNP is not None:
        print(sys.exit("You can use -theta_map or -thetaPerSNP but not both."))

    if args.onewin:
        if args.fold:
            output.write("Position\tBeta1*_std\n")
        elif args.B2:
            output.write("Position\tBeta2_std\n")
        else:
            output.write("Position\tBeta1_std\n")
    elif not args.std and args.fold:
        output.write("Position\tBeta1*\n")
    elif args.std and args.fold:
        output.write("Position\tBeta1*\tBeta1*_std\n")
    elif args.std and not args.B2:
        output.write("Position\tBeta1\tBeta1_std\n")
    elif not args.B2:
        output.write("Position\tBeta1\n")
    elif args.B2 and not args.std:
        output.write("Position\tBeta2\n")
    else:
        output.write("Position\tBeta2\tBeta2_std\n")

    if not args.B2 and np.any(SNPs[:, 1] == SNPs[:, 2]):
        SNPs = SNPs[(SNPs[:, 1] != SNPs[:, 2]) & (SNPs[:, 1] != 0)]

    prev_start_i = 0
    prev_end_i = 0
    var_dic = {}  # records variance calculations so don't need to be recalculated
    theta_map = None
    if args.theta_map is not None:
        theta_map = np.loadtxt(args.theta_map, dtype=float)
    elif args.thetaPerSNP is not None:
        theta_map = np.loadtxt(args.thetaPerSNP, dtype=float)

    curr_theta_i = 0

    if args.onewin:
        theta = calc_thetaw_unfolded(SNPs[:, 1:], int(np.mean(SNPs[:, 2])))
        for snp_i in range(len(SNPs)):
            loc = SNPs[snp_i, 0]
            if len(SNPs) == 1:
                T = 0
                output.write(str(loc)+"\t"+str(round(T, 6))+"\n")
                break

            freqCount = float(SNPs[snp_i, 1])
            sample_n = int(SNPs[snp_i, 2])
            freq = freqCount/sample_n
            SNPSet = np.delete(SNPs, snp_i, axis=0)[:, 1:]
            if int(freqCount) != sample_n and freq < 1.0-args.m and freq > args.m and sample_n > 3:
                if args.fold:
                    T = calc_t_fold(SNPSet, freqCount, sample_n, args.p, theta, var_dic)
                elif args.B2:
                    T = calc_t_b2(SNPSet, freqCount, args.DivTime, sample_n, args.p, theta, var_dic)
                else:
                    T = calc_t_unfolded(SNPSet, freqCount, sample_n, args.p, theta, var_dic)
                output.write(str(loc)+"\t"+str(round(T, 6))+"\n")
            elif freq > 1.0 or freq < 0:
                print(sys.exit("Error: Input file contains SNP of invalid frequency on line "+str(snp_i)+"."))
            elif freq < 1.0-args.m and freq > args.m and sample_n <= 3:
                print(sys.exit("Error: Sample size must be greater than 3 haploid individuals to make inference,\
                                or else theta_beta will always equal theta_watterson's. You may wish to increase\
                                the m paramter value to exclude this SNP from being a core SNP."))
    else:
        for snp_i in range(len(SNPs)):
            loc = int(SNPs[snp_i, 0])
            freqCount = float(SNPs[snp_i, 1])
            sample_n = int(SNPs[snp_i, 2])
            freq = freqCount/sample_n

            if int(freqCount) != sample_n and freq < 1.0-args.m and freq > args.m and sample_n > 3:
                sI, endI = find_win_indx(prev_start_i, prev_end_i, snp_i, SNPs, args.w)
                prev_start_i = sI
                prev_end_i = endI
                B = None
                T = None
                if endI > sI:

                    SNPSet = np.take(SNPs, list(range(sI, snp_i))+list(range(snp_i+1, endI+1)), axis=0)[:, 1:]
                    if args.fold:
                        B = calc_beta_folded(SNPSet, freqCount/sample_n, sample_n, args.p)
                    elif not args.fold and not args.B2:
                        B = calc_beta_unfolded(SNPSet, freqCount/sample_n, sample_n, args.p)
                    elif args.B2:
                        B = calc_beta_2(SNPSet, args.DivTime, sample_n, freqCount/sample_n, args.p)

                    if args.theta_map is not None or args.thetaPerSNP is not None:
                        theta = None
                        if args.thetaPerSNP is not None:
                            theta = theta_map[np.where(theta_map[:, 0] == int(loc)), 1] 
                            if len(theta[0]) == 1:
                                theta = float(theta)
                            elif len(theta[0]) > 1:
                                theta = float(theta[0][0])
                            else:
                                print(sys.exit("SNP at location "+str(loc)+" is not in thetaPerSNP file or is found \
                                               more than once"))
                        else:
                            theta, curr_theta_i = find_local_theta(theta_map, curr_theta_i, loc)
                        if args.fold:
                            T = calc_t_fold(SNPSet, freqCount, sample_n, args.p, theta * args.w, var_dic)
                        elif args.B2:
                            T = calc_t_b2(SNPSet, freqCount, args.DivTime, sample_n, args.p, theta*args.w, var_dic)
                        else:
                            T = calc_t_unfolded(SNPSet, freqCount, sample_n, args.p, theta*args.w, var_dic)
                    elif args.std:
                        if args.fold:
                            T = calc_t_fold(SNPSet, freqCount, sample_n, args.p, args.theta * args.w, var_dic)
                        elif args.B2:
                            T = calc_t_b2(SNPSet, freqCount, args.DivTime, sample_n, args.p, args.theta*args.w, var_dic)
                        else:
                            T = calc_t_unfolded(SNPSet, freqCount, sample_n, args.p, args.theta * args.w, var_dic)

                if endI == sI:
                    B = 0
                    T = 0
                if not args.std:
                    output.write(str(loc)+"\t"+str(round(B, 6))+"\n")  # Remove thetas
                else:
                    output.write(str(loc)+"\t"+str(round(B, 6))+"\t"+str(round(T, 6))+"\n")
            elif freq > 1.0 or freq < 0:
                print(sys.exit("Error: Input file contains SNP of invalid frequency on line "+str(snp_i)+"."))
            elif freq < 1.0-args.m and freq > args.m and sample_n <= 3:
                print(sys.exit("Error: Sample size must be greater than 3 haploid individuals to make inference, \
                               or else theta_beta will always equal theta_watterson's. You may wish to increase the \
                               m paramter value to exclude this SNP from being a core SNP."))


if __name__ == "__main__":
    main()
