"""
Classes implementing Lp and Lr and methods
to estimate relationships.
"""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

from math import exp, log, factorial, log1p
from operator import itemgetter
import inflect
from ersa.chisquare import LL_ratio_p, likelihood_ratio_CI
from ersa.mask import total_masked

class Background:
    """
    Class Background
    ----------------
    Represents a calculator for the null hypothesis, Lp.

    Created with empirical parameters, t, theta, lambba_.

    Provides the likelihood function, L(n, s), and log-likelihood
    function, LL(n, s).

    Parameters
    ----------
    t : float
        minimum threshold for segment lengths (in cM)

    theta : float
        mean shared segment length (in cm) in the population
        for all segments of size >t and <h

    lambda_ : float
        the sample mean of the number of segments shared
        in the population

    Notes
    -----
    Lp(n,s|t) = Np(n|t) * Sp(s|t)

        Np(n|t) = Poisson distribution with mean equal to sample mean
                  of the number of segments shared in the population (mu)
        Sp(s|t) = Product over i in s of Fp(i|t)
        Fp(i|t) = (e^(i-t)/theta)/theta

        n = number of shared segments
        s = list of shared segment lengths (in cM)
        t = minimum threshold for segment lengths (in cM)
        h = maximum threshold for segment lengths (in cM)
        lambda_ = the sample mean of the number of segments shared
                  in the population
        theta = mean shared segment length (in cm) in the population
                for all segments of size >t and <h

    Example
    -------
    n = 5                  # number of shared segments
    s = [3, 5, 4, 9, 4.5]  # list of shared segment lengths

    h0 = Background(2.5, 3.12, 15)
    h0.LL(n, s)
    """
    def __init__(self, t, theta, lambda_):
        self.t = t
        self.theta = theta
        self.lambda_ = lambda_

    def _Fp(self, i):
        assert i >= self.t
        l_prob = -(i - self.t) / self.theta - log(self.theta)
        return l_prob

    def _Sp(self, s):
        result = 0
        for i in s:
            result += self._Fp(i)
        return result

    def _Np(self, n):
        l_prob = n * log(self.lambda_) - self.lambda_ - log(factorial(n))
        return l_prob

    def LL(self, n, s):
        """
        Compute the log-likelihood of the background.

        Parameters
        ----------
        n : int
            number of shared segments

        s : list[float]
            list of segment lengths (in cM)

        Returns
        -------
        ret: float
           log-likelihood
        """
        ret = 0
        ret += self._Np(n)
        ret += self._Sp(s)
        return ret


class Relation(Background):
    """
    Class Relation
    --------------
    Represents a calculator for the alternative hypothesis, Lr.

    Created with empirical parameters c, r, t, theta, and lambda_. The
    latter three empirical parameters are described in Class Background.

    Provides the a function returning the maximum log-likelihood over
    all possible np, MLL(n, s).  Create this class for each d evaluated.

    Parameters
    ----------
    r : float
        expected number of recombination events per haploid genome
        per generation
    c : int
        number of autosomes

    t : float
        minimum threshold for segment lengths (in cM)

    theta : float
        mean shared segment length (in cm) in the population
        for all segments of size >t and <h

    lambda_ : float
        the sample mean of the number of segments shared
        in the population

    first_deg_adj : bool
        Controls whether first degree adjustment equations should be used

    avuncular_adj : bool
        Use Li et al (2014) Equation 9 for (a=2, d=3) relationships

    Notes
    -----
    Lr = La(na, sa | d, a, t) * Lp(np, sp | t)

        Lp = follows description of the null hypothesis (see
             Class Background)
        La(na, sa | d, a, t) = Na(n | d, a, t) * Sa(sa | d, t)
        Sa(s | d, t) = product over i in s of Fa(i | t)
        Fa(i | d, t) = (e^(-d(i-t)/100)) / (100/d0)
        Na(n | d, a, t) = <eq. 8 in Huff et. al. 2011>

        a = number of ancestors shared (set to 2 based on supp. mat.)
        d = combined number of generations separating the individuals
            from their ancestor(s)
        n = number of shared segments
        na = number of shared segments inherited from recent ancestors
        np = number of shared segments shared due to the population background
        s = list of shared segment lengths (in cM)
        sa and sp = mutually exclusive subsets of s, with sa equal to a list
                    of segment lengths inherited from recent ancestors and
                    sp equal to a list of segment lengths shared due to the
                    background

        r = expected number of recombination events per haploid genome
            per generation
        c = number of autosomes
        p(t) = probability that a shared segment is longer than t

    Example
    -------
    n = 5                  # number of shared segments
    s = [3, 5, 4, 9, 4.5]  # list of shared segment lengths

    ha = Relation(22, 35.3, 2.5, 3.12, 14)
    ha.MLL(n, s)

    See also
    -------
    Class Background
    """
    def __init__(self, c, r, t, theta, lambda_, first_deg_adj=False, nomask=False, avuncular_adj=False):
        super(Relation, self).__init__(t, theta, lambda_)
        self.c = c
        if nomask:
            self.r = r
        else:
            m = total_masked()
            self.r = r - m / 100
        self.a = 2  # see Huff et al 2011 supplemental material
        self.first_deg_adj = first_deg_adj
        self.avuncular_adj = avuncular_adj

    def _Fa(self, i, d, k_max=200):
        assert i >= self.t
        if self.first_deg_adj and d == 2:
            # Li et al (2014) Eqn 6
            # k = 1
            # l_a1 = k * log(1/2) + (k-1) * log(l) -
            #        (d * l) / 100 - k * log(100/d) - log(factorial(k - 1))
            l_a1 = log(1/2) - (d * i) / 100 - log(100/d)

            sum = 0
            for k in range(2, k_max + 1):
                # Todo speed improvement
                # l_ai = k * log(1/2) + (k-1) * log(l) - (d * l) / 100
                #      - k * log(100/d) - log(factorial(k - 1))
                diff = (k - 1) * log(i / 100)
                term = 0
                for j in range(2, k):
                    term += log(j)
                diff -= term
                sum += exp(diff)
            y = log1p(sum)
            l_prob = l_a1 + y
        else:
            l_prob = (-d * (i - self.t) / 100)
            l_prob += -log(100 / d)
        return l_prob

    def _Sa(self, s, d):
        result = 0
        for i in s:
            fa = self._Fa(i, d)
            result += fa
        return result

    def _p(self, d):
        return exp((-d * self.t) / 100)

    def _Na(self, n, d):
        if self.first_deg_adj and d == 2:
            # Huff et al (2011) Eqn S1 & Li et al (2014) Eqn 7
            lambda_ = (3/4) * self.c + 2 * d * self.r * (3/4) * (1/4)
        elif self.avuncular_adj and d == 3:
            # Li et al (2014) Eqn 9
            lambda_ = (3/4) * self.c + 4 * self.r * ((3/4) * (1/4))
        else:
            lambda_ = (self.a * (self.r * d + self.c) * self._p(d)) / (2 ** (d - 1))
        l_prob = n * log(lambda_) - lambda_ - log(factorial(n))
        return l_prob

    def _LLr(self, np, na, s, d):
        """
        Compute Lr given d and np.

        Requires s to be pre-sorted from smallest to largest.
        """
        result = 0
        result += self._Np(np)
        result += self._Na(na, d)
        result += self._Sp(s[:np])
        sa = self._Sa(s[np:], d)
        result += sa
        return result

    def MLL(self, n, s, d):
        """
        For a given d, return the maximum log-likelihood (MLL) and
        corresponding number of shared segments attributed to
        background. Requires s to be sorted smallest to largest.

        Parameters
        ----------
        n : int
            number of shared segments

        s : list[float]
            list of segment lengths (in cM)

        d : int
            combined number of generations separating the
            individuals from their ancestor(s)

        Returns
        -------
        max_np, max_mll : (int, float)
           number of segments attributed to background, log-likelihood
        """
        max_mll, max_np = None, None
        for np in range(n + 1):
            mll = self._LLr(np, n - np, s, d)
            if max_mll is None:
                max_mll, max_np = mll, np
            elif max_mll < mll:
                max_mll, max_np = mll, np
        return max_np, max_mll


class Estimate:
    """
    Structure to hold results from estimate_relation
    """
    def __init__(self, pair, dob, d, reject, null_LL, max_LL, lower_d, upper_d, alts, s, np, p):
        self.pair = pair
        self.indv1, self.indv2 = pair.split(':')
        self.dob = dob
        self.reject = reject
        self.alts = alts
        self.s = s
        self.null_LL = null_LL
        self.max_LL = max_LL
        self.lower_d = lower_d
        self.upper_d = upper_d
        self.np = np
        self.p = p
        if reject:
            years = [dob[0], dob[1]]
            if dob[0] is None or dob[1] is None:
                if d % 2 == 0:
                    years[0], years[1] = 0, 0
                else:
                    years[0], years[1] = 0, 31
            self.rel_est1, self.rel_est2 = \
                potential_relationship(d, self.indv1, self.indv2, years[0], years[1])
        else:
            self.rel_est1, self.rel_est2 = None, None
        # "collapse" d from number of meiosis to
        # relationship degree. Note that for d > 1
        # this is just a shift, but for d = 1
        # the new d becomes 0, which in reality should
        # still be d = 1 since the algorithm
        # does not test for MZ twins.
        self.d = d - 1
        self.cm = sum(s)
        self.LLs = {}
        for i in range(len(self.alts)):
            alt = self.alts[i]
            self.LLs[alt[0] - 1] = alt[2]


def estimate_relation(pair, dob, n, s, h0, ha, max_d, alpha, ci=False):
    """
    Tests a pair of individuals for a relation.  Requires s to be a
    pre-sorted list of shared segments, from smallest to largest.

    Parameters
    ----------
    pair : str
        Identifier for pair of individuals being tested in the format
        "indv1:indv2"

    dob : (int, int) | (None, None)
        Tuple for years of birth for (indv1, indv2)

    n : int
        Number of shared segments

    s : list[float]
        List of shared segment lengths

    h0 : Background

    ha : Relation

    max_d : int
        Maximum d to test

    alpha : float
        Significance level for likelihood ratio test

    ci : bool
        Controls whether confidence intervals are calculated

    Returns
    -------
    est : Estimate
        Contains parameters related to the overall maximum likelihood
        model, and the log-likelihoods of maximum models for each d.
    """
    assert isinstance(h0, Background)
    assert isinstance(ha, Relation)
    for i in range(1, len(s)):
        assert s[i - 1] <= s[i]

    null_LL = h0.LL(n, s)

    alts = []
    for d in range(1, max_d + 1):
        alt_np, alt_MLL = ha.MLL(n, s, d)
        alts.append((d, alt_np, alt_MLL))
    max_alt = max(alts, key=itemgetter(2))
    d, np, max_LL = max_alt[0], max_alt[1], max_alt[2]
    p = LL_ratio_p(max_LL, null_LL)
    reject = True if p < alpha else False
    lower_d, upper_d = None, None
    if ci and reject:
        lower_d, upper_d = likelihood_ratio_CI(alts, max_LL, alpha)

    est = Estimate(pair, dob, d, reject, null_LL, max_LL, lower_d, upper_d, alts, s, np, p)
    return est


def _static_vars(**kwargs):
    """
    Decorator for adding static variables to functions.

    Reference
    ---------
    Answer by Claudiu and ony at:
    http://stackoverflow.com/questions/279561/what-is-the-python-equivalent-of-static-variables-inside-a-function
    """
    def decorate(func):
        for k in kwargs:
            setattr(func, k, kwargs[k])
        return func
    return decorate


def _n_to_ord(n):
    """
    Converts an integer n to an ordinal number (string).

    Parameters
    ----------
    n : int

    Returns
    -------
    str(n) + suffix : str
    """
    assert type(n) == int
    suffix = "th"
    suffixes = {1: "st", 2: "nd", 3: "rd"}
    i = n if n < 20 else n % 10
    if i in suffixes:
        suffix = suffixes[i]
    return str(n) + suffix


def _n_to_w(n, capitalize=True):
    """
    Converts an integer n to a (capitalized) word.

    Parameters
    ----------
    n : int

    capitalize : bool
        Controls whether to capitalize the word

    Returns
    -------
    s : str
    """
    inflect_eng = inflect.engine()
    s = inflect_eng.number_to_words(n)
    excepts = {1: "once", 2: "twice", 3: "thrice"}
    if n in excepts:
        s = excepts[n]
    if capitalize:
        s = s.capitalize()
    return s


def _build_rel_map(dmax=20):
    """
    Helper function for potential_relationship() that
    builds the static relationship/consanguinity map.

    Returns
    -------
    rel_map : dict[int, dict[int, str]]

    References
    ----------
    https://en.wikipedia.org/wiki/File:Table_of_Consanguinity_showing_degrees_of_relationship.png
    """
    rel_map = {0: {0: "Identical Twins or Duplication"},
               1: {-1: "Parent", 1: "Child"},
               2: {-2: "Grandparent", 0: "Sibling", 2: "Grandchild"},
               3: {-3: "Great Grandparent", -1: "Aunt/Uncle", 1: "Niece/Nephew", 3: "Great Grandchild"}}
    for d in range(4, dmax + 1):
        gen_bin = {}
        if d % 2:
            k = 1
        else:
            k = 2
            gen_bin[0] = _n_to_ord(d // 2 - 1) + " Cousin"
        for i in range(k, d + 1, 2):
            name = ""
            if i == d:
                name = _n_to_ord(i - 2)
                name += " Great Grand"
                name2 = name + "child"
                name += "parent"
            elif i == d - 2:
                if i - 2 > 1:
                    name = _n_to_ord(i - 2)
                    name += " "
                name += "Great "
                if i - 2 > 0:
                    name += "Grand "
                name2 = name + "Niece/Nephew"
                name += "Aunt/Uncle"
            else:
                name = _n_to_ord(d // 2 - 1 - i // 2)
                name += " Cousin "
                name += "{} ".format(_n_to_w(i))
                if i > 3:
                    name += "Times"
                name += "Removed"
                name2 = name
            gen_bin[-i] = name
            gen_bin[i] = name2
        rel_map[d] = gen_bin

    return rel_map


@_static_vars(rel_map=_build_rel_map(dmax=20))
def potential_relationship(d_est, indv1, indv2, dob1, dob2):
    """
    Estimates a potential consanguinity between two individuals,
    based on assumption 30 years between generations, centered
    at dob1.

    Parameters
    ----------
    d_est : int
        Point estimate for combined number of generations separating
        the individuals

    indv1 : str
        Identifier for individual 1

    indv2 : str
        Identifier for individual 2

    dob1 : int
        Year of birth for individual 1

    dob2 :
        Year of birth for individual 2

    Returns
    -------
    bin_map[gen_bin], bin_map[-gen_bin] : (str, str) | None
        Estimated potential consanguinity tuple with the first value
        from the perspective of indv1 and the second value from the
        perspective of indv2. If a relationship cannot be estimated,
        None is returned.
    """
    yr_per_gen = 30
    delta = dob2 - dob1
    if d_est == 0:
        gen_bin = 0
    else:
        gen_bin = (delta + yr_per_gen / 2) // yr_per_gen
    if d_est not in potential_relationship.rel_map:
        return None
    if gen_bin not in potential_relationship.rel_map[d_est]:
        return None
    else:
        bin_map = potential_relationship.rel_map[d_est]
    if gen_bin % 2 == 1 and d_est < 5:
        ret = bin_map[gen_bin]
        ret += " or "
        ret += bin_map[-gen_bin]
        return ret, ret
    return bin_map[gen_bin], bin_map[-gen_bin]


