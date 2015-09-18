"""Classes implementing Lp and Lr."""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

from math import exp, log, factorial
from operator import itemgetter
from ersa.chisquare import LL_ratio_test, likelihood_ratio_CI
import inflect

class Background:
    """
    Class Background
    ----------------
    Represents a calculator for the null hypothesis, Lp.

    Created with empirical parameters, t, theta, lambba_.

    Provides the likelihood function, L(n, s), and log-likelihood
    function, LL(n, s).


    Example
    -------
    n = 5                  # number of shared segments
    s = [3, 5, 4, 9, 4.5]  # list of shared segment lengths

    h0 = Background(2.5, 3.12, 15)
    h0.LL(n, s)


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
    """
    def __init__(self, t, theta, lambda_):
        self.t = t
        self.theta = theta
        self.lambda_ = lambda_

    def _Fp(self, i):
        assert i >= self.t
        prob = exp(-(i - self.t) / (self.theta - self.t)) / (self.theta - self.t)
        return log(prob)

    def _Sp(self, s):
        result = 0
        for i in s:
            result += self._Fp(i)
        return result

    def _Np(self, n):
        l_prob = n * log(self.lambda_) - self.lambda_ - log(factorial(n))
        return l_prob

    def LL(self, n, s):
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


    Example
    -------
    n = 5                  # number of shared segments
    s = [3, 5, 4, 9, 4.5]  # list of shared segment lengths

    ha = Relation(22, 35.3, 2.5, 3.12, 14)
    ha.MLL(n, s)


    See also
    -------
    Class Background


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
    """
    def __init__(self, c, r, t, theta, lambda_):
        super(Relation, self).__init__(t, theta, lambda_)
        self.c = c
        self.r = r
        self.a = 2  # see Huff et al 2011 supplemental material

    def _Fa(self, i, d):
        assert i >= self.t
        l_prob = (-d * (i-self.t) / 100) - log(100 / d)
        return l_prob

    def _Sa(self, s, d):
        result = 0
        for i in s:
            result += self._Fa(i, d)
        return result

    def _p(self, d):
        return exp((-d * self.t) / 100)

    def _Na(self, n, d):
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
        result += self._Sa(s[np:], d)
        return result

    def MLL(self, n, s, d):
        """
        For a given d, return the maximum log-likelihood (MLL).
        Requires s to be sorted smallest to largest.
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
    def __init__(self, pair, dob, d, reject, null_LL, max_LL, lower_d, upper_d, alts, s):
        self.indv1, self.indv2 = pair.split(':')
        self.dob = dob
        self.d = d
        self.reject = reject
        self.alts = alts
        self.s = s
        self.null_LL = null_LL
        self.max_LL = max_LL
        self.lower_d = lower_d
        self.upper_d = upper_d
        if reject and dob[0] and dob[1]:
            self.rel_est = potential_relationship(self.d, self.indv1, self.indv2, dob[0], dob[1])
        else:
            self.rel_est = None


def estimate_relation(pair, dob, n, s, h0, ha, max_d, alpha):
    """
    Tests a pair of individuals for a relation and
    returns a structure holding relevant parameters.

    Requires s to a pre-sorted list of shared segments,
    from smallest to largest.  h0 and ha must be Background and
    Relation objects, respectively.
    """
    assert isinstance(h0, Background)
    assert isinstance(ha, Relation)
    for i in range(1, len(s)):
        assert s[i - 1] <= s[i]

    null_LL = h0.LL(n, s)

    alts = []
    for d in range(1, max_d + 1):
        alt_np, alt_MLL = ha.MLL(n, s, d)
        alts.append((d - 1, alt_np, alt_MLL))  # subtract one from d since a = 2
    max_alt = max(alts, key=itemgetter(2))
    d = max_alt[0]
    max_LL = max_alt[2]

    reject = LL_ratio_test(max_LL, null_LL, alpha)
    lower_d, upper_d = 0, 0
    if reject:
        lower_d, upper_d = likelihood_ratio_CI(alts, max_LL, alpha)

    est = Estimate(pair, dob, d, reject, null_LL, max_LL, lower_d, upper_d, alts, s)
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
    """
    inflect_eng = inflect.engine()
    s = inflect_eng.number_to_words(n)
    excepts = {1: "once", 2: "twice", 3: "thrice"}
    if n in excepts:
        s = excepts[n]
    if capitalize:
        s = s.capitalize()
    return s


def _build_rel_map():
    """
    Helper function for potential_relationship() that
    builds the static relationship/consanguinity map.

    References
    ----------
    https://en.wikipedia.org/wiki/File:Table_of_Consanguinity_showing_degrees_of_relationship.png
    """
    rel_map = {1: {-1: "Parent", 1: "Child"},
               2: {-2: "Grandparent", 0: "Sibling", 2: "Grandchild"},
               3: {-3: "Great Grandparent", -1: "Aunt/Uncle", 1: "Niece/Nephew", 3: "Great Grandchild"}}
    for d in range(4, 7 + 1):
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


@_static_vars(rel_map=_build_rel_map())
def potential_relationship(d_est, indv1, indv2, dob1, dob2):
    """
    Estimates a potential consanguinity between two individuals,
    based on assumption 30 years between generations, centered
    at dob1.

    Returns
    -------
    A tuple with the first value from the perspective
    of indv1 and the second value from the perspective of
    indv2.

    If a relationship cannot be estimated, None is returned.
    """
    yr_per_gen = 30
    delta = dob2 - dob1
    gen_bin = (delta + yr_per_gen / 2) // yr_per_gen
    if d_est not in potential_relationship.rel_map:
        return None
    if gen_bin not in potential_relationship.rel_map[d_est]:
        return None
    else:
        bin_map = potential_relationship.rel_map[d_est]
    return bin_map[gen_bin], bin_map[-gen_bin]


