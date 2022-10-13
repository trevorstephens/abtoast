"""
abtoast - A/B testing tools.
"""

# Author: Trevor Stephens <trevorstephens.com>
#
# License: BSD 3 clause

import math
from scipy.stats import norm

__all__ = ['power_prop_test']


def power_prop_test(p1, p2, power=0.8, significance=0.05, tside=2, k=1.0):
    # Bernard Rosner - Fundamentals of Biostatistics.
    # Where n2 is k times as large as n1
    # Returns required sample size in *each* group
    q1, q2 = 1. - p1, 1. - p2
    p_bar = (p1 + (k * p2)) / (1. + k)
    q_bar = 1. - p_bar
    delta = abs(p2 - p1)
    z_alpha = norm.isf(1 - (significance / tside))
    z_beta = norm.isf(power) # power = 1 - beta
    lhs = math.sqrt(p_bar * q_bar * (1. + (1. / k))) * z_alpha
    rhs = math.sqrt(p1 * q1 + (p2 * q2 / k)) * z_beta
    n = ((lhs + rhs) ** 2) / delta ** 2
    if k == 1:
        return math.ceil(n)
    return math.ceil(n), math.ceil(k * n)
