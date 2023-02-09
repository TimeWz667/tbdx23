import numpy as np
from numba import njit

__author__ = 'Chu-Chang Ku'
__all__ = ['calc_dy', 'extract_tr']


@njit
def calc_dy0(y, frs, tos, rates):
    dy = np.zeros_like(y)

    for i in range(len(frs)):
        fr, to, rate = frs[i], tos[i], rates[i]
        if rate > 0:
            n_tr = rate * y[fr]
            dy[fr] -= n_tr
            dy[to] += n_tr
    return dy


def calc_dy(y, trs):
    frs = np.array([tr[0] for tr in trs], int)
    tos = np.array([tr[1] for tr in trs], int)
    rates = np.array([tr[2] for tr in trs])
    return calc_dy0(y, frs, tos, rates)


@njit
def extract_tr0(y, frs, rates):
    n_tr = 0
    for i in range(len(frs)):
        fr, rate = frs[i], rates[i]
        if rate > 0:
            n_tr += rate * y[fr]
    return n_tr


def extract_tr(y, trs, fil):
    trs = [tr for tr in trs if fil(tr)]
    frs = np.array([tr[0] for tr in trs], int)
    rates = np.array([tr[2] for tr in trs])
    return extract_tr0(y, frs, rates)

