import numpy.random as rd
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['Cascade', 'RepoCascade']


class Cascade:
    def __init__(self, src, year0=2010):
        self.Source = src
        self.Year0 = year0
        self.YearSurveyed = src['Year0']
        self.Prev = self.calc_prev()

        self.R_Onset = src['r_sym']
        self.R_CSI = src['r_aware0']
        self.Det = src['r_det0'] * self.Prev['ExCS']

        self.R_SelfCure = src['r_sc']
        self.R_Death_A = src['r_death_a']
        self.R_Death_S = src['r_death_s']

        self.R_ReCSI = src['r_det0']

        self.RR_Det_t = src['rr_det_t']
        self.PrDx0 = 0
        self.PrDx1 = 1

        self.PrUnder = src['p_under']
        self.PPV = src['ppv']

    def calc_k_det(self, t):
        if t >= self.Year0:
            return np.power(self.RR_Det_t, t - self.YearSurveyed)
        else:
            return np.power(self.RR_Det_t, self.Year0 - self.YearSurveyed)

    def calc_prev(self):
        src = self.Source
        adr = src['adr']
        r_sc, r_death_a, r_death_s, r_death_bg = src['r_sc'], src['r_death_a'], src['r_death_s'], src['r_death_bg']
        r_sym = src['r_sym']
        r_aware = src['r_aware0'] * np.power(src['rr_det_t'], self.Year0 - src['Year0'])
        r_det = src['r_det0'] * np.power(src['rr_det_t'], self.Year0 - src['Year0'])

        inc = src['inc0'] * np.exp(- src['adr'] * (self.Year0 - src['Year0']))

        ra = r_sc + r_death_a + r_death_bg
        rs = r_sc + r_death_s + r_death_bg
        rc = r_sc + r_death_s + r_death_bg

        prv_a = inc / (r_sym + ra - adr)
        prv_s = prv_a * r_sym / (r_aware + rs - adr)
        prv_c = prv_s * r_aware / (r_det + rc - adr)

        return {
            'Asym': prv_a,
            'Sym': prv_s,
            'ExCS': prv_c
        }

    def reform(self, pdx0=0, pdx1=1):
        pdx0 = min(pdx0, self.Det / (self.R_CSI * self.Prev['Sym']))

        self.PrDx0, self.PrDx1 = pdx0, pdx1

        det0 = self.R_CSI * pdx0 * self.Prev['Sym']
        det1 = self.Det - det0
        self.R_ReCSI = det1 / (pdx1 * self.Prev['ExCS'])

    def revert(self):
        self.reform(0, 1)

    def __repr__(self):
        return f'Cascade({self.R_Onset:.3f} -> {self.R_CSI:.3f}({self.PrDx0:.0%}) -> {self.R_ReCSI:.3f}({self.PrDx1:.0%}))'


class RepoCascade:
    def __init__(self, src, year0=2010):
        self.Cascades = [Cascade(dict(row), year0=year0) for _, row in src.iterrows()]

    def sample(self):
        return rd.choice(self.Cascades, 1)[0]


if __name__ == '__main__':
    import pandas as pd

    pars_cs = pd.read_csv('../../results/pars_IND.csv')
    repo_cs = RepoCascade(pars_cs)

    cas = repo_cs.sample()
    print(cas)

    for t in [2005, 2010, 2015, 2020, 2025]:
        print(t, ': ', cas.calc_k_det(t))
