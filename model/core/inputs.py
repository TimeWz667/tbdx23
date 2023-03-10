from scipy.interpolate import interp1d
from functools import lru_cache
from collections import namedtuple
import json
import pandas as pd
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['load_inputs', 'load_k_covid']


Inputs = namedtuple("Inputs", ('Demography', 'TxOut', 'DyHIV'))


class Demography:
    def __init__(self, src, year0=2010):
        self.Source = src
        self.Years = src['years']
        self.Year0 = year0

        self.YearRange = [min(self.Years), max(self.Years)]

        self.RateDeath = interp1d(self.Years, src['dr'])
        self.RateBirth = interp1d(self.Years, src['br'])
        self.N = interp1d(self.Years, src['N'])
        self.N0 = self.N(year0)

    @lru_cache(maxsize=1024)
    def __call__(self, time):
        if time < self.YearRange[0]:
            time = self.YearRange[0]
        elif time > self.YearRange[1]:
            time = self.YearRange[1]

        return {
            'Year': time,
            'r_birth': float(self.RateBirth(time)),
            'r_die': float(self.RateDeath(time))
        }


class TxOut:
    def __init__(self, src):
        self.PrTxDie = src['Death']
        self.PrTxSucc = src['Succ']
        self.PrTxLTFU = 1 - self.PrTxDie - self.PrTxSucc

    def __repr__(self):
        return f'Successful: {self.PrTxSucc:.0%}, LTFU: {self.PrTxLTFU:.0%}, Die: {self.PrTxDie:.0%}'


class HIV:
    def __init__(self, src):
        self.SRC = src
        self.R_ART0, self.R_ART1 = src['r_art0'], src['r_art1']
        self.RT_ART, self.T0_ART = src['rt_art'], src['t0_art']

        self.R_HIV0, self.RT_HIV = src['r_hiv0'], src['rt_hiv']

        self.DR_HIV0, self.DR_HIV1, self.DRT_HIV = src['dr_hiv0'], src['dr_hiv1'], src['drt_hiv']

        self.HIV0 = np.array(src['y0'])

    def __call__(self, time):
        if time < 2010:
            time = 2010

        return {
            'r_art': self.R_ART0 + (self.R_ART1 - self.R_ART0) / (1 + np.exp(- self.RT_ART * (time - self.T0_ART))),
            'r_hiv': self.R_HIV0 * np.exp(- self.RT_HIV * (time - 2010)),
            'r_die_hiv': self.DR_HIV0 * np.exp(- self.DRT_HIV * (time - 2010)) + self.DR_HIV1
        }


def load_k_covid(filepath):
    df = pd.read_csv(filepath)
    return interp1d(df.Time, df.k_covid, kind='nearest', bounds_error=False, fill_value=1)


def load_inputs(root):
    with open(f'{root}/pars_pop.json', 'r') as f:
        src = json.load(f)

    demo = Demography(src)

    with open(f'{root}/pars_tx.json', 'r') as f:
        src = json.load(f)


    tos = {k: TxOut(s) for k, s in src.items()}
    to_all = TxOut(src['All'])

    try:
        with open(f'{root}/pars_hiv.json', 'r') as f:
            src = json.load(f)
        hiv = HIV(src)
    except FileNotFoundError:
        hiv = None

    return Inputs(demo, tos, hiv)


if __name__ == '__main__':
    with open('../../pars/IND/pars_pop.json', 'r') as f:
        src = json.load(f)

    demo = Demography(src)

    print('Year: 1990')
    print(demo(1990))

    print('Year: 2030')
    print(demo(2030))

    with open('../../pars/IND/pars_tx.json', 'r') as f:
        src = json.load(f)

    inp = load_inputs('../../pars/ZAF')
    print(inp)

    with open('../../pars/ZAF/pars_hiv.json', 'r') as f:
        src = json.load(f)

    hiv = HIV(src)

    for t in range(2010, 2030):
        s = hiv(t)
        print(t, s['r_hiv'], s['r_art'])
