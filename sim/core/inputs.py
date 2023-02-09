from scipy.interpolate import interp1d
from functools import lru_cache
from collections import namedtuple
import json

__author__ = 'Chu-Chang Ku'
__all__ = ['load_inputs']


Inputs = namedtuple("Inputs", ('Demography', 'TxOut', 'TxOutHIV'))


class Demography:
    def __init__(self, src):
        self.Source = src
        self.Years = src['years']
        self.Year0 = src['Year0']

        self.YearRange = [min(self.Years), max(self.Years)]

        self.N0 = src['N0']

        self.RateDeath = interp1d(self.Years, src['dr'])
        self.RateBirth = interp1d(self.Years, src['br'])

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
        self.PrTxLTFU = src['LTFU']
        self.PrTxSucc = src['Succ']

    def __repr__(self):
        return f'Successful: {self.PrTxSucc:.0%}, LTFU: {self.PrTxLTFU:.0%}, Die: {self.PrTxDie:.0%}'


def load_inputs(root):
    with open(f'{root}/pars_pop.json', 'r') as f:
        src = json.load(f)

    demo = Demography(src)

    with open(f'{root}/pars_tx.json', 'r') as f:
        src = json.load(f)

    to_all = TxOut(src['All'])

    if 'PLHIV' in src:
        to_hiv = {
            'PLHIV': TxOut(src['PLHIV']),
            'NonHIV': TxOut(src['NonHIV'])
        }
    else:
        to_hiv = None

    return Inputs(demo, to_all, to_hiv)


if __name__ == '__main__':
    import json

    with open('../../data/pars/IND/pars_pop.json', 'r') as f:
        src = json.load(f)

    demo = Demography(src)

    print('Year: 1990')
    print(demo(1990))

    print('Year: 2030')
    print(demo(2030))

    with open('../../data/pars/IND/pars_tx.json', 'r') as f:
        src = json.load(f)

    inp = load_inputs('../../data/pars/ZAF')
    print(inp)
