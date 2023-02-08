from scipy.interpolate import interp1d
from functools import lru_cache

__author__ = 'Chu-Chang Ku'
__all__ = ['Demography']


class Demography:
    def __init__(self, src):
        self.Source = src
        self.Years = src['years']
        self.Year0 = src['Year0']

        self.YearRange = [min(self.Years), max(self.Years)]

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


if __name__ == '__main__':
    import json

    with open('../../data/IND/pars_pop.json', 'r') as f:
        src = json.load(f)

    demo = Demography(src)

    print('Year: 1990')
    print(demo(1990))

    print('Year: 2030')
    print(demo(2030))
