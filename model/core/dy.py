from abc import ABCMeta, abstractmethod
import numpy as np
from scipy.integrate import solve_ivp
import pandas as pd

__all__ = ['AbsModel']


class AbsModel(metaclass=ABCMeta):
    def __init__(self, inputs):
        self.Inputs = inputs
        self.Year0 = inputs.Demography.Year0
        self.YearBaseline = 2022

    @abstractmethod
    def get_y0(self, pars) -> np.ndarray:
        pass

    @abstractmethod
    def update_parameters(self, pars):
        pass

    @abstractmethod
    def calc_demography(self, t, y, pars, dy, da, calc):
        pass

    @abstractmethod
    def calc_transmission(self, t, y, pars, dy, da, calc):
        pass

    @abstractmethod
    def calc_progression(self, t, y, pars, dy, da, calc):
        pass

    @abstractmethod
    def calc_cascade(self, t, y, pars, dy, da, calc):
        pass

    @abstractmethod
    def structure_ya(self, y):
        pass

    def __call__(self, t, y, pars):
        y, aux = self.structure_ya(y)

        dy = np.zeros_like(y)
        da = np.zeros_like(aux)

        calc = {'n': y.sum()}
        self.calc_transmission(t, y, pars, dy, da, calc)
        self.calc_progression(t, y, pars, dy, da, calc)
        self.calc_cascade(t, y, pars, dy, da, calc)
        self.calc_demography(t, y, pars, dy, da, calc)

        return np.concatenate([dy.reshape(-1), da])

    @abstractmethod
    def measure(self, t, y, pars):
        pass

    @staticmethod
    def reform_ms(ms):
        ms = pd.DataFrame(ms).set_index('Year')
        ns = ms.N.rolling(2).mean()

        ms1 = pd.DataFrame({
            'N': ns,
            'IncR_apx': ms.IncR.rolling(2).mean().shift(-1),
            'MorR_apx': ms.MorR.rolling(2).mean().shift(-1),
            'CNR_apx': ms.CNR.rolling(2).mean().shift(-1),
            'CumInc': ms.CumInc.rolling(2).mean().shift(-1),
            'CumMor': ms.CumMor.rolling(2).mean().shift(-1),
            'CumACF': ms.CumACF.rolling(2).mean().shift(-1),
            'IncR': (ms.CumInc.diff() / ns).shift(-1),
            'MorR': (ms.CumMor.diff() / ns).shift(-1),
            'DetR': (ms.CumDet.diff() / ns).shift(-1),
            'CNR': (ms.CumNoti.diff() / ns).shift(-1),
            'PrLat': ms.LTBI.rolling(2).mean().shift(-1),
            'PrRecent': (ms.CumIncRecent.diff() / ms.CumInc.diff()).shift(-1),
            'Prev': ms.Prev.rolling(2).mean().shift(-1),
            'PrevA': ms.PrevA.rolling(2).mean().shift(-1),
            'PrevS': ms.PrevS.rolling(2).mean().shift(-1),
            'PrevC': ms.PrevC.rolling(2).mean().shift(-1),
        }).iloc[:-1, ]

        ms1['PrA'] = ms1['PrevA'] / ms1['Prev']
        ms1['PrS'] = ms1['PrevS'] / ms1['Prev']
        ms1['PrC'] = ms1['PrevC'] / ms1['Prev']

        return ms1

    def __sim(self, p, y0, t0, t1, t_eval=None):
        try:
            if t_eval is not None:
                ys = solve_ivp(self, [t0, t1], y0, t_eval=t_eval, args=(p, ))
            else:
                ys = solve_ivp(self, [t0, t1], y0, args=(p, ))
        except ValueError:
            return None, None, {'succ': False}

        if t_eval is not None:
            ms = [self.measure(t_eval[i], ys.y.T[i], p) for i in range(t_eval.shape[0])]
            ms = self.reform_ms(ms)
        else:
            ms = None

        return ys, ms, {'succ': True, 't0': t0, 't1': t1}

    def simulate_to_fit(self, p, t_eval=np.linspace(2014, 2020, 7)):
        p = self.update_parameters(p)
        y0 = self.get_y0(p).reshape(-1)
        return self.__sim(p, y0, self.Year0 - 300, max(t_eval), t_eval)

    def simulate_to_preCOVID(self, p, t_end=2020):
        p = self.update_parameters(p)
        y0 = self.get_y0(p).reshape(-1)
        return self.__sim(p, y0, self.Year0 - 300, t_end)

    def simulate_to_postCOVID(self, p, ys0, t_eval=np.linspace(2020, 2023, 13)):
        p = self.update_parameters(p)
        t0, y0 = ys0.t[-1], ys0.y[:, -1]
        return self.__sim(p, y0, t0, max(t_eval), t_eval)

    def simulate_to_baseline(self, p):
        p = self.update_parameters(p) if 'sus' not in p else p
        y0 = self.get_y0(p).reshape(-1)
        return self.__sim(p, y0, self.YearBaseline - 300, self.YearBaseline)

    def simulate_onward(self, p, ys0, t_end=2036, dt=1, intv=None):
        if intv is not None:
            p['intv'] = intv
        elif 'intv' in p:
            del p['intv']

        t_eval = np.linspace(self.YearBaseline, t_end, int((t_end - self.YearBaseline) / dt) + 1)
        p = self.update_parameters(p)
        t0, t1 = ys0.t[-1], t_end
        y0 = ys0.y[:, -1]
        return self.__sim(p, y0, t0, t1, t_eval)
