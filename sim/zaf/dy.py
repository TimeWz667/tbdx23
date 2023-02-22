from sim.core.dy import Model

__author__ = 'Chu-Chang Ku'
__all__ = ['ModelZAF']


def simulate_intv(self, p, ys0, intv=None):
    if intv is not None:
        p['intv'] = intv
    elif 'intv' in p:
        del p['intv']

    ys, ms, msg = self.simulate_onward(p, ys0, 2036, 1)
    return ys, ms, msg


class ModelZAF(Model):
    def update_parameters(self, pars):
        pars = Model.update_parameters(self, pars)
        pars['cas'].reform(pdx0=0.45, pdx1=0.45)
        return pars
