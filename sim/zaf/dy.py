from sim.core.dy import Model

__author__ = 'Chu-Chang Ku'
__all__ = ['ModelZAF']


class ModelZAF(Model):
    def update_parameters(self, pars):
        pars = Model.update_parameters(self, pars)
        pars['cas'].reform(pdx0=0.4, pdx1=0.7)
        return pars
