from sim.core.inputs import load_inputs
from sim.core.cascade import RepoCascade
from sim.core.dy import Model
from sim.zaf.dhiv import ModelHIV
from sims_pars import bayes_net_from_script
from sims_pars.fitting import AbsObjectiveSimBased
import numpy as np
import pandas as pd
import json
import time

__author__ = 'Chu-Chang Ku'
__all__ = ['Objective', 'load_objective']


class Objective(AbsObjectiveSimBased):
    def __init__(self, targets, bn, model, cs: RepoCascade, exo=None):
        AbsObjectiveSimBased.__init__(self, bn, exo)
        self.Model = model
        self.RepoCascade = cs
        self.Targets = targets

    def simulate(self, pars):
        time.sleep(0.001)
        if 'cas' not in pars:
            pars['cas'] = self.RepoCascade.sample()
        return self.Model.simulate_to_fit(pars)

    def link_likelihood(self, sim):
        _, ms, msg = sim
        if not msg['succ']:
            return - np.Inf

        tar = self.Targets

        dist = ((ms.IncR - tar.IncR_mu) ** 2 / tar.IncR_eps).sum()
        dist += ((ms.MorR - tar.MorR_mu) ** 2 / tar.MorR_eps).sum()
        dist += ((ms.CNR - tar.CNR_mu) ** 2 / tar.CNR_eps).sum()
        # dist = ((ms.IncR / tar.IncR_mu - 1) ** 2).sum()
        # dist += ((ms.MorR / tar.MorR_mu - 1) ** 2).sum()
        # dist += ((ms.CNR / tar.CNR_mu - 1) ** 2).sum()
        return - dist


class ObjectiveHIV(Objective):
    def link_likelihood(self, sim):
        _, ms, msg = sim
        if not msg['succ']:
            return - np.Inf

        tar = self.Targets['All']

        dist = ((ms.IncR - tar.IncR_mu) ** 2 / tar.IncR_eps).sum()
        dist += ((ms.MorR - tar.MorR_mu) ** 2 / tar.MorR_eps).sum()
        dist += ((ms.CNR - tar.CNR_mu) ** 2 / tar.CNR_eps).sum()

        tar = self.Targets['HIV']
        tar = tar.loc[tar.Group == 'HIV']
        dist += ((ms.IncR_PLHIV - tar.IncR_mu) ** 2 / tar.IncR_eps).sum()
        # dist += ((ms.MorR - tar.MorR_mu) ** 2 / tar.MorR_eps).sum()
        # dist += ((ms.CNR - tar.CNR_mu) ** 2 / tar.CNR_eps).sum()

        tar = self.Targets['HIV']
        tar = tar.loc[tar.Group == 'NonHIV']
        dist += ((ms.IncR_NonHIV - tar.IncR_mu) ** 2 / tar.IncR_eps).sum()

        return - dist


def load_objective(folder_inputs, file_cs, exo=None, hiv=False):
    pars_cs = pd.read_csv(file_cs)
    repo_cs = RepoCascade(pars_cs)

    inp = load_inputs(folder_inputs)

    with open(f'{folder_inputs}/prior.txt', 'r') as f:
        bn = bayes_net_from_script(f.read())

    if not hiv:
        m = Model(inp)
        targets = pd.read_csv(f'{folder_inputs}/targets.csv').set_index('Year')
        obj = Objective(targets, bn=bn, model=m, cs=repo_cs, exo=exo)
    else:
        m = ModelHIV(inp)
        targets = {
            'All': pd.read_csv(f'{folder_inputs}/targets.csv').set_index('Year'),
            'HIV': pd.read_csv(f'{folder_inputs}/targets_hiv.csv').set_index('Year')
        }
        obj = ObjectiveHIV(targets, bn=bn, model=m, cs=repo_cs, exo=exo)

    return obj


if __name__ == '__main__':
    obj = load_objective('../../data/pars/ZAF', '../../results/pars_ZAF.csv', hiv=False)

    for do in obj.Domain:
        print(do)

    p0 = obj.sample_prior()

    print(obj.calc_likelihood(p0))

    # With HIV
    obj = load_objective('../../data/pars/ZAF', '../../results/pars_ZAF.csv', hiv=True)

    for do in obj.Domain:
        print(do)

    p0 = obj.sample_prior()

    print(obj.calc_likelihood(p0))
