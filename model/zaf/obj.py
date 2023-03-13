from model.core.inputs import load_inputs
from model.zaf.dy import ModelZAF
from sims_pars import bayes_net_from_script
from sims_pars.fitting import AbsObjectiveSimBased
import numpy as np
import pandas as pd
import time

__author__ = 'Chu-Chang Ku'
__all__ = ['Objective']


def flatten_targets(tars):
    res = dict()
    for _, tar in tars.iterrows():
        res[f'{tar["Index"]}_{tar["Year"]}_{tar["Tag"]}'] = tar
    return res


class Objective(AbsObjectiveSimBased):
    def __init__(self, targets, targets_hiv, targets_prev, bn, model, exo=None):
        AbsObjectiveSimBased.__init__(self, bn, exo)
        self.Model = model
        self.Targets = targets
        self.TargetsHIV = targets_hiv
        self.TargetsPrev = targets_prev

    def simulate(self, pars):
        time.sleep(0.001)
        return self.Model.simulate_to_fit(pars)

    def link_likelihood(self, sim):
        _, ms, msg = sim
        if not msg['succ']:
            return - np.Inf

        dist0 = 0

        for index in ['MorR', 'IncR', 'CNR']:
            for yr in range(2015, 2020):
                est = ms[index][yr]
                obs = self.Targets[f'{index}_{yr:d}_All']
                dist0 += (est - obs.m) ** 2 / obs.eps

        for yr in range(2015, 2020):
            for tag in ['NonHIV', 'PLHIV']:
                est = ms[f'IncR_{tag}'][yr]
                obs = self.TargetsHIV[f'IncR_{yr:d}_{tag}']
                dist0 += (est - obs.m) ** 2 / obs.eps

        yr = 2018
        for index in ['PrevA', 'PrevS', 'PrevC']:
            est = ms[index][yr]
            obs = self.TargetsPrev[f'{index}_{yr}_All']
            dist0 += (est - obs.m) ** 2 / obs.eps

        return - dist0

    @staticmethod
    def load(pars_folder, exo=None):
        inp = load_inputs(pars_folder)

        with open(f'{pars_folder}/prior.txt', 'r') as f:
            bn = bayes_net_from_script(f.read())

        m = ModelZAF(inp)
        targets = pd.read_csv(f'{pars_folder}/targets.csv')
        targets = flatten_targets(targets)

        targets_hiv = pd.read_csv(f'{pars_folder}/targets_hiv.csv')
        targets_hiv = flatten_targets(targets_hiv)
        targets_prev = pd.read_csv(f'{pars_folder}/targets_prev.csv')
        targets_prev = flatten_targets(targets_prev)
        obj = Objective(targets, targets_hiv, targets_prev, bn=bn, model=m, exo=exo)

        return obj


if __name__ == '__main__':
    obj = Objective.load('../../pars/ZAF')
    p0 = obj.sample_prior()
    s0 = obj.simulate(p0)
    li0 = obj.link_likelihood(s0)
    print(li0)
