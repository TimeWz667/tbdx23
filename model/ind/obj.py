from model.core.inputs import load_inputs
from model.ind.dy import ModelIND
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
    def __init__(self, targets, targets_pupr, targets_prev, bn, model, exo=None):
        AbsObjectiveSimBased.__init__(self, bn, exo)
        self.Model = model
        self.Targets = targets
        self.TargetsPuPr = targets_pupr
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
            for yr in [2017, 2018, 2019]:
                est = ms[index][yr]
                obs = self.Targets[f'{index}_{yr:d}_All']
                dist0 += (est - obs.m) ** 2 / obs.eps

        for index in ['PrevA', 'PrevS', 'PrevC']:
            est = ms[index][2019]
            obs = self.TargetsPrev[f'{index}_2019_All']
            dist0 += (est - obs.m) ** 2 / obs.eps

        for yr in [2017, 2018, 2019]:
            est = ms['CNR_Pub'][yr]
            obs = self.TargetsPuPr[f'CNR_{yr:d}_Pub']
            dist0 += (est - obs.m) ** 2 / obs.eps

            est = ms['CNR_Pri'][yr]
            obs = self.TargetsPuPr[f'CNR_{yr:d}_Pri']
            dist0 += (est - obs.m) ** 2 / obs.eps

        return - dist0

    @staticmethod
    def load(pars_folder, exo=None):
        inp = load_inputs(pars_folder)

        with open(f'{pars_folder}/prior.txt', 'r') as f:
            bn = bayes_net_from_script(f.read())

        m = ModelIND(inp)
        targets = pd.read_csv(f'{pars_folder}/targets.csv')
        targets = targets[targets.Year >= 2017]
        targets = flatten_targets(targets)

        targets_pupr = pd.read_csv(f'{pars_folder}/targets_pupr.csv')
        targets_pupr = targets_pupr[targets_pupr.Year >= 2017]
        targets_pupr = flatten_targets(targets_pupr)
        targets_prev = pd.read_csv(f'{pars_folder}/targets_prev.csv')
        targets_prev = flatten_targets(targets_prev)
        obj = Objective(targets, targets_pupr, targets_prev, bn=bn, model=m, exo=exo)

        return obj


if __name__ == '__main__':
    obj = Objective.load('../../pars/IND')
    p0 = obj.sample_prior()
    s0 = obj.simulate(p0)
    li0 = obj.link_likelihood(s0)
    print(li0)
