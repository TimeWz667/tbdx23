import numpy as np
import sim.core.keys as I
from sim.util import calc_dy, extract_tr
from scipy.integrate import solve_ivp
import pandas as pd


__author__ = 'Chu-Chang Ku'
__all__ = ['Model']


def blend(x0, x1, wt):
    return x0 + wt * (x1 - x0)


class Model:
    def __init__(self, inputs):
        self.Inputs = inputs
        self.Year0 = inputs.Demography.Year0
        self.YearBaseline = 2022

    def get_y0(self, pars) -> np.ndarray:
        y0 = np.zeros((I.N_State_TB, I.N_Strata))

        n0 = self.Inputs.Demography.N0

        y0[I.Asym, 0] = n0 * pars['cas'].Prev['Asym']
        y0[I.Sym, 0] = n0 * pars['cas'].Prev['Sym']
        y0[I.ExCS, 0] = n0 * pars['cas'].Prev['ExCS']
        y0[I.SLat, 0] = n0 * 0.5
        y0[I.U, 0] = n0 - y0.sum(0, keepdims=True)

        return np.concatenate([y0.reshape(-1), np.zeros(I.N_Aux)])

    def update_parameters(self, pars):
        pars = dict(pars)

        pars['sus'] = sus = np.zeros((I.N_State_TB, I.N_Strata))

        sus[I.U] = 1
        sus[I.SLat] = pars['rr_sus_slat']
        sus[I.RLow] = pars['rr_sus_rec']
        sus[I.RHigh] = pars['rr_sus_rec']
        sus[I.RSt] = pars['rr_sus_rec']

        pars['trans'] = trans = np.zeros((I.N_State_TB, I.N_Strata))
        trans[I.Asym] = pars['rr_inf_asym']
        trans[I.Sym] = 1
        trans[I.ExCS] = pars['rr_inf_cs']

        to = self.Inputs.TxOut
        pars['r_ts'] = r_ts = 1 / pars['dur_succ']
        pars['r_tl'] = r_ts * to.PrTxLTFU / to.PrTxSucc
        pars['r_death_tx'] = r_ts * to.PrTxDie / to.PrTxSucc

        if 'k_covid' not in pars:
            pars['k_covid'] = lambda t: 1

        return pars

    def calc_demography(self, t, y, pars, dy, da, calc):
        rates = self.Inputs.Demography(max(t, self.Year0))

        n = calc['n']

        mu = rates['r_die']

        dr_tb = np.zeros_like(y)
        dr_tb[I.Asym] = pars['cas'].R_Death_A
        dr_tb[I.Sym] = pars['cas'].R_Death_S
        dr_tb[I.ExCS] = pars['cas'].R_Death_S
        dr_tb[I.Tx] = pars['r_death_tx']

        mu = mu - (dr_tb * y).sum() / y.sum()

        calc['deaths'] = deaths = mu * y
        calc['deaths_tb'] = deaths_tb = dr_tb * y

        # Demography
        dy[I.U, 0] += rates['r_birth'] * n
        dy -= deaths + deaths_tb

        da[I.A_Mor] += (deaths + deaths_tb)[I.PTB].sum()

    def calc_transmission(self, t, y, pars, dy, da, calc):
        adr = pars['adr']
        if t < self.Year0:
            adj = 1
        else:
            adj = np.exp(- adr * (t - self.Year0))
        # elif t < 2026:
        #     adj = np.exp(- adr * (2023 - self.Year0))
        #     wt = 0.5 * blend(1, 0.5, (t - 2023) / (2026 - 2023))
        #     adj *= np.exp(- adr * wt)
        # else:
        #     adj = np.exp(- adr * (2023 - self.Year0))
        #     adj *= np.exp(- adr * np.exp(- adr * (2026 - 2023)))
        #     adj *= np.exp(- adr / 2 * (t - 2026))

        foi = adj * pars['beta'] * (pars['trans'] * y).sum() / calc['n']

        calc['infection'] = infection = foi * pars['sus'] * y

        dy -= infection
        dy[I.FLat] += infection.sum(0)

    def calc_progression(self, t, y, pars, dy, da, calc):
        r_lat, r_stab = pars['r_lat'], pars['r_stab']
        r_act, r_react = pars['r_act'], pars['r_react']
        r_rel, r_rel_td, r_rel_tc = pars['r_relapse'], pars['r_relapse_td'], pars['r_relapse_tc']
        r_sc, r_clear = pars['r_sc'], pars['r_clear']
        r_onset = pars['cas'].R_Onset

        trs = [
            (I.FLat, I.Asym, r_act, 'inc_recent'),
            (I.SLat, I.Asym, r_react, 'inc_remote'),
            (I.RLow, I.Asym, r_rel_tc, 'inc_remote'),
            (I.RHigh, I.Asym, r_rel_td, 'inc_remote'),
            (I.RSt, I.Asym, r_rel, 'inc_remote'),

            (I.FLat, I.SLat, r_lat, 'stab'),
            (I.RLow, I.RSt, r_stab, 'stab'),
            (I.RHigh, I.RSt, r_stab, 'stab'),

            (I.Asym, I.Sym, r_onset, 'onset'),

            (I.Asym, I.SLat, r_sc, 'self_cure'),
            (I.Sym, I.SLat, r_sc, 'self_cure'),
            (I.ExCS, I.SLat, r_sc, 'self_cure'),

            (I.SLat, I.U, r_clear, 'self_clear'),
            (I.RSt, I.U, r_clear, 'self_clear'),
        ]

        calc['inc_recent'] = 0
        calc['inc_remote'] = 0

        for i in range(I.N_Strata):
            dy[:, i] += calc_dy(y[:, i], trs)
            calc['inc_recent'] += extract_tr(y[:, i], trs, fil=lambda x: x[3] == 'inc_recent')
            calc['inc_remote'] += extract_tr(y[:, i], trs, fil=lambda x: x[3] == 'inc_remote')

        calc['inc'] = calc['inc_recent'] + calc['inc_remote']

        da[I.A_Inc] += calc['inc']
        da[I.A_IncRecent] += calc['inc_recent']
        da[I.A_IncRemote] += calc['inc_remote']

    def calc_cascade(self, t, y, pars, dy, da, calc):
        cas = pars['cas']

        try:
            k_covid = pars['k_covid'](t)
        except KeyError:
            k_covid = 1

        pdx0, pdx1 = cas.PrDx0, cas.PrDx1
        r_csi, r_recsi = cas.R_CSI, cas.R_ReCSI

        if t > 2023 and 'intv' in pars:
            wt = (t - 2023) / (2026 - 2023) if t < 2026 else 1
            intv = pars['intv']
            pdx0 = blend(pdx0, intv['pdx0'], wt)
            pdx1 = blend(pdx1, intv['pdx1'], wt)
            r_csi_acf = blend(0, intv['r_csi_acf'], wt) * pars['p_dx_pub']
            r_recsi_acf = blend(0, intv['r_recsi_acf'], wt) * pars['p_dx_pub']
        else:
            r_csi_acf = 0
            r_recsi_acf = 0

        r_det_s = k_covid * r_csi * pdx0
        r_fn_s = k_covid * r_csi * (1 - pdx0)
        r_det_c = k_covid * r_recsi * pdx1

        trs = [
            (I.Sym, I.Tx, r_det_s, 'det'),
            (I.Sym, I.ExCS, r_fn_s, 'fn'),
            (I.ExCS, I.Tx, r_det_c, 'det'),
            (I.Sym, I.Tx, r_csi_acf, 'acf'),
            (I.ExCS, I.Tx, r_recsi_acf, 'acf'),
            (I.Tx, I.RLow, pars['r_ts'], 'txs'),
            (I.Tx, I.RHigh, pars['r_tl'], 'txl'),
        ]

        calc['det'] = 0
        calc['acf'] = 0
        for i in range(I.N_Strata):
            dy[:, i] += calc_dy(y[:, i], trs)
            calc['det'] += extract_tr(y[:, i], trs, fil=lambda x: x[3] == 'det')
            calc['acf'] += extract_tr(y[:, i], trs, fil=lambda x: x[3] == 'acf')

        da[I.A_Det] += calc['det']
        da[I.A_ACF] += calc['acf']
        da[I.A_Noti] += calc['det'] * cas.PrReport(t) / cas.PPV

    def __call__(self, t, y, pars):
        y, aux = y[:- I.N_Aux], y[- I.N_Aux:]
        y = y.reshape((I.N_State_TB, I.N_Strata))

        dy = np.zeros_like(y)
        da = np.zeros_like(aux)

        calc = {'n': y.sum()}
        self.calc_transmission(t, y, pars, dy, da, calc)
        self.calc_progression(t, y, pars, dy, da, calc)
        self.calc_cascade(t, y, pars, dy, da, calc)
        self.calc_demography(t, y, pars, dy, da, calc)

        if t <= self.Year0:
            dy -= y / y.sum(0, keepdims=True) * dy.sum(0, keepdims=True)

        return np.concatenate([dy.reshape(-1), da])

    def measure(self, t, y, pars):
        y, aux = y[:- I.N_Aux], y[- I.N_Aux:]
        y = y.reshape((I.N_State_TB, I.N_Strata))

        dy = np.zeros_like(y)
        da = np.zeros_like(aux)

        n = y.sum()
        calc = {'n': n}
        self.calc_transmission(t, y, pars, dy, da, calc)
        self.calc_progression(t, y, pars, dy, da, calc)
        self.calc_cascade(t, y, pars, dy, da, calc)
        self.calc_demography(t, y, pars, dy, da, calc)

        mor = (calc['deaths_tb'] + calc['deaths'])[I.PTB].sum()

        cas = pars['cas']
        mea = {
            'Time': t,
            'N': n,
            'Prev': y[I.PTB].sum() / n,
            'PrevA': y[I.Asym].sum() / n,
            'PrevS': y[I.Sym].sum() / n,
            'PrevC': y[I.ExCS].sum() / n,
            'IncR': calc['inc'] / n,
            'MorR': mor / n,
            'LTBI': y[I.LTBI].sum() / n,
            'DetR': calc['det'] / n,
            'CNR': calc['det'] / n * (cas.PrReport(t) / cas.PPV),
            'CumInc': aux[I.A_Inc],
            'CumIncRecent': aux[I.A_IncRecent],
            'CumIncRemote': aux[I.A_IncRemote],
            'CumMor': aux[I.A_Mor],
            'CumDet': aux[I.A_Det],
            'CumNoti': aux[I.A_Noti],
            'CumACF': aux[I.A_ACF]
        }
        return mea

    def dfe(t, y, pars):
        return y[I.Sym].min()

    dfe.terminal = True
    dfe.direction = -1

    def __sim(self, p, y0, t0, t1, t_eval=None):
        if t_eval is not None:
            ys = solve_ivp(self, [t0, t1], y0, t_eval=t_eval, args=(p, ))
        else:
            ys = solve_ivp(self, [t0, t1], y0, args=(p, ))

        # if len(ys0.t_events[0]) > 0 or not ys0.success:
        #     return None, None, {'succ': False, 'res': 'DFE reached'}
        if t_eval is not None:
            ms = [self.measure(t_eval[i], ys.y.T[i], p) for i in range(t_eval.shape[0])]
            ms = self.reform_ms(ms)
        else:
            ms = None

        return ys, ms

    @staticmethod
    def reform_ms(ms):
        ms = pd.DataFrame(ms).set_index('Time')
        ns = ms.N.rolling(2).mean()
        ms1 = pd.DataFrame({
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
        return ms1

    def simulate_to_fit(self, p, t_eval=np.linspace(2014, 2020, 7)):
        p = self.update_parameters(p) if 'sus' not in p else p
        y0 = self.get_y0(p)

        try:
            ys, ms = self.__sim(p, y0, self.Year0 - 300, max(t_eval), t_eval)
        except ValueError:
            return None, None, {'succ': False}

        return ys, ms, {'succ': True}

    def simulate_to_preCOVID(self, p, t_end=2020):
        p = self.update_parameters(p) if 'sus' not in p else p
        y0 = self.get_y0(p)

        try:
            ys, ms = self.__sim(p, y0, self.Year0 - 300, t_end)
        except ValueError:
            return None, None, {'succ': False}
        return ys, ms, {'succ': True}

    def simulate_to_postCOVID(self, p, ys0, t_eval=np.linspace(2020, 2023, 13)):
        p = self.update_parameters(p) if 'sus' not in p else p
        t0, y0 = ys0.t[-1], ys0.y[:, -1]

        try:
            ys, ms = self.__sim(p, y0, t0, max(t_eval), t_eval)
        except ValueError:
            return None, None, {'succ': False}
        return ys, ms, {'succ': True}

    def simulate_to_baseline(self, p):
        p = self.update_parameters(p) if 'sus' not in p else p

        t0, t1 = self.Year0 - 300, self.YearBaseline
        y0 = self.get_y0(p).reshape(-1)

        try:
            ys, ms = self.__sim(p, y0, self.YearBaseline - 300, self.YearBaseline)
        except ValueError:
            return None, {'succ': False}

        return ys, {'succ': True, 't0': t0, 't1': t1}

    def simulate_onward(self, p, ys0, t_end=2036, dt=1):
        t_eval = np.linspace(self.YearBaseline, t_end, int((t_end - self.YearBaseline) / dt) + 1)
        p = self.update_parameters(p) if 'sus' not in p else p
        t0, y0, t1 = ys0.t[-1], ys0.y[:, -1], max(t_eval)

        try:
            ys, ms = self.__sim(p, y0, t0, t1, t_eval)
        except ValueError:
            return None, None, {'succ': False}

        return ys, ms, {'succ': True, 't0': t0}


if __name__ == '__main__':
    from sim.core.inputs import load_inputs
    from sims_pars import bayes_net_from_script, sample
    from sim.core.cascade import RepoCascade
    import matplotlib.pyplot as plt
    import numpy.random as rd

    rd.seed(1167)

    with open('../../data/pars/IND/prior.txt', 'r') as f:
        r_prior = bayes_net_from_script(f.read())

    pars_cs = pd.read_csv('../../results/pars_IND.csv')
    repo_cs = RepoCascade(pars_cs, 2010)

    pars = dict(sample(r_prior))
    pars['beta'] = 6.9
    pars['rr_inf_asym'] = 1
    pars['adr'] = 0.005
    pars['cas'] = repo_cs.sample()

    inp = load_inputs('../../data/pars/IND')
    model = Model(inp)
    pars = model.update_parameters(pars)

    y0 = model.get_y0(pars)

    ys, ms, _ = model.simulate_to_fit(pars, t_eval=np.linspace(2000, 2030, 31))

    fig, axes = plt.subplots(2, 2)

    ms.Prev.plot(ax=axes[0, 0])
    ms.PrevA.plot(ax=axes[0, 0])
    ms.PrevS.plot(ax=axes[0, 0])
    ms.PrevC.plot(ax=axes[0, 0])
    axes[0, 0].set_title('Prevalence')
    ms.CNR_apx.plot(ax=axes[0, 1])
    ms.DetR.plot(ax=axes[0, 1])
    axes[0, 1].set_title('CNR')
    ms.IncR.plot(ax=axes[1, 0])
    axes[1, 0].set_title('Incidence')
    ms.MorR.plot(ax=axes[1, 1])
    axes[1, 1].set_title('Mortality')

    fig.tight_layout()
    plt.show()
