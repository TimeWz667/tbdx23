import numpy as np
import pandas as pd
from model.core.dy import AbsModel
from model.util import calc_dy, extract_tr, blend
import model.ind.keys as I

__author__ = 'Chu-Chang Ku'
__all__ = ['ModelIND']


class ModelIND(AbsModel):
    def update_parameters(self, pars):
        pars = dict(pars)

        pars['sus'] = sus = np.zeros(I.N_State_TB)

        sus[I.U] = 1
        sus[I.SLat] = pars['rr_sus_slat']
        sus[I.RLow] = pars['rr_sus_rec']
        sus[I.RHigh] = pars['rr_sus_rec']
        sus[I.RSt] = pars['rr_sus_rec']

        pars['trans'] = trans = np.zeros(I.N_State_TB)
        trans[I.Asym] = pars['rr_inf_asym']
        trans[I.Sym] = 1
        trans[I.ExCS] = pars['rr_inf_cs']

        pars['r_ts'] = r_ts = 1 / pars['dur_succ']

        to = self.Inputs.TxOut['Pub']
        pars['r_tl_pub'] = r_ts * to.PrTxLTFU / to.PrTxSucc
        pars['r_death_tx_pub'] = r_ts * to.PrTxDie / to.PrTxSucc

        to = self.Inputs.TxOut['Pri']
        pars['r_tl_pri'] = r_ts * to.PrTxLTFU / to.PrTxSucc
        pars['r_death_tx_pri'] = r_ts * to.PrTxDie / to.PrTxSucc

        if 'k_covid' not in pars:
            pars['k_covid'] = lambda t: 1

        return pars

    def get_y0(self, pars) -> np.ndarray:
        y0 = np.zeros(I.N_State_TB)

        n0 = self.Inputs.Demography.N0

        y0[I.Asym] = n0 * 0.001
        y0[I.Sym] = n0 * 0.001
        y0[I.ExCS] = n0 * 0.001
        y0[I.SLat] = n0 * 0.5
        y0[I.U] = n0 - y0.sum()

        return np.concatenate([y0.reshape(-1), np.zeros(I.N_Aux)])

    def calc_demography(self, t, y, pars, dy, da, calc):
        # rates = self.Inputs.Demography(max(t, self.Year0))
        rates = self.Inputs.Demography(max(t, self.Year0))
        n = calc['n']
        mu = rates['r_die']

        dr_tb = np.zeros_like(y)
        dr_tb[I.Asym] = pars['r_die_ut'] * pars['rr_die_a']
        dr_tb[I.Sym] = pars['r_die_ut']
        dr_tb[I.ExCS] = pars['r_die_ut']

        p_cs_pub, p_dx_pub, p_dx_pri = pars['p_cs_pub'], pars['p_dx_pub'], pars['p_dx_pri']

        pdx = p_cs_pub * p_dx_pub + (1 - p_cs_pub) * p_dx_pri

        r_death_tx = pars['r_death_tx_pub'] * p_cs_pub * p_dx_pub + pars['r_death_tx_pri'] * (1 - p_cs_pub) * p_dx_pri
        r_death_tx /= pdx

        dr_tb[I.Tx] = r_death_tx

        mu = mu - (dr_tb * y).sum() / y.sum()

        calc['deaths'] = deaths = mu * y
        calc['deaths_tb'] = deaths_tb = dr_tb * y

        # Demography
        dy[I.U] += rates['r_birth'] * n
        # dy[I.U] += (deaths + deaths_tb).sum()
        dy -= deaths + deaths_tb

        da[I.A_Mor] += (deaths + deaths_tb)[I.PTB].sum()

        if t <= self.Year0:
            dy -= y / y.sum() * dy.sum()

    def calc_transmission(self, t, y, pars, dy, da, calc):
        foi = pars['beta'] * (pars['trans'] * y).sum() / calc['n']

        calc['infection'] = infection = foi * pars['sus'] * y

        dy -= infection
        dy[I.FLat] += infection.sum()

    def calc_progression(self, t, y, pars, dy, da, calc):
        r_lat, r_stab = pars['r_lat'], pars['r_stab']
        r_act, r_react = pars['r_act'], pars['r_react']
        r_rel, r_rel_td, r_rel_tc = pars['r_relapse'], pars['r_relapse_td'], pars['r_relapse_tc']
        r_sc, r_clear = pars['r_sc'], pars['r_clear']
        r_onset = pars['r_onset']

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

        dy += calc_dy(y, trs)

        calc['inc_recent'] = extract_tr(y, trs, fil=lambda x: x[3] == 'inc_recent')
        calc['inc_remote'] = extract_tr(y, trs, fil=lambda x: x[3] == 'inc_remote')

        calc['inc'] = calc['inc_recent'] + calc['inc_remote']

        da[I.A_Inc] += calc['inc']
        da[I.A_IncRecent] += calc['inc_recent']
        da[I.A_IncRemote] += calc['inc_remote']

    def calc_cascade(self, t, y, pars, dy, da, calc):
        try:
            k_covid = pars['k_covid'](t)
        except KeyError:
            k_covid = 1

        p_cs_pub, p_dx_pub, p_dx_pri = pars['p_cs_pub'], pars['p_dx_pub'], pars['p_dx_pri']

        pdx0 = pdx1 = p_cs_pub * p_dx_pub + (1 - p_cs_pub) * p_dx_pri

        k_cs = np.exp(pars['rt_cs'] * (max(t, 2010) - 2021)) if t < 2021 else 1
        # k_cs = 1
        r_csi, r_recsi = pars['r_csi'], pars['r_recsi']

        if t > 2023 and 'intv' in pars:
            wt = (t - 2023) / (2026 - 2023) if t < 2026 else 1
            intv = pars['intv']
            pdx0 = blend(pdx0, intv['pdx0'], wt)
            pdx1 = blend(pdx1, intv['pdx1'], wt)
            r_csi_acf = blend(0, intv['r_csi_acf'], wt) * pdx0
            r_recsi_acf = blend(0, intv['r_recsi_acf'], wt) * pdx1
            r_asym_acf = blend(0, intv['r_asym_acf'], wt)
        else:
            r_csi_acf = 0
            r_recsi_acf = 0
            r_asym_acf = 0

        r_det_s = k_cs * k_covid * r_csi * pdx0
        r_fn_s = k_cs * k_covid * r_csi * (1 - pdx0)
        r_det_c = k_cs * k_covid * r_recsi * pdx1

        r_tl = pars['r_tl_pub'] * p_cs_pub * p_dx_pub + pars['r_tl_pri'] * (1 - p_cs_pub) * p_dx_pri

        trs = [
            (I.Asym, I.Tx, r_asym_acf, 'acf_a'),
            (I.Sym, I.Tx, r_det_s, 'det'),
            (I.Sym, I.ExCS, r_fn_s, 'fn'),
            (I.ExCS, I.Tx, r_det_c, 'det'),
            (I.Sym, I.Tx, r_csi_acf, 'acf'),
            (I.ExCS, I.Tx, r_recsi_acf, 'acf'),
            (I.Tx, I.RLow, pars['r_ts'], 'txs'),
            (I.Tx, I.RHigh, r_tl, 'txl'),
        ]

        dy += calc_dy(y, trs)
        calc['det'] = extract_tr(y, trs, fil=lambda x: x[3] == 'det')
        calc['acf'] = extract_tr(y, trs, fil=lambda x: x[3] == 'acf')
        calc['acf_a'] = extract_tr(y, trs, fil=lambda x: x[3] == 'acf_a')

        da[I.A_Det] += calc['det']
        da[I.A_ACF] += calc['acf']
        da[I.A_Noti] += calc['det'] / pars['ppv']

    def structure_ya(self, y):
        y, aux = y[:- I.N_Aux], y[- I.N_Aux:]
        y = y.reshape(I.N_State_TB)
        return y, aux

    def measure(self, t, y, pars):
        y, aux = self.structure_ya(y)

        dy = np.zeros_like(y)
        da = np.zeros_like(aux)

        n = y.sum()
        calc = {'n': n}
        self.calc_transmission(t, y, pars, dy, da, calc)
        self.calc_progression(t, y, pars, dy, da, calc)
        self.calc_cascade(t, y, pars, dy, da, calc)
        self.calc_demography(t, y, pars, dy, da, calc)

        mor = (calc['deaths_tb'] + calc['deaths'])[I.PTB].sum()

        mea = {
            'Year': t,
            'N': n,
            'Prev': y[I.PTB].sum() / n,
            'PrevA': y[I.Asym].sum() / n,
            'PrevS': y[I.Sym].sum() / n,
            'PrevC': y[I.ExCS].sum() / n,
            'IncR': calc['inc'] / n,
            'MorR': mor / n,
            'LTBI': y[I.LTBI].sum() / n,
            'DetR': calc['det'] / n,
            'CNR': calc['det'] / n * (1 / pars['ppv']),
            'CumInc': aux[I.A_Inc],
            'CumIncRecent': aux[I.A_IncRecent],
            'CumIncRemote': aux[I.A_IncRemote],
            'CumMor': aux[I.A_Mor],
            'CumDet': aux[I.A_Det],
            'CumNoti': aux[I.A_Noti],
            'CumNotiAcf': aux[I.A_CNR_Acf],
            'CumNotiPub': aux[I.A_CNR_Pub],
            'CumNotiPri': aux[I.A_CNR_Pri],
            'CumACF': aux[I.A_ACF]
        }
        return mea

    @staticmethod
    def reform_ms(ms):
        ms1 = AbsModel.reform_ms(ms)

        ms = pd.DataFrame(ms).set_index('Year')
        ns = ms.N.rolling(2).mean()

        ms2 = pd.DataFrame({
            'CNR_Acf': (ms.CumNotiAcf.diff() / ns).shift(-1),
            'CNR_Pub': (ms.CumNotiPub.diff() / ns).shift(-1),
            'CNR_Pri': (ms.CumNotiPri.diff() / ns).shift(-1)
        }).iloc[:-1, ]

        return ms1.join(ms2)


if __name__ == '__main__':
    from model.core.inputs import load_inputs
    from sims_pars import bayes_net_from_script, sample
    import matplotlib.pyplot as plt
    import numpy.random as rd

    rd.seed(1167)

    with open('../../pars/IND/prior.txt', 'r') as f:
        r_prior = bayes_net_from_script(f.read())

    pars = dict(sample(r_prior))
    pars['beta'] = 11
    pars['rr_inf_asym'] = 1
    pars['r_recsi'] = 3
    pars['r_csi'] = 3
    pars['rt_cs'] = 0.01

    inp = load_inputs('../../pars/IND')
    model = ModelIND(inp)
    pars = model.update_parameters(pars)

    y0 = model.get_y0(pars)

    ys, ms, _ = model.simulate_to_fit(pars, t_eval=np.linspace(2000, 2030, 31))

    fig, axes = plt.subplots(2, 3)

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
    ms.MorR_apx.plot(ax=axes[1, 1])
    axes[1, 1].set_title('Mortality')

    ms.PrA.plot(ax=axes[0, 2])
    ms.PrS.plot(ax=axes[0, 2])
    ms.PrC.plot(ax=axes[0, 2])
    axes[1, 2].set_title('PrPrev')

    ms.PrLat.plot(ax=axes[1, 2])
    ms.PrRecent.plot(ax=axes[1, 2])
    axes[1, 2].set_title('Prop')

    fig.tight_layout()
    plt.show()
