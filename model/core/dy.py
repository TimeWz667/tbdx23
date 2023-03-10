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
        rates = self.Inputs.Demography(max(t, self.Year0))
        r_hiv = self.Inputs.DyHIV(t)

        n = calc['n']
        mu = rates['r_die']

        dr_tb = np.zeros_like(y)
        dr_tb[I.Asym] = pars['cas'].R_Death_A
        dr_tb[I.Sym] = pars['cas'].R_Death_S
        dr_tb[I.ExCS] = pars['cas'].R_Death_S
        dr_tb[I.Tx] = pars['r_death_tx']

        dr_hiv = np.zeros_like(y)
        dr_hiv[:, 1] = r_hiv['r_die_hiv']

        mu = mu - ((dr_tb + dr_hiv) * y).sum() / y.sum()

        calc['deaths'] = deaths = mu * y
        calc['deaths_tb'] = deaths_tb = dr_tb * y
        calc['deaths_hiv'] = deaths_hiv = dr_hiv * y

        # Demography
        dy[I.U, 0] += rates['r_birth'] * n
        dy -= deaths + deaths_tb + deaths_hiv

        rh, ra = r_hiv['r_hiv'], r_hiv['r_art']

        dy[:, I.U] += - rh * y[:, I.U]
        dy[:, I.HIV] += rh * y[:, I.U] - ra * y[:, I.HIV]
        dy[:, I.ART] += ra * y[:, I.HIV]
        da[I.A_Mor] += (deaths + deaths_tb)[I.PTB].sum()

    @abstractmethod
    def calc_transmission(self, t, y, pars, dy, da, calc):
        adr = 0 # pars['adr']
        if t < self.Year0:
            adj = 1
        else:
            adj = np.exp(- adr * (t - self.Year0))

        foi = adj * pars['beta'] * (pars['trans'] * y).sum() / calc['n']

        calc['infection'] = infection = foi * pars['sus'] * y

        dy -= infection
        dy[I.FLat] += infection.sum(0)

    @abstractmethod
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

        irr = pars['irr_hiv']
        trs_hiv = [(fr, to, rate * irr if tag.startswith('inc') else rate, tag) for fr, to, rate, tag in trs]

        irr = pars['irr_art']
        trs_art = [(fr, to, rate * irr if tag.startswith('inc') else rate, tag) for fr, to, rate, tag in trs]
        trs = [trs, trs_hiv, trs_art]

        calc['inc_recent'] = np.zeros(I.N_Strata)
        calc['inc_remote'] = np.zeros(I.N_Strata)

        for i in range(I.N_Strata):
            dy[:, i] += calc_dy(y[:, i], trs[i])

            calc['inc_recent'][i] = extract_tr(y[:, i], trs[i], fil=lambda x: x[3] == 'inc_recent')
            calc['inc_remote'][i] = extract_tr(y[:, i], trs[i], fil=lambda x: x[3] == 'inc_remote')

        calc['inc'] = calc['inc_recent'] + calc['inc_remote']

        da[I.A_Inc] += calc['inc'].sum()
        da[I.A_IncRecent] += calc['inc_recent'].sum()
        da[I.A_IncRemote] += calc['inc_remote'].sum()
        da[I.A_Inc_NonHIV] += calc['inc'][0]
        da[I.A_Inc_PLHIV] += calc['inc'][1:].sum()

    @abstractmethod
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
            r_csi_acf = blend(0, intv['r_csi_acf'], wt) * pdx0
            r_recsi_acf = blend(0, intv['r_recsi_acf'], wt) * pdx1
            r_asym_acf = blend(0, intv['r_asym_acf'], wt)
        else:
            r_csi_acf = 0
            r_recsi_acf = 0
            r_asym_acf = 0

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

        trs_hiv = [list(trs) for _ in range(3)]
        trs_hiv[1].append((I.Asym, I.Tx, r_asym_acf, 'acf_a'))
        trs_hiv[2].append((I.Asym, I.Tx, r_asym_acf, 'acf_a'))

        calc['det'] = 0
        calc['acf'] = 0
        calc['acf_a'] = 0

        for i in range(I.N_Strata):
            dy[:, i] += calc_dy(y[:, i], trs_hiv[i])
            calc['det'] += extract_tr(y[:, i], trs_hiv[i], fil=lambda x: x[3] == 'det')
            calc['acf'] += extract_tr(y[:, i], trs_hiv[i], fil=lambda x: x[3] == 'acf')
            calc['acf_a'] += extract_tr(y[:, i], trs_hiv[i], fil=lambda x: x[3] == 'acf_a')

        da[I.A_Det] += calc['det']
        da[I.A_ACF] += calc['acf']
        da[I.A_Noti] += calc['det'] * cas.PrReport(t) / cas.PPV

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
