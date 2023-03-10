from sim.core.dy import Model

__author__ = 'Chu-Chang Ku'
__all__ = ['ModelIND', 'get_intv']


def get_intv(p, pdx0_pub=None, pdx0_eng=None, pdx1_pub=None, pdx1_eng=None, ppm=None, rr_csi=1, rr_recsi=1, rd_csi=0, rd_recsi=0, r_asym_acf=0):
    cas = p['cas']

    pdx_pub = p['p_dx_pub']
    entry_pub = p['p_cs_pub']
    pdx_all = cas.PrDx0
    pdx_pri = (pdx_all - pdx_pub * entry_pub) / (1 - entry_pub)
    ppm0 = 1 - (1 - cas.PrReport(2024)) * pdx_all / (1 - entry_pub) / pdx_pri
    ppm0 = max(ppm0, 0)
    ppm = max(ppm, ppm0) if ppm is not None else ppm0

    pdx0_pub = pdx_pub if pdx0_pub is None else pdx0_pub
    pdx0_eng = pdx_pri if pdx0_eng is None else pdx0_eng
    pdx0 = entry_pub * pdx0_pub + (1 - entry_pub) * ppm * pdx0_eng + (1 - entry_pub) * (1 - ppm) * pdx_pri

    pdx1_pub = pdx_pub if pdx1_pub is None else pdx1_pub
    pdx1_eng = pdx_pri if pdx1_eng is None else pdx1_eng
    pdx1 = entry_pub * pdx1_pub + (1 - entry_pub) * ppm * pdx1_eng + (1 - entry_pub) * (1 - ppm) * pdx_pri

    r_csi1 = rr_csi * cas.R_CSI + rd_csi
    r_recsi1 = rr_recsi * cas.R_ReCSI + rd_recsi

    r_asym_acf = 0 if r_asym_acf is None else r_asym_acf

    return {
        'pdx0': pdx0,
        'pdx1': pdx1,
        'r_csi_acf': r_csi1 - cas.R_CSI,
        'r_recsi_acf': r_recsi1 - cas.R_ReCSI,
        'r_asym_acf': r_asym_acf,
        'ppm': ppm
    }


class ModelIND(Model):
    def update_parameters(self, pars):
        pars = Model.update_parameters(self, pars)
        pars['cas'].reform(pdx0=0.45, pdx1=0.45)
        return pars

    def simulate_intv(self, p, ys0, intv=None):
        if intv is not None:
            p['intv'] = intv
        elif 'intv' in p:
            del p['intv']

        ys, ms, msg = self.simulate_onward(p, ys0, 2036, 1)
        return ys, ms, msg
