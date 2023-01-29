<template>
  <div>
    <div>Setting: {{Country}} ({{Setting}})</div>
    <div>Sample size: {{Pars.length}}</div>
    <div>{{Results.Cas0}}</div>
    <div>{{Results.Cas1}}</div>
    <div>{{Results.Delay0}}</div>
    <div>{{Results.Delay1}}</div>
  </div>
</template>

<script>
import pars from "../public/pars.json";
import {quantileSeq, exp} from "mathjs";

export default {
  name: "cascade",
  props: {
    Setting: String
  },
  computed: {
    Country() {
      switch (this.Setting) {
        case "ZAF":
          return "South Africa";
        case "IND":
          return "India";
        default:
          return "None"
      }
    }
  },
  data() {
    return {
      Pars: pars[this.Setting],
      ReformedPars0: [],
      ReformedPars1: [],
      Inputs: {
        pdx0: 0.4,
        pdx1: 0.7,
        or_pdx0: 1,
        or_pdx1: 1.2,
        rr_csi: 1.2,
        rr_recsi: 1
      },
      Results: {
        Delay0: {
          Pat: [0, 0, 0],
          Sys: [0, 0, 0],
          Tot: [0, 0, 0]
        },
        Delay1: {
          Pat: [0, 0, 0],
          Sys: [0, 0, 0],
          Tot: [0, 0, 0]
        },
        Cas0: {
          Detect: [0, 0, 0],
          Report: [0, 0, 0],
        },
        Cas1: {
          Det: [0, 0, 0],
          Report: [0, 0, 0],
        }
      }
    }
  },
  mounted() {
    this.reform_pars0();
    this.update_cascade0();

    this.reform_pars1();
    this.update_cascade1();
  },
  methods: {
    reform_pars0() {
      this.ReformedPars0 = this.Pars.map(p => {
        const r = {};

        r.prv = p.prv0 * exp(- p.adr * (2023 - p.Year0));
        r.ra = p.r_sc + p.r_death_a + p.r_death_bg;
        r.rs = p.r_sc + p.r_death_s + p.r_death_bg;
        r.rc = p.r_sc + p.r_death_s + p.r_death_bg;

        const a0 = (r.rs + p.r_aware - p.adr) / p.r_sym;
        const c0 = p.r_aware / (r.rc + p.r_det - p.adr);

        const pr_a = a0 / (a0 + 1 + c0);
        const pr_s = 1 / (a0 + 1 + c0);
        const pr_c = c0 / (a0 + 1 + c0);

        r.prv_a = r.prv * pr_a;
        r.prv_s = r.prv * pr_s;
        r.prv_c = r.prv * pr_c;

        const det = p.r_det * r.prv_c;
        r.pdx0 = this.Inputs.pdx0;
        r.pdx1 = this.Inputs.pdx1;
        r.r_csi = p.r_aware;
        const det0 = r.r_csi * r.pdx0 * r.prv_s;
        const det1 = det - det0;
        const fn0 = r.r_csi * (1 - r.pdx0) * r.prv_s;
        r.r_recsi = det1 / (r.pdx1 * r.prv_c);
        r.r_sym = p.r_sym;
        r.adr = p.adr;
        r.p_under = p.p_under;
        return r;
      })
    },
    update_cascade0() {
      const calc = this.ReformedPars0.map(p => {
        const dur_a = 1 / (p.r_sym + p.ra);
        const dur_s = 1 / (p.r_csi + p.rs);
        const dur_c = 1 / (p.r_recsi * p.pdx1 + p.rc);
        const NoSys = p.r_csi * p.pdx0 * dur_s;

        const delays = {
          Pat: dur_s,
          Sys: (1 - NoSys) * dur_c
        };
        console.log(delays);
        delays.Tot = delays.Pat + delays.Sys;

        const cas = {
          Detect: p.r_sym * dur_a * (p.r_csi * dur_s * (p.pdx0 + (1 - p.pdx0) * p.r_recsi * p.pdx1 * dur_c))
        };
        cas.Report = cas.Detect * (1 - p.p_under);

        return {delays, cas};
      });

      this.Results.Delay0.Pat = quantileSeq(calc.map(c => c.delays.Pat), [0.25, 0.5, 0.75]);
      this.Results.Delay0.Sys = quantileSeq(calc.map(c => c.delays.Sys), [0.25, 0.5, 0.75]);
      this.Results.Delay0.Tot = quantileSeq(calc.map(c => c.delays.Tot), [0.25, 0.5, 0.75]);
      this.Results.Cas0.Detect = quantileSeq(calc.map(c => c.cas.Detect), [0.25, 0.5, 0.75]);
      this.Results.Cas0.Report = quantileSeq(calc.map(c => c.cas.Report), [0.25, 0.5, 0.75]);
    },
    reform_pars1() {
      this.ReformedPars1 = this.ReformedPars0.map(p => {
        const r = Object.assign({}, p);
        let x;
        x = r.pdx0 / (1 - r.pdx0) * this.Inputs.or_pdx0;
        r.pdx0 = x / (1 + x);

        x = r.pdx1 / (1 - r.pdx1) * this.Inputs.or_pdx1;
        r.pdx1 = x / (1 + x);
        r.r_csi *= this.Inputs.rr_csi;
        r.r_recsi *= this.Inputs.rr_recsi;
        return r;
      })
    },
    update_cascade1() {
      const calc = this.ReformedPars1.map(p => {
        const dur_a = 1 / (p.r_sym + p.ra);
        const dur_s = 1 / (p.r_csi + p.rs);
        const dur_c = 1 / (p.r_recsi * p.pdx1 + p.rc);
        const NoSys = p.r_csi * p.pdx0 * dur_s;

        const delays = {
          Pat: dur_s,
          Sys: (1 - NoSys) * dur_c
        };
        console.log(delays);
        delays.Tot = delays.Pat + delays.Sys;

        const cas = {
          Detect: p.r_sym * dur_a * (p.r_csi * dur_s * (p.pdx0 + (1 - p.pdx0) * p.r_recsi * p.pdx1 * dur_c))
        };
        cas.Report = cas.Detect * (1 - p.p_under);

        return {delays, cas};
      });

      this.Results.Delay1.Pat = quantileSeq(calc.map(c => c.delays.Pat), [0.25, 0.5, 0.75]);
      this.Results.Delay1.Sys = quantileSeq(calc.map(c => c.delays.Sys), [0.25, 0.5, 0.75]);
      this.Results.Delay1.Tot = quantileSeq(calc.map(c => c.delays.Tot), [0.25, 0.5, 0.75]);
      this.Results.Cas1.Detect = quantileSeq(calc.map(c => c.cas.Detect), [0.25, 0.5, 0.75]);
      this.Results.Cas1.Report = quantileSeq(calc.map(c => c.cas.Report), [0.25, 0.5, 0.75]);
    },
  }
}
</script>

<style scoped>

</style>