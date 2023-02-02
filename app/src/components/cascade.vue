<template>
  <div>
    <div class="container">
      <h2>Setting: {{Country}} ({{Setting}})</h2>

      <div class="row">
        <div class="col-md-6">
          <h4>Setting</h4>
          <div class="form-check form-check-inline">
            <input class="form-check-input" type="radio" name="inlineRadioOptions" id="ind" value="IND" v-model="Setting" />
            <label class="form-check-label" for="ind">IND</label>
          </div>

          <div class="form-check form-check-inline">
            <input class="form-check-input" type="radio" name="inlineRadioOptions" id="zaf" value="ZAF" v-model="Setting" />
            <label class="form-check-label" for="zaf">ZAF</label>
          </div>


          <h4>Baseline</h4>
          <label for="r_prdx0" class="form-label">PrDx at first care-seeking: {{`${(Inputs.pdx0 * 100).toFixed()}%`}}</label>
          <input type="range" class="form-range" id="r_prdx0" v-model="Inputs.pdx0" min="0.1" max="0.95" step="0.01">
          <label for="r_prdx1" class="form-label">PrDx at revisiting care-seeking: {{`${(Inputs.pdx1 * 100).toFixed()}%`}}</label>
          <input type="range" class="form-range" id="r_prdx1" v-model="Inputs.pdx1" min="0.1" max="0.95" step="0.01">
        </div>
        <div class="col-md-6">
          <h4>TB care improvements</h4>
          <label for="r_r0" class="form-label">Relative rate of initial care-seeking: {{`${(Inputs.rr_csi * 1).toFixed(1)}`}}</label>
          <input type="range" class="form-range" id="r_r0" v-model="Inputs.rr_csi" min="1" max="10" step="0.1">
          <label for="r_r1" class="form-label">Relative rate of revisiting: {{`${(Inputs.rr_recsi * 1).toFixed(1)}`}}</label>
          <input type="range" class="form-range" id="r_r1" v-model="Inputs.rr_recsi" min="1" max="10" step="0.1">
          <label for="r_dx0" class="form-label">Odds ratio of PrDx at initial care-seeking: {{`${(Inputs.or_pdx0 * 1).toFixed(1)}`}}</label>
          <input type="range" class="form-range" id="r_dx0" v-model="Inputs.or_pdx0" min="1" max="20" step="0.1">
          <label for="r_dx1" class="form-label">Odds ratio of PrDx at revisiting: {{`${(Inputs.or_pdx1 * 1).toFixed(1)}`}}</label>
          <input type="range" class="form-range" id="r_dx1" v-model="Inputs.or_pdx1" min="1" max="20" step="0.1">
        </div>
      </div>
      <div class="row">
        <div class="col-md-4">
          <ul class="list-group">
            <li class="list-group-item">PrDx at first care-seeking: <br>
              {{`${(Stats[0].pdx0 * 100).toFixed()}%`}} -> {{`${(Stats[1].pdx0 * 100).toFixed()}%`}}</li>
            <li class="list-group-item">PrDx at re-visting: <br>
              {{`${(Stats[0].pdx1 * 100).toFixed()}%`}} -> {{`${(Stats[1].pdx1 * 100).toFixed()}%`}}</li>
            <li class="list-group-item">Duration between visits: <br>
              {{`${(Stats[0].dur * 12).toFixed(1)}`}} months -> {{`${(Stats[1].dur * 12).toFixed(1)}`}} months</li>
          </ul>
        </div>
        <div class="col-md-8">
          <b-tabs content-class="mt-3">
            <b-tab title="Visualisation" active>
              <vis :results="Results"></vis>
            </b-tab>
            <b-tab title="Values">
              <div class="row">
                <div class="col-md-6">
                  <div class="card">
                    <div class="card-header">
                      Cascade
                    </div>
                    <div class="card-body">
                      <h5 class="card-title">Baseline</h5>
                      <p>Detection from first care-seeking</p>
                      <div>{{fmt_ps(Results.Cas0.FromDx0)}}</div>
                      <p>Case detection</p>
                      <div>{{fmt_ps(Results.Cas0.Detect)}}</div>
                      <p>Case notification</p>
                      <div>{{fmt_ps(Results.Cas0.Report)}}</div>
                      <h5 class="card-title">Improved TB care</h5>
                      <p>Detection from first care-seeking</p>
                      <div>{{fmt_ps(Results.Cas1.FromDx0)}}</div>
                      <p>Case detection</p>
                      <div>{{fmt_ps(Results.Cas1.Detect)}}</div>
                      <p>Case notification</p>
                      <div>{{fmt_ps(Results.Cas1.Report)}}</div>
                    </div>
                  </div>
                </div>
                <div class="col-md-6">
                  <div class="card">
                    <div class="card-header">
                      Delays
                    </div>
                    <div class="card-body">
                      <h5 class="card-title">Baseline</h5>
                      <p>Patient Delay</p>
                      <div>{{fmt_durs(Results.Delay0.Pat)}}</div>
                      <p>System Delay</p>
                      <div>{{fmt_durs(Results.Delay0.Sys)}}</div>
                      <p>Total Delay</p>
                      <div>{{fmt_durs(Results.Delay0.Tot)}}</div>
                      <h5 class="card-title">Improved TB care</h5>
                      <p>Patient Delay</p>
                      <div>{{fmt_durs(Results.Delay1.Pat)}}</div>
                      <p>System Delay</p>
                      <div>{{fmt_durs(Results.Delay1.Sys)}}</div>
                      <p>Total Delay</p>
                      <div>{{fmt_durs(Results.Delay1.Tot)}}</div>

                    </div>
                  </div>
                </div>
              </div>
            </b-tab>
          </b-tabs>
        </div>
      </div>
    </div>

  </div>
</template>

<script>
import pars from "../assets/pars.json";
import {quantileSeq, exp, median} from "mathjs";
import vis from "./visual.vue";


export default {
  name: "v-cascade",
  components: {
    vis
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
    const sts = {
      pdx0: 0.4,
      pdx1: 0.7,
      dur: 1 /// quantileSeq(p_sel.map(p => p.r_recsi), 0.5)
    };

    return {
      Setting: 'IND',
      Pars: pars,
      Pars_sel: pars['IND'],
      ReformedPars0: [],
      ReformedPars1: [],
      Inputs: {
        pdx0: 0.4,
        pdx1: 0.7,
        or_pdx0: 1.2,
        or_pdx1: 1.2,
        rr_csi: 1,
        rr_recsi: 1
      },
      Stats: [sts, sts],
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
          FromDx0: [0, 0, 0],
          Detect: [0, 0, 0],
          Report: [0, 0, 0],
        },
        Cas1: {
          FromDx0: [0, 0, 0],
          Detect: [0, 0, 0],
          Report: [0, 0, 0],
        }
      }
    }
  },
  watch: {
    Inputs: {
      handler() {
        this.update_values();
      },
      deep: true
    },
    Setting: {
      handler() {
        this.Pars_sel = this.Pars[this.Setting]
        this.update_values();
      },
      deep: true
    }
  },
  mounted() {
    this.update_values();
  },
  methods: {
    update_values() {
      this.reform_pars0();
      this.update_cascade0();

      this.reform_pars1();
      this.update_cascade1();
    },
    reform_pars0() {
      this.ReformedPars0 = this.Pars_sel.map(p => {
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
        r.pdx0 = parseFloat(this.Inputs.pdx0);
        r.pdx1 = parseFloat(this.Inputs.pdx1);
        r.r_csi = p.r_aware;
        const det0 = r.r_csi * r.pdx0 * r.prv_s;
        const det1 = det - det0;
        // const fn0 = r.r_csi * (1 - r.pdx0) * r.prv_s;
        r.r_recsi = det1 / (r.pdx1 * r.prv_c);
        r.r_sym = p.r_sym;
        r.adr = p.adr;
        r.p_under = p.p_under;
        return r;
      }).filter(r => r.r_recsi > 0);

      this.Stats[0] = {
        pdx0: median(this.ReformedPars0.map(p => p.pdx0)),
        pdx1: median(this.ReformedPars0.map(p => p.pdx1)),
        dur: median(this.ReformedPars0.map(p => 1 / p.r_recsi)),
      }
    },
    update_cascade0() {
      const calc = this.ReformedPars0.map(p => {
        const dur_a = 1 / (p.r_sym + p.ra);
        const dur_s = 1 / (p.r_csi + p.rs);
        const dur_c = 1 / (p.r_recsi * p.pdx1 + p.rc);

        let pd0 = p.r_csi * p.pdx0 * dur_s;
        let pd1 = p.r_recsi * p.pdx1 * dur_c;
        const k = pd0 + pd1;

        const delays = {
          Pat: dur_s,
          Sys: dur_c * pd1 / k
        };

        delays.Tot = delays.Pat + delays.Sys;

        const cas = {
          Detect: p.r_sym * dur_a * (p.r_csi * dur_s * (p.pdx0 + (1 - p.pdx0) * p.r_recsi * p.pdx1 * dur_c))
        };
        cas.Report = cas.Detect * (1 - p.p_under);
        cas.FromDx0 = pd0 / k;

        return {delays, cas};
      });

      this.Results.Delay0.Pat = quantileSeq(calc.map(c => c.delays.Pat), [0.25, 0.5, 0.75]);
      this.Results.Delay0.Sys = quantileSeq(calc.map(c => c.delays.Sys), [0.25, 0.5, 0.75]);
      this.Results.Delay0.Tot = quantileSeq(calc.map(c => c.delays.Tot), [0.25, 0.5, 0.75]);
      this.Results.Cas0.FromDx0 = quantileSeq(calc.map(c => c.cas.FromDx0), [0.25, 0.5, 0.75]);
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
      }).filter(r => r.r_recsi > 0);

      this.Stats[1] = {
        pdx0: median(this.ReformedPars1.map(p => p.pdx0)),
        pdx1: median(this.ReformedPars1.map(p => p.pdx1)),
        dur: median(this.ReformedPars1.map(p => 1 / p.r_recsi))
      }
    },
    update_cascade1() {
      const calc = this.ReformedPars1.map(p => {
        const dur_a = 1 / (p.r_sym + p.ra);
        const dur_s = 1 / (p.r_csi + p.rs);
        const dur_c = 1 / (p.r_recsi * p.pdx1 + p.rc);

        let pd0 = p.r_csi * p.pdx0 * dur_s;
        let pd1 = p.r_recsi * p.pdx1 * dur_c;
        const k = pd0 + pd1;

        const delays = {
          Pat: dur_s,
          Sys: dur_c * pd1 / k
        };

        delays.Tot = delays.Pat + delays.Sys;

        const cas = {
          Detect: p.r_sym * dur_a * (p.r_csi * dur_s * (p.pdx0 + (1 - p.pdx0) * p.r_recsi * p.pdx1 * dur_c))
        };
        cas.Report = cas.Detect * (1 - p.p_under);
        cas.FromDx0 = pd0 / k;

        return {delays, cas};
      });

      this.Results.Delay1.Pat = quantileSeq(calc.map(c => c.delays.Pat), [0.25, 0.5, 0.75]);
      this.Results.Delay1.Sys = quantileSeq(calc.map(c => c.delays.Sys), [0.25, 0.5, 0.75]);
      this.Results.Delay1.Tot = quantileSeq(calc.map(c => c.delays.Tot), [0.25, 0.5, 0.75]);
      this.Results.Cas1.FromDx0 = quantileSeq(calc.map(c => c.cas.FromDx0), [0.25, 0.5, 0.75]);
      this.Results.Cas1.Detect = quantileSeq(calc.map(c => c.cas.Detect), [0.25, 0.5, 0.75]);
      this.Results.Cas1.Report = quantileSeq(calc.map(c => c.cas.Report), [0.25, 0.5, 0.75]);
    },
    fmt_durs(x) {
      return `${(x[1] * 12).toFixed(1)}` +
          ` (${(x[0] * 12).toFixed(1)} - ` +
          `${(x[2] * 12).toFixed(1)})` +
          'months';
    },
    fmt_ps(x) {
      return `${(x[1] * 100).toFixed()}%` +
          ` (${(x[0] * 100).toFixed()}% - ` +
          `${(x[2] * 100).toFixed()}%)`;
    }
  }
}
</script>

<style scoped>

</style>