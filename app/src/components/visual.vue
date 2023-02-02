<template>
  <div>
    <div class="row">
      <div class="col-md-6">
        <div class="card">
          <div class="card-header">
            Cascade
          </div>
          <div class="card-body">
            <p>Detection from first care-seeking</p>
            <bar vid="v_cas_dx0" height="50px" :chart-data="Cas.FromDx0" upper="1"></bar>
            <p>Case detection</p>
            <bar vid="v_cas_det" height="50px" :chart-data="Cas.Detect" upper="1"></bar>
            <p>Case notification</p>
            <bar vid="v_cas_rep" height="50px" :chart-data="Cas.Report" upper="1"></bar>
          </div>
        </div>
      </div>
      <div class="col-md-6">
        <div class="card">
          <div class="card-header">
            Delays
          </div>
          <div class="card-body">
            <p>Patient Delay</p>
            <bar vid="v_del_pat" height="50px" :chart-data="Delay.Pat" :upper="results.Delay0.Tot[1]"></bar>
            <p>System Delay</p>
            <bar vid="v_del_sys" height="50px" :chart-data="Delay.Sys" :upper="results.Delay0.Tot[1]"></bar>
            <p>Total Delay</p>
            <bar vid="v_del_tot" height="50px" :chart-data="Delay.Tot" :upper="results.Delay0.Tot[1]"></bar>
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import bar from "./bar.vue";

export default {
  name: "visual-res",
  components: {
    bar
  },
  props: {
    results: {
      type: Object,
      required: true
    }
  },
  data() {
    return {
      Delay: {
        Pat: [{value: 0, text: "0"}, {value: 0, text: "0"}],
        Sys: [{value: 0, text: "0"}, {value: 0, text: "0"}],
        Tot: [{value: 0, text: "0"}, {value: 0, text: "0"}]
      },
      Cas: {
        FromDx0: [{value: 0, text: "0"}, {value: 0, text: "0"}],
        Detect: [{value: 0, text: "0"}, {value: 0, text: "0"}],
        Report: [{value: 0, text: "0"}, {value: 0, text: "0"}],
      }
    }
  },
  watch: {
    results: {
      handler() {
        this.update();
      },
      deep: true
    }
  },
  mounted() {
    this.update()
  },
  methods: {
    update() {
      ['Pat', 'Sys', 'Tot'].forEach(k => {
        let v = this.results.Delay0[k][1];
        console.log(v)
        this.Delay[k][0].value = v;
        this.Delay[k][0].text = this.fmt_dur(v) + "months";

        v = this.results.Delay1[k][1]
        this.Delay[k][1].value = v;
        this.Delay[k][1].text = this.fmt_dur(v) + "months";
      });

      ['FromDx0', 'Detect', 'Report'].forEach(k => {
        let v = this.results.Cas0[k][1];
        this.Cas[k][0].value = v;
        this.Cas[k][0].text = this.fmt_p(v);

        v = this.results.Cas1[k][1]
        this.Cas[k][1].value = v;
        this.Cas[k][1].text = this.fmt_p(v);
      })
    },
    fmt_dur(x) {
      return `${(x * 12).toFixed(1)}`;
    },
    fmt_durs(x) {
      return `${this.fmt_dur(x[1])} ( ${this.fmt_dur(x[0])} - ${this.fmt_dur(x[2])}) months`

    },
    fmt_p(x) {
      return `${(x * 100).toFixed()}%`
    },
    fmt_ps(x) {
      return `${this.fmt_p(x[1])} ( ${this.fmt_p(x[0])} - ${this.fmt_p(x[2])})`
    }
  }
}
</script>

<style scoped>

</style>