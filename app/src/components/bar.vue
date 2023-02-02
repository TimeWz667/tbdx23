<template>
<div>
  <svg :id="vid"></svg>
</div>
</template>

<script>
import * as d3 from "d3";

export default {
  name: "vis-bar",
  props: {
    vid: {
      type: String,
      required: true
    },
    chartData: {
      required: true,
      type: Array
    },
    upper: {
      default: 12
    },
    height: {
      default: "20px",
      type: String
    }
  },
  data() {
    return {
      svg: null,
      g: null,
      bars: null,
      x: null,
      y: null
    }
  },
  mounted() {
    this.initialise_vis();
    window.addEventListener("resize", this.resize);
    this.$nextTick(this.update_vis);
  },
  watch: {
    chartData: {
      handler() {
        this.update_vis();
      },
      deep: true
    },
    upper: {
      handler() {
        this.update_vis();
      },
      deep: true
    }
  },
  methods: {
    initialise_vis() {
      this.svg = d3
          .select("#" + this.vid)
          .style("width", "100%")
          .style("height", this.height);

      this.width = this.svg.node().parentNode.clientWidth;

      this.g = this.svg.append("g")

      const ht = this.svg.node().parentNode.clientHeight;

      this.x = d3.scaleLinear()
          .domain([0, this.upper])
          .range([0, this.width * 0.7])
      this.y = d3.scaleBand()
          .range([0, ht])
          .domain([0, 1])
          .padding(.1);
    },
    update_vis() {
      this.x.domain([0, this.upper]);

      this.svg
        .selectAll("g.bar")
        .data(this.chartData)
      .join(
          enter => {
            enter.append("g").attr("class", "bar")
            .call(g => {
              g.append("rect")
                  .attr("x", 0)
                  .attr("y", (d, i) => this.y(i))
                  .attr("width", d => this.x(d.value))
                  .attr("height", this.y.bandwidth())
                  .attr("fill", (d, i) => (i === 0 ? "red" : "green"));

              g.append("text")
                  .attr("x", d => this.x(d.value) + 10)
                  .attr("y", (d, i) => this.y(i) + this.y.bandwidth() * 0.7)
                  .text(d => d.text)
                  .attr("text-anchor", "right")
                  .attr("alignment-baseline", "middle")
            })
          },
          update => {
            update.call(g => {
              g.select("rect")
                .attr("width", d => this.x(d.value))

              g.select("text")
                  .attr("x", d => this.x(d.value) + 10)
                  .text(d => d.text)
            })
          }
      )
    },
    resize() {
      const width = this.svg.node().parentNode.clientWidth;
      const height = this.svg.node().parentNode.clientHeight;
      this.width = width;

      if (this.x !== null) {
        this.x.range([0, width]);
      }
      if (this.y !== null) {
        this.y.range([0, height]);
      }
    }
  }
}
</script>

<style scoped>

</style>