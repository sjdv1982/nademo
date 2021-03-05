seamless_read_cells = {
  "text": [
    "plot",
    "df_stackings"
  ],
  "json": [
    "pdb_code",
    "na_chain",
    "protein_chain",
    "na_resid",
    "protein_resid"
  ]
}
seamless_write_cells = {
  "text": [],
  "json": [
    "pdb_code",
    "na_chain",
    "protein_chain",
    "na_resid",
    "protein_resid"
  ]
}

ctx = connect_seamless()
ctx.self.onsharelist = function(sharelist) {
    sharelist.forEach(element => {
        if (seamless_read_cells["json"].indexOf(element) >= 0) {
            ctx[element].onchange = function() {
                const v = JSON.parse(this.value)
                vm[element] = v
            }
        }
        else if (seamless_read_cells["text"].indexOf(element) >= 0) {
            ctx[element].onchange = function() {
                vm[element] = this.value
            }
        }
    })
}

function seamless_update(cell, value, encoding) {
  if (!ctx) return
  if (!ctx.self.sharelist) return
  if (ctx.self.sharelist.indexOf(cell) < 0) return
  if (encoding == "json") {
    ctx[cell].set(JSON.stringify(value))
  }
  else if (encoding == "text") {
    ctx[cell].set(value)
  }
}

const app = new Vue({
  vuetify: new Vuetify(),
  data() {
    return {
      "pdb_code": "",
      "na_chain": "",
      "protein_chain": "",
      "na_resid": 0,
      "protein_resid": 0,
      "plot": "",
      "df_stackings": ""
    }
  },
  watch: {
    pdb_code: function (value) {
    seamless_update("pdb_code", value, "json")
    },
    na_chain: function (value) {
    seamless_update("na_chain", value, "json")
    },
    protein_chain: function (value) {
    seamless_update("protein_chain", value, "json")
    },
    na_resid: function (value) {
    seamless_update("na_resid", value, "json")
    },
    protein_resid: function (value) {
    seamless_update("protein_resid", value, "json")
    },
    
  }
})

vm = app.$mount('#app')
