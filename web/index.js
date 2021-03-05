seamless_read_cells = {
  "text": [
    "plot"
  ],
  "json": [
    "pdb_code",
    "na_chain",
    "protein_chain",
    "na_resid",
    "protein_resid",
    "stackings_annotated"
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
      "stackings_annotated": {}
    }
  },
  watch: {
    pdb_code: function (value) {
      seamless_update("pdb_code", value, "json")
      loadNGL()
    },
    na_chain: function (value) {
      seamless_update("na_chain", value, "json")
      loadNGL()
    },
    protein_chain: function (value) {
      seamless_update("protein_chain", value, "json")
      loadNGL()
    },
    na_resid: function (value) {
      seamless_update("na_resid", value, "json")
      loadNGL()
    },
    protein_resid: function (value) {
      seamless_update("protein_resid", value, "json")
      loadNGL()
    },

  }
})

vm = app.$mount('#app')

stage = new NGL.Stage("ngl-viewer");
// Handle window resizing
window.addEventListener( "resize", function( event ){
  stage.handleResize();
}, false );

function loadNGL() {
  if (!app.pdb_code) return;
  if (!app.na_resid) return;
  if (!app.na_chain) return;
  if (!app.protein_resid) return;
  if (!app.protein_chain) return;
  stage.removeAllComponents()
  var pdb_url = "rcsb://" + app.pdb_code
  var selection=`(${app.na_resid} and :${app.na_chain}) or (${app.protein_resid} and :${app.protein_chain})`
  stage.loadFile( pdb_url, { defaultRepresentation: true } ).then( function( pdb ){
    pdb.addRepresentation('ball+stick', {"sele":selection, "color":'blue'})
  });
}
