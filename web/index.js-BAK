seamless_read_cells = {
  "text": [],
  "json": []
}
seamless_write_cells = {
  "text": [],
  "json": []
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
    return {}
  },
  watch: {
    
  }
})

vm = app.$mount('#app')
