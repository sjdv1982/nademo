from seamless.highlevel import Context, Cell, Transformer
from seamless.stdlib  import stdlib
import json

ctx = Context()
ctx.include(stdlib.map.map_dict)
ctx.load_vault("vault")

ctx.sub = Context()
ctx.translate()

""" # not yet implemented...
graph = json.load(open("graph/nademo.seamless"))
ctx.sub.set_graph(graph, mounts=True, shares=False)
"""

graph = json.load(open("graph/nademo.seamless"))
ctx2 = Context()
ctx2.load_vault("vault")
ctx2.set_graph(graph, mounts=True, shares=False)
ctx.sub = ctx2


del ctx.sub.na_resid
del ctx.sub.protein_resid
del ctx.sub.get_all_properties
del ctx.sub.all_properties
del ctx.sub.get_df_pdb
del ctx.sub.df_pdb
del ctx.sub.get_df_x3dna
del ctx.sub.df_x3dna
del ctx.sub.get_df_stackings
del ctx.sub.df_stackings
del ctx.sub.get_plot
del ctx.sub.plot

ctx.compute(0.1)
ctx.sub.inp = Cell("mixed").set({
    "pdb_code": ctx.sub.pdb_code.value,
    "na_chain": ctx.sub.na_chain.value,
    "protein_chain": ctx.sub.protein_chain.value,
})
ctx.sub.inp2 = Cell()
ctx.sub.inp2 = ctx.sub.inp
ctx.translate()
ctx.sub.pdb_code = ctx.sub.inp2.pdb_code
ctx.sub.na_chain = ctx.sub.inp2.na_chain
ctx.sub.protein_chain = ctx.sub.inp2.protein_chain
ctx.sub.result = ctx.sub.stackings
ctx.sub.result.celltype = "mixed"

ctx.compute()
print(ctx.status)
print(ctx.sub.stackings.value)

data = {
    "1B7F": ["P", "A"],
    "1A9N": ["Q", "B"]
}

inp = {}
for k,v in data.items():
    inp[k] = {
        "pdb_code": k,
        "na_chain": v[0],
        "protein_chain": v[1],
    }

ctx.inp = inp
ctx.result = Cell()
ctx.keyorder = Cell("plain")
ctx.map_pdbs = ctx.lib.map_dict(
    context_graph=ctx.sub,
    inp = ctx.inp,
    keyorder0 = [],
    keyorder = ctx.keyorder,
    result = ctx.result,
    elision = True,
    elision_chunksize = 2
)
# Elision and keyorder are overkill for few PDBs
# But for tens or hundreds of thousands, they save many minutes
#  because the construction of the dependency graph itself is elided

# Execution proceeds in parallel
# You should configure your cluster or cloud to handle Seamless transformations
# efficiently
# => execute the code where the data is, not the other way around
# But this is moot in this case where the data is downloaded from the PDB website

ctx.compute()
print(ctx.status)

# Only do this if the data fits (easily) in memory
# If it does not, define a hash pattern for ctx.result
from pprint import pprint
pprint(ctx.result.value.unsilk)

# Modify the calculation interactively
ctx.sub.calc_properties.mount("calc-properties.py")
ctx.translate()
# do ctx.translate(force=True) whenever you modify calc-properties.py
# => The PDBs do not get re-downloaded

# Modify the input PDB codes+chains interactively
# => Only the new PDBs get recomputed
# ...