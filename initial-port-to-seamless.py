from seamless.highlevel import Context, Cell, Transformer, Module
ctx = Context()

ctx.pdb_code = "1B7F"
ctx.na_chain = "P"
ctx.protein_chain = "A"
ctx.na_resid = 5
ctx.protein_resid = 256

#####

import nglview
widget = nglview.NGLWidget()

def show_ngl(*args, **kwargs):
    pdb_code = ctx.pdb_code.value.unsilk
    na_chain = ctx.na_chain.value.unsilk
    protein_chain = ctx.protein_chain.value.unsilk
    na_resid = ctx.na_resid.value.unsilk
    protein_resid = ctx.protein_resid.value.unsilk

    widget.clear()
    widget.add_component("rcsb://" + pdb_code)
    selection='({0} and :{1}) or ({2} and :{3})'.format(na_resid, na_chain, protein_resid, protein_chain)
    widget.add_representation('ball+stick', selection=selection, color='blue')
    widget.center(selection)

await ctx.computation()
show_ngl()
display(widget)

ctx.pdb_code.traitlet().observe(show_ngl)
ctx.na_chain.traitlet().observe(show_ngl)
ctx.na_resid.traitlet().observe(show_ngl)
ctx.protein_chain.traitlet().observe(show_ngl)
ctx.protein_resid.traitlet().observe(show_ngl)

#####

def download_pdb(pdb_code):
    import urllib
    pdb_data = urllib.request.urlopen("https://files.rcsb.org/download/{}.pdb".format(pdb_code)).read().decode()
    return pdb_data

ctx.download_pdb = download_pdb
ctx.download_pdb.pdb_code = ctx.pdb_code
ctx.pdb_data = ctx.download_pdb

####

ctx.execute_x3dna = Transformer()
ctx.execute_x3dna.language = "docker"
ctx.execute_x3dna.docker_image = "x3dna"
ctx.execute_x3dna.code = "x3dna-dssr -i=pdb_data --json -o=RESULT"
ctx.execute_x3dna.pdb_data = ctx.pdb_data
ctx.x3dna_analysis = ctx.execute_x3dna

await ctx.computation()
ctx.x3dna_analysis.output()

####

def get_x3dna_nucleotides(x3dna_analysis, na_chain):
    return [nt for nt in x3dna_analysis["nts"] if nt["chain_name"] == na_chain]
ctx.get_x3dna_nucleotides = get_x3dna_nucleotides
ctx.get_x3dna_nucleotides.x3dna_analysis = ctx.x3dna_analysis
ctx.get_x3dna_nucleotides.na_chain = ctx.na_chain

####

ctx.x3dna_nucleotides = ctx.get_x3dna_nucleotides
ctx.x3dna_nucleotides.celltype = "plain"
await ctx.computation()

def get_df_x3dna(x3dna_nucleotides):
    import pandas as pd
    df_x3dna = pd.DataFrame(x3dna_nucleotides)
    return df_x3dna.to_html()

ctx.get_df_x3dna = get_df_x3dna
ctx.get_df_x3dna.x3dna_nucleotides = ctx.x3dna_nucleotides
ctx.get_df_x3dna.pins.x3dna_nucleotides.celltype = "plain"

ctx.df_x3dna = ctx.get_df_x3dna
await ctx.translation()
ctx.df_x3dna.mimetype = "text/html"
await ctx.translation()
display(ctx.df_x3dna.output())

####

ctx.parse_pdb = Module()
ctx.parse_pdb.mount("parse_pdb.py")

def get_parsed_pdb(pdb_data):
    parsed_pdb = parse_pdb.parse_pdb(pdb_data)
    return parsed_pdb

ctx.get_parsed_pdb = get_parsed_pdb
ctx.get_parsed_pdb.parse_pdb = ctx.parse_pdb
ctx.get_parsed_pdb.pdb_data = ctx.pdb_data
ctx.parsed_pdb = ctx.get_parsed_pdb
ctx.parsed_pdb.celltype = "binary"

await ctx.computation()
import parse_pdb
display(parse_pdb.print_atom(ctx.parsed_pdb.value[:2]))

####

def get_df_pdb(parsed_pdb):
    import numpy as np
    import pandas as pd
    df_pdb = pd.DataFrame(parsed_pdb)
    for col, dtype in df_pdb.dtypes.items():
        if dtype == np.object:  # Only process byte object columns.
            df_pdb[col] = df_pdb[col].apply(lambda x: x.decode("utf-8"))
    return df_pdb.to_html()

ctx.get_df_pdb = get_df_pdb
ctx.get_df_pdb.parsed_pdb = ctx.parsed_pdb
ctx.get_df_pdb.pins.parsed_pdb.celltype = "binary"
ctx.df_pdb = ctx.get_df_pdb

await ctx.translation()
ctx.df_pdb.mimetype = "text/html"

await ctx.translation()
display(ctx.df_pdb.output())

def select_chains(parsed_pdb, protein_chain, na_chain):
    protein_atoms = parsed_pdb[parsed_pdb["chain"]==protein_chain.encode()]
    na_atoms = parsed_pdb[parsed_pdb["chain"]==na_chain.encode()]
    return {
        "protein_atoms": protein_atoms,
        "na_atoms": na_atoms
    }
ctx.select_chains = select_chains
ctx.select_chains.parsed_pdb = ctx.parsed_pdb
ctx.select_chains.pins.parsed_pdb.celltype = "binary"
ctx.select_chains.protein_chain = ctx.protein_chain
ctx.select_chains.na_chain = ctx.na_chain
ctx.selected_chains = ctx.select_chains
ctx.protein_atoms = ctx.selected_chains.protein_atoms
ctx.na_atoms = ctx.selected_chains.na_atoms

def calculate_stacking_properties(protein_atoms, protein_resid, na_atoms, na_resid):
    import numpy as np
    import scipy.spatial.distance
    res_protein = protein_atoms[protein_atoms["resid"]==protein_resid]
    assert len(res_protein)
    aa = res_protein[0]["resname"].decode().strip()
    res_na = na_atoms[na_atoms["resid"]==na_resid]
    assert len(res_na)
    nuc = res_na[0]["resname"].decode().strip()[-1] # one-letter
    coor_res_protein = np.stack((res_protein["x"], res_protein["y"], res_protein["z"])).T
    coor_res_na = np.stack((res_na["x"], res_na["y"], res_na["z"])).T

    result = {}
    dist = scipy.spatial.distance.cdist(coor_res_protein, coor_res_na)
    result["closest_distance"] = dist.min()

    sidechains = {
        "PHE": ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']
    }
    sidechain_mask = np.isin(res_protein["name"], [name.encode() for name in sidechains[aa]])
    bases = {
        "U": ['C2', 'C4', 'C5', 'C6', 'N1', 'N3']
    }
    base_mask = np.isin(res_na["name"], [name.encode() for name in bases[nuc]])
    stacking_dist = dist[sidechain_mask][:,base_mask]
    result["mean_stacking_dist"] = stacking_dist.mean()
    result["std_stacking_dist"] = stacking_dist.std()

    return result

def calculate_all_properties(protein_atoms, protein_resid, na_atoms, na_resid, x3dna_nucleotides):
    stacking_properties = calculate_stacking_properties(protein_atoms, protein_resid, na_atoms, na_resid)
    x3dna_nucl = [nucl for nucl in x3dna_nucleotides if nucl["nt_resnum"] == na_resid]
    assert len(x3dna_nucl) == 1
    nucl_props = ["gamma", "delta", "chi"]
    result = {}
    for prop in nucl_props:
        result[prop] = x3dna_nucl[0][prop]
    result.update(stacking_properties)
    return result

ctx.calc_properties = Module()
import inspect
src = inspect.getsource(calculate_stacking_properties) \
+ "\n" \
+ inspect.getsource(calculate_all_properties)
ctx.calc_properties.code = src

def get_all_properties(protein_atoms, protein_resid, na_atoms, na_resid, x3dna_nucleotides):
    return calc_properties.calculate_all_properties(
        protein_atoms, protein_resid, na_atoms, na_resid, x3dna_nucleotides
    )
ctx.get_all_properties = get_all_properties
ctx.get_all_properties.calc_properties = ctx.calc_properties
ctx.get_all_properties.protein_atoms = ctx.protein_atoms
ctx.get_all_properties.pins.protein_atoms.celltype = "binary"
ctx.get_all_properties.protein_resid = ctx.protein_resid
ctx.get_all_properties.na_atoms = ctx.na_atoms
ctx.get_all_properties.pins.na_atoms.celltype = "binary"
ctx.get_all_properties.na_resid = ctx.na_resid
ctx.get_all_properties.x3dna_nucleotides = ctx.x3dna_nucleotides
ctx.get_all_properties.pins.x3dna_nucleotides.celltype = "plain"
ctx.all_properties = ctx.get_all_properties
ctx.all_properties.celltype = "plain"

def get_stackings(protein_atoms, na_atoms, x3dna_nucleotides):
    import numpy as np
    print
    from .calc_properties import calculate_all_properties
    all_protein_resids = np.unique(protein_atoms["resid"])
    all_na_resids = np.unique(na_atoms["resid"])
    stackings = []
    for curr_na_resid in all_na_resids:
        for curr_protein_resid in all_protein_resids:
            try:
                properties = calculate_all_properties(
                    protein_atoms, curr_protein_resid,
                    na_atoms, curr_na_resid,
                    x3dna_nucleotides
                )
            except (KeyError, AssertionError):
                continue
            properties["na_resid"] = int(curr_na_resid)
            properties["protein_resid"] = int(curr_protein_resid)
            stackings.append(properties)
    return stackings
display(ctx.all_properties.output())

ctx.get_stackings = get_stackings
ctx.get_stackings.calc_properties = ctx.calc_properties
ctx.get_stackings.protein_atoms = ctx.protein_atoms
ctx.get_stackings.pins.protein_atoms.celltype = "binary"
ctx.get_stackings.na_atoms = ctx.na_atoms
ctx.get_stackings.pins.na_atoms.celltype = "binary"
ctx.get_stackings.x3dna_nucleotides = ctx.x3dna_nucleotides
ctx.get_stackings.pins.x3dna_nucleotides.celltype = "plain"
ctx.stackings = ctx.get_stackings
ctx.stackings.celltype = "plain"

def get_df_stackings(stackings):
    import pandas as pd
    df_stackings = pd.DataFrame(stackings)
    return df_stackings.to_html()

ctx.get_df_stackings = get_df_stackings
ctx.get_df_stackings.stackings = ctx.stackings
ctx.get_df_stackings.pins.stackings.celltype = "plain"
ctx.df_stackings = ctx.get_df_stackings
await ctx.translation()
ctx.df_stackings.mimetype = "text/html"
await ctx.translation()
display(ctx.df_stackings.output())

def get_plot(stackings):
    from matplotlib import pyplot as plt
    import mpld3
    fig, ax = plt.subplots()
    ax.scatter(
        [stacking["chi"] for stacking in stackings],
        [stacking["closest_distance"] for stacking in stackings],
    )
    ax.set_xlabel('Chi')
    ax.set_ylabel('Closest distance')
    return mpld3.fig_to_html(fig)

ctx.get_plot = get_plot
ctx.get_plot.stackings = ctx.stackings
ctx.get_plot.pins.stackings.celltype = "plain"
ctx.plot = ctx.get_plot
await ctx.translation()
ctx.plot.mimetype = "text/html"
await ctx.translation()
display(ctx.plot.output())

####

ctx.save_graph("initial-port.seamless")
ctx.save_zip("initial-port.zip")