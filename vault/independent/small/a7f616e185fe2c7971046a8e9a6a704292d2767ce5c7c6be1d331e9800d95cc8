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
