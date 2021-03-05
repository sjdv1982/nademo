from sympy import Plane, Point3D
from math import acos, sqrt
import numpy as np

def calculate_stacking_properties(protein_atoms, protein_resid, na_atoms, na_resid):
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
    result["aa"] = aa
    dist = scipy.spatial.distance.cdist(coor_res_protein, coor_res_na)
    result["closest_distance"] = dist.min()

    if result["closest_distance"] > 6:
        result["stacking_angle"] = 0
        result["min_stacking_distance"] = 10
        return result

    sidechains = {
        "PHE": ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        "TRP": ['CD', 'CZ2', 'CZ2', 'CE2', 'CE3', 'CH'],
        "TYR": ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        "HIS": ['CG', 'ND1', 'CD2', 'CE1', 'NE2'],
        "GLN": ['CD', 'OE1', 'NE2'],
        "ASP": ['OD1', 'OD2', 'CG'],
        "ARG": ['CZ', 'NH1', 'NH2'],
    }
    sidechain_mask = np.isin(res_protein["name"], [name.encode() for name in sidechains[aa]])
    bases = {
        "U": ['C2', 'C4', 'C5', 'C6', 'N1', 'N3'],
        "C": ['C2', 'C4', 'C5', 'C6', 'N1', 'N3'],
        "G": ['C2', 'C4', 'C5', 'C6', 'N1', 'N3', "N7", "C8","N9"],
        "A": ['C2', 'C4', 'C5', 'C6', 'N1', 'N3', "N7", "C8","N9"],
    }
    base_mask = np.isin(res_na["name"], [name.encode() for name in bases[nuc]])
    stacking_dist = dist[sidechain_mask][:,base_mask]
    result["min_stacking_distance"] = float(stacking_dist.min())
    #result["std_stacking_dist"] = stacking_dist.std()

    # 10.1261/rna.054924.115 (Suppl S8,S9) define stacking
    # by angle and min heavy-atom dist between the
    # cycles (F, Y, W, H) or cycles-triangles (Q, D, R):
    #   min_dist in [2.7-4.3A],
    #   angle in [0 - 20] for parall stacking
    #   angle in [70 - 90] for T-shape stacking.
    cycle_na = coor_res_na[base_mask]
    cycle_protein = coor_res_protein[sidechain_mask]

    plane_na = Plane(Point3D(cycle_na[0]), Point3D(cycle_na[1]), Point3D(cycle_na[2]))
    plane_protein = Plane(Point3D(cycle_protein[0]), Point3D(cycle_protein[1]), Point3D(cycle_protein[2]))

    rad = plane_na.angle_between(plane_protein).evalf()
    deg = round(rad*180/3.14159)
    d = min(deg, 180-deg)
    result["stacking_angle"] = float(d)

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
