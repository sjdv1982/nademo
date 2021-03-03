"""
Routines to parse the atoms of a PDB into a structured Numpy array
Biopython is used for the initial parsing

The resulting array is C-aligned, so could be passed as a C struct to a compiled language

Sjoerd de Vries, 2020-2021
"""

import warnings
import Bio.PDB
from Bio.PDB.StructureBuilder import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)
import numpy as np


atomic_dtype = [
    ("model", 'uint16'),
    ("hetero", "S1"),
    ("name", "S4"),
    ("altloc","S1"),
    ("resname", "S3"),
    ("chain","S1"),
    ("index", 'uint32'),
    ("icode", "S1"),
    ("resid", 'uint16'),
    ("x", 'float32'),
    ("y", 'float32'),
    ("z", 'float32'),
    ("occupancy", 'float32'),
    ("bfactor", 'float32'),
    ("segid", "S4"),
    ("element", "S2")
]

atomic_dtype = [tuple(x) for x in atomic_dtype]
atomic_dtype = np.dtype(atomic_dtype, align=True)

def print_atom(atom):
    if atom.ndim == 1:
        return print_atoms(atom)
    if atom.ndim > 0:
        raise TypeError
    fields = atom.dtype.fields
    result = {}
    for field in fields:
        result[field] = atom[field]
    return result

def print_atoms(atoms):
    if atoms.ndim == 0:
        return print_atom(atoms)
    if atoms.ndim > 1:
        raise TypeError
    result = []
    for atom in atoms:
        result.append(print_atom(atom))
    return result

def parse_pdb(pdbdata):
    from io import StringIO
    pdb_obj = StringIO(pdbdata)

    p = Bio.PDB.PDBParser()
    obj = "pdb"
    struc = p.get_structure(obj, pdb_obj)
    natoms = len(list(struc.get_atoms()))
    atomstate = np.zeros(natoms,dtype=atomic_dtype)

    count = 0
    for modelnr, model in enumerate(struc.get_models()):
        atomlist = list(model.get_atoms())
        atomlist.sort(key=lambda atom: atom.serial_number)
        for atom in atomlist:
            residue = atom.get_parent()
            hetero, resid, icode = residue.get_id()
            segid = residue.segid
            resname = residue.resname
            chainid = residue.get_parent().id
            aa = atomstate[count]
            aa["model"] = modelnr + 1
            aa["hetero"] = hetero
            aa["name"] = atom.name
            aa["altloc"] = atom.altloc
            aa["resname"] = resname
            aa["chain"] = chainid
            aa["index"] = atom.serial_number
            aa["icode"] = icode
            aa["resid"] = resid
            aa["x"] = atom.coord[0]
            aa["y"] = atom.coord[1]
            aa["z"] = atom.coord[2]
            occ = atom.occupancy
            if occ is None or occ < 0:
                occ = 0
            aa["occupancy"] = occ
            aa["segid"] = segid
            aa["element"] = atom.element
            aa["bfactor"] = atom.bfactor
            count += 1
    return atomstate

if __name__ == "__main__":
    import sys
    pdbfile = sys.argv[1]
    npyfile = sys.argv[2]
    pdbdata = open(pdbfile).read()
    result = parse_pdb(pdbdata)
    np.save( npyfile,result)