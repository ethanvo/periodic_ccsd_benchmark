#!/usr/bin/env python

from pyscf.pbc import gto
from pyscf.pbc.tools import lattice, pyscf_ase
from pyscf.gto.basis import parse_nwchem
import json

def load(filename):
    with open(filename, 'r') as fin:
        data = json.loads(fin.read())
    return data
filename = 'c_cc-pvdz-lc.json'

material = load(filename)
# Create cell
cell = gto.Cell()
formula = material['formula']
ase_atom = lattice.get_ase_atom(formula)
cell.atom = pyscf_ase.ase_atoms_to_pyscf(ase_atom)
cell.a = ase_atom.cell
cell.unit = 'B'
cell.basis = material['basis']
cell.pseudo = 'gthhfrev'
cell.verbose = 7
cell.build()

# Calculate number of molecular orbitals and virtuals
nmo = cell.nao
nocc = cell.nelectron // 2
nvir = nmo - nocc

print(nmo)
print(nvir)
