#!/usr/bin/env python

from pyscf.gto.basis import parse_nwchem
import json

def dump(data, filename):
    with open(filename, 'w') as fout:
        fout.write(json.dumps(data, indent=4))

def get_basis(foru

def write_basis(formula, material_components, basis_file):
    basis = '{'
    for element_index, element in enumerate(material_components[formula]):
        if element_index > 0:
            basis = basis + ', '
        basis = basis + '\'' + element + '\': parse_nwchem.load(\'' + basis_file + '\', \'' + element + '\')'
    basis = basis + '}'
    return basis

def calculate_norbs(formula, basis):
    # Create cell
    cell = gto.Cell()
    ase_atom = lattice.get_ase_atom(formula)
    cell.atom = pyscf_ase.ase_atoms_to_pyscf(ase_atom)
    cell.a = ase_atom.cell
    cell.unit = 'B'
    cell.basis = basis
    cell.pseudo = 'gthhfrev'
    cell.verbose = 7
    cell.build()

    # Calculate number of molecular orbitals and virtuals
    nmo = cell.nao
    nocc = cell.nelectron // 2
    nvir = nmo - nocc

material_components = {
        'c':['C'], 
        'si':['Si'], 
        'sic':['Si', 'C'], 
        'bn':['B', 'N'], 
        'bp':['B', 'P'], 
        'aln':['Al', 'N'], 
        'alp':['Al', 'P'], 
        'mgo':['Mg', 'O'], 
        'mgs':['Mg', 'S'], 
        'lih':['Li', 'H'], 
        'lif':['Li', 'F'], 
        'licl':['Li', 'Cl']}

vb_scaled_centers = {
        'c':[0., 0., 0.], 
        'si':[0., 0., 0.], 
        'sic':[0., 0., 0.], 
        'bn':[0., 0., 0.], 
        'bp':[0., 0., 0.], 
        'aln':[0., 0., 0.], 
        'alp':[0., 0., 0.], 
        'mgo':[0., 0., 0.], 
        'mgs':[0., 0., 0.], 
        'lih':[0.5, 0.,  0.5], 
        'lif':[0., 0., 0.], 
        'licl':[0.0530303,  0.0530303,  0.10606061]}

cb_scaled_centers = {
        'c':[0.36363636, 0.,         0.36363636], 
        'si':[0.3989899, 0.,        0.3989899], 
        'sic':[0.5, 0.,  0.5], 
        'bn':[0.5, 0.,  0.5], 
        'bp':[0.38888889, 0.,         0.38888889], 
        'aln':[0.0199005, 0.0199005, 0.       ], 
        'alp':[0.5, 0.,  0.5], 
        'mgo':[0., 0., 0.], 
        'mgs':[0.5, 0.,  0.5], 
        'lih':[0.5, 0.,  0.5], 
        'lif':[0., 0., 0.], 
        'licl':[0., 0., 0.]}

vb_nroots = {
        'c':3, 
        'si':3, 
        'sic':3, 
        'bn':3, 
        'bp':3, 
        'aln':1, 
        'alp':3, 
        'mgo':3, 
        'mgs':3, 
        'lih':1, 
        'lif':3, 
        'licl':1}

cb_nroots = {
        'c':1, 
        'si':1, 
        'sic':1, 
        'bn':1, 
        'bp':1, 
        'aln':1, 
        'alp':1, 
        'mgo':1, 
        'mgs':1, 
        'lih':1, 
        'lif':1, 
        'licl':1}

ccpvdz = '/burg/berkelbach/users/eav2136/builds/ccgto/basis/gth-hf-rev/cc-pvdz-lc.dat'
ccpvtz = '/burg/berkelbach/users/eav2136/builds/ccgto/basis/gth-hf-rev/cc-pvdz-lc.dat'
ccpvqz = '/burg/berkelbach/users/eav2136/builds/ccgto/basis/gth-hf-rev/cc-pvdz-lc.dat'

basis_sets = [ccpvdz, ccpvtz, ccpvqz]

materials = ['c', 'si', 'sic', 'bn', 'bp', 'aln', 'alp', 'mgo', 'mgs', 'lih', 'lif', 'licl']

for formula in materials:
    for basis in basis_sets:
        

material = dict([
('formula', 'c'),
('basis', {'C': parse_nwchem.load('/burg/berkelbach/users/eav2136/builds/ccgto/basis/gth-hf-rev/cc-pvdz-lc.dat', 'C')})
])
filename = 'data/c_cc-pvdz-lc.json'
dump(material, filename)

