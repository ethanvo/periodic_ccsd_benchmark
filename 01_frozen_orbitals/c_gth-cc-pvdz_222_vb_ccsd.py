#!/usr/bin/env python
import numpy as np
from pyscf.gto.basis import parse_nwchem
from pyscf.pbc import gto, scf, cc
from pyscf.pbc.cc import eom_kccsd_rhf
from pyscf.pbc.tools import pyscf_ase, lattice
from pyscf.pbc.cc.eom_kccsd_rhf import _IMDS
from pyscf.pbc.scf import chkfile
import os
import h5py

prefix = os.getcwd()
imds = None
    
##############################
# Create a "Cell"
##############################

cell = gto.Cell()
# Candidate formula of solid: c, si, sic, bn, bp, aln, alp, mgo, mgs, lih, lif, licl
formula = 'c'
scaled_center = [0.0, 0.0, 0.0]
ase_atom = lattice.get_ase_atom(formula)
cell.atom = pyscf_ase.ase_atoms_to_pyscf(ase_atom)
cell.a = ase_atom.cell
cell.unit = 'B'
cell.basis = {'C': parse_nwchem.load('/burg/berkelbach/users/eav2136/builds/ccgto/basis/gth-hf-rev/cc-pvdz-lc.dat', 'C')}
cell.pseudo = 'gthhfrev'
cell.verbose = 7
cell.build()
    
##############################
#  K-point SCF 
##############################

kdensity = 2
kmesh = [kdensity, kdensity, kdensity] 
kpts = cell.make_kpts(kmesh, scaled_center=scaled_center)
mymf = scf.KRHF(cell, kpts=kpts, exxdiv='ewald')
mymf = mymf.density_fit()

# save ERIs if not exist
cderi_h5name = 'c_gth-cc-pvdz_222_vb' + '_cderi.h5'
cderi_h5name = os.path.join(prefix, cderi_h5name)
mymf.with_df._cderi_to_save = cderi_h5name

if os.path.isfile(cderi_h5name):
    mymf.with_df._cderi = cderi_h5name

chkfile_name = 'c_gth-cc-pvdz_222_vb' + '.chk'
chkfile_name = os.path.join(prefix, chkfile_name)
mymf.chkfile = chkfile_name

if os.path.isfile(chkfile_name):
    chk_cell, chk_scf = chkfile.load_scf(chkfile_name)
    mymf.mo_coeff = chk_scf['mo_coeff']
    mymf.mo_energy = chk_scf['mo_energy']
    mymf.mo_occ = chk_scf['mo_occ']
    mymf.e_tot = chk_scf['e_tot']
else:
    ekrhf = mymf.kernel()
    
##############################
# K-point CCSD
##############################

cc_h5name = 'c_gth-cc-pvdz_222_vb' + '_cc.h5'
cc_h5name = os.path.join(prefix, cc_h5name)

frozen = None
mycc = cc.KRCCSD(mymf, frozen=frozen)

if kdensity > 3:
    mycc.max_cycle = 2

if os.path.isfile(cc_h5name):
    print('Reading t1, t2 from {}'.format(cc_h5name))
    t1 = None
    t2 = None
    with h5py.File(cc_h5name, 'r') as fin:
        t1 = fin['t1'][:]
        t2 = fin['t2'][:]
    mycc.t1 = t1
    mycc.t2 = t2

ekrcc, t1, t2 = mycc.kernel()

if os.path.isfile(cc_h5name):
    print('Saving t1, t2 to {}'.format(cc_h5name))
    with h5py.File(cc_h5name, 'a') as fout:
        fout['t1'][:] = t1
        fout['t2'][:] = t2
else:
    print('Saving t1, t2 to {}'.format(cc_h5name))
    with h5py.File(cc_h5name, 'w') as fout:
        fout.create_dataset('t1', data=t1)
        fout.create_dataset('t2', data=t2)
    
