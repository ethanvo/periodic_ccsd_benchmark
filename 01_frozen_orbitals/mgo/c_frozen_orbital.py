#!/usr/bin/env python
import numpy as np
from pyscf.gto.basis import parse_nwchem
from pyscf.pbc import gto, scf, cc
from pyscf.pbc.cc import eom_kccsd_rhf
from pyscf.pbc.tools import pyscf_ase, lattice
from pyscf.pbc.cc.eom_kccsd_rhf import _IMDS
from pyscf.pbc.scf import chkfile
import h5py

au2ev = 27.211386245988

##############################
# Create a "Cell"
##############################

cell = gto.Cell()
# Candidate formula of solid: c, si, sic, bn, bp, aln, alp, mgo, mgs, lih, lif, licl
formula = "c"
g_vbmax = 3
g_cbmin = 1
vb_scaled_center = [0., 0., 0.]
cb_scaled_center = [0.36363636, 0.,         0.36363636]
ase_atom = lattice.get_ase_atom(formula)
cell.atom = pyscf_ase.ase_atoms_to_pyscf(ase_atom)
cell.a = ase_atom.cell
cell.unit = 'B'
cell.basis = {'C': parse_nwchem.load('/burg/berkelbach/users/eav2136/builds/ccgto/basis/gth-hf-rev/cc-pvdz-lc.dat', 'C')}
cell.pseudo = 'gthhfrev'
cell.verbose = 7
cell.build()

nmo = cell.nao
nvir = cell.nao - cell.nelectron//2
kdensity = 2
kmesh = [kdensity, kdensity, kdensity]

ccsd_vb_array = []
ccsd_cb_array = []
eip_array = []
eea_array = []
band_gap_array = []

for nfrozen in range(nvir):

    ##############################
    #  K-point SCF 
    ##############################

    kpts = cell.make_kpts(kmesh, scaled_center=vb_scaled_center)
    mymf = scf.KRHF(cell, kpts=kpts, exxdiv="ewald")
    mymf = mymf.density_fit()
    ekrhf = mymf.kernel()

    ##############################
    # K-point CCSD
    ##############################

    frozen = list(range(nmo - 1, nmo - 1 - nfrozen, -1))
    if nfrozen == 0:
        frozen = None
    print(frozen)
    mycc = cc.KRCCSD(mymf, frozen=frozen)
    ekrcc, t1, t2 = mycc.kernel()
    ccsd_vb_array.append(ekrcc)

    ##############################
    # EOM-IP/EA-KRCCSD
    ##############################

    # number of roots requested
    # index(indices) of targetted k point(s) 
    kptlist = [0]

    eip, vip = mycc.ipccsd(nroots=g_vbmax, kptlist=kptlist)
    print('EIP Energy: {} eV'.format(np.amin(eip) * au2ev)) 
    eip_array.append(eip)

    ##############################
    # New scaled center
    ##############################

    kpts = cell.make_kpts(kmesh, scaled_center=cb_scaled_center)
    mymf = scf.KRHF(cell, kpts=kpts, exxdiv="ewald")
    mymf = mymf.density_fit()
    ekrhf = mymf.kernel()

    ##############################
    # K-point CCSD
    ##############################

    mycc = cc.KRCCSD(mymf, frozen=frozen)
    ekrcc, t1, t2 = mycc.kernel()
    ccsd_cb_array.append(ekrcc)

    ##############################
    # EOM-IP/EA-KRCCSD
    ##############################

    # number of roots requested
    # index(indices) of targetted k point(s) 
    kptlist = [0]

    eea, vea = mycc.eaccsd(nroots=g_cbmin, kptlist=kptlist)
    print('EEA Energy: {} eV'.format(np.amin(eea) * au2ev)) 
    eea_array.append(eea)

    band_gap = (np.amin(eip) + np.amin(eea)) * au2ev 
    band_gap_array.append(band_gap)
    print('Band gap:', band_gap, 'eV')

filename = formula + '_frozen_data.h5py'
with h5py.File(filename, 'w') as fout:
    fout.create_dataset('ccsd_vb_array', data=ccsd_vb_array)
    fout.create_dataset('ccsd_cb_array', data=ccsd_vb_array)
    fout.create_dataset('eip_array', data=eip_array)
    fout.create_dataset('eea_array', data=eea_array)
    fout.create_dataset('band_gap_array', data=band_gap_array)
