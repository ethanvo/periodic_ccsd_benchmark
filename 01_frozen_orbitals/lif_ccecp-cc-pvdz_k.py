#!/usr/bin/env python
import numpy as np
from pyscf.pbc import gto, scf, cc
from pyscf.pbc.tools import pyscf_ase, lattice
import sys

##############################
# Create a "Cell"
##############################

cell = gto.Cell()
# Candidate formula of solid: c, si, sic, bn, bp, aln, alp, mgo, mgs, lih, lif, licl
formula = "lif"
g_vbmax = 3
g_cbmin = 1
vb_scaled_center = [0., 0., 0.]
cb_scaled_center = [0., 0., 0.]
ase_atom = lattice.get_ase_atom(formula)
cell.atom = pyscf_ase.ase_atoms_to_pyscf(ase_atom)
cell.a = ase_atom.cell
cell.unit = 'B'
cell.basis = 'ccecp-cc-pvdz'
cell.pseudo = "ccecp"
cell.verbose = 7
cell.build()

##############################
#  K-point SCF 
##############################

kdensity = int(sys.argv[1])
kmesh = [kdensity, kdensity, kdensity] 
kpts = cell.make_kpts(kmesh, scaled_center=vb_scaled_center)
mymf = scf.KRHF(cell, kpts=kpts, exxdiv="ewald")
mymf = mymf.density_fit()
ekrhf = mymf.kernel()

##############################
# K-point CCSD
##############################

frozen = None
mycc = cc.KRCCSD(mymf, frozen=frozen)
ekrcc, t1, t2 = mycc.kernel()

##############################
# EOM-IP/EA-KRCCSD
##############################

# number of roots requested
# index(indices) of targetted k point(s) 
kptlist = [0]

eip, vip = mycc.ipccsd(nroots=g_vbmax, kptlist=kptlist)

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

frozen = None
mycc = cc.KRCCSD(mymf, frozen=frozen)
ekrcc, t1, t2 = mycc.kernel()

##############################
# EOM-IP/EA-KRCCSD
##############################

# number of roots requested
# index(indices) of targetted k point(s) 
kptlist = [0]

eea, vea = mycc.eaccsd(nroots=g_cbmin, kptlist=kptlist)

band_gap = (np.amin(eip) + np.amin(eea)) * 27.211386245988
print('Band gap:', band_gap, 'eV')
