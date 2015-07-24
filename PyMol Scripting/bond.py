# Basic Python Script using PyMol API to bond all the alpha carbons in a protein chain.
# Scripted to visualize knots and slipknots in a protein chain.
__author__ = 'Shashank'
from pymol import cmd, stored

# An array list used to store all the residue numbers which will be fetched.
residues = []
# set stick_radius=0.25
# iterate command used to fetch all the alpha carbon residue numbers from the given PDB
cmd.iterate('(n. CA)', 'residues.append(resi)')

n = residues.__len__()
i = 0
# Bond command in the API is used to connect all the consecutive CA atoms in the chain
while i<n-1:
    cmd.bond(residues[i]+'/CA',residues[i+1]+'/CA')
    print "Bonding "+residues[i]+" And "+residues[i+1]
    i = i+1
# Coloring the bonds using the Spectrum command
cmd.load("1ALK_A_test.pdb")
cmd.spectrum("1ALK_A_test")
# cmd.label('(n. CA)', "resi")
