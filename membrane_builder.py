"""
    Author: Halil ibrahim özdemir
    Loc: Marmara University / Bioengineering
"""

from pdbfixer import PDBFixer
from simtk.openmm import app
from simtk import unit
import os
import parmed as pmd


def add_membrane(pdb_path, membrane_lipid_type='POPC', out_as=None):
    """
        Make a lipid bilayer for your protein easy.

            Parameters
            ----------

            pdb_path: Give your pdb whole path to this parameter

            membrane_lipid_type : Add POPC or POPE lipid membranes to your system.

            out_as: Give and extension list like ['psf', 'crd', 'gro']

            Example
            ----------

            add_membrane('protein.pdb', 'POPC', ['crd', 'gro'])

        """
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    print('\nAdding membrane:', membrane_lipid_type)
    app.PDBFile.writeFile(fixer.topology, fixer.positions, open("fixed.pdb", 'w'))
    fixer.addMembrane(lipidType=membrane_lipid_type,
                      membraneCenterZ=0 * unit.nanometer,
                      minimumPadding=1 * unit.nanometer,
                      positiveIon="Na+",
                      negativeIon="Cl-",
                      ionicStrength=0.0 * unit.molar)
    app.PDBFile.writeFile(fixer.topology, fixer.positions, open("fixed_membrane.pdb", 'w'), keepIds=True)

    if out_as is not None:
        struct = pmd.load_file('fixed_membrane.pdb')
        for i in out_as:
            try:
                print("Savind *.%s extension File" % i)
                struct.save('fixed_membrane.%s' % i)
            except:
                pass


