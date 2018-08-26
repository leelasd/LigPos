"""
RelativeFEP_LigAlign.py: 

This code is created as a tool to facilitate Relative FEP's with NAMD
Initial and Final states of ligand need to aligned to do Relative FEP
This code aligns the ligand in which ever requested. 
Initial and Final states of Ligands can be generated from LigParGen Server
http://zarbi.chem.yale.edu/ligpargen/
 
@author: Leela Sriram Dodda
@email:  leela.dodda@yale.edu

MIT License

Copyright (c) 2017 Leela S. Dodda

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from rdkit.Chem import AllChem
import os
from rdkit.Chem import rdMolAlign
import warnings
import argparse
from distutils import spawn


def GetAtomNames(ali_pdb):
    '''
    Get the names of atoms and connection info from the pdb file
    ali_pdb : pdb file name 
              type = string 
    '''
    from Bio.PDB import PDBParser
    from Bio import BiopythonWarning
    warnings.simplefilter('ignore', BiopythonWarning)
    parser_p = PDBParser(PERMISSIVE=1)
    lig = parser_p.get_structure('LIG', ali_pdb)
    anames = [at.get_id() for at in lig.get_atoms()]
    connects = []
    for line in open(ali_pdb, 'r').readlines():
        if 'CONECT' in line:
            connects.append(line.strip())
    return (anames, connects)


def pdb_prep(atoms, coos, connects, opdb, resid='UNK'):
    '''
    Create PDB from atoms, coordinates, connection info
    residue name and output filename
    atoms: string array
    coos : dictionary({SerialNo:[X,Y,Z]})
    connects: string array
    resid: string
    '''
    opdb = open(opdb, 'w+')
    opdb.write('REMARK LIGPARGEN GENERATED PDB FILE\n')
    num = 0
    for (i, j) in zip(atoms, coos):
        num += 1
        opdb.write('%-6s%5d %4s %3s  %4d    %8.3f%8.3f%8.3f\n' %
                   ('ATOM', num, i, resid, 1, coos[j][0], coos[j][1], coos[j][2]))
    opdb.write('TER \n')
    for i in range(len(connects)):
        opdb.write('%s \n' % connects[i])
    opdb.write('END')
    opdb.close()
    return None


def ReadMolFile(molfile):
    '''
    Read mol file and collect coordinates
    molfile: string,  Mol File of the desired molecule
    '''
    mollines = open(molfile, 'r').readlines()
    [nats, nbonds] = map(int, mollines[3].split()[0:2])
    cooslines = mollines[4:4 + nats]
    coos = {}
    atypes = {}
    for i in range(nats):
        els = cooslines[i].split()
        coos[i + 1] = [float(e) for e in els[0:3]]
        atypes[i + 1] = els[3]
    return (coos)


def AlignMolecules(ref_mol, ali_mol):
    '''
    Code to Align a random molecule with a reference molecule
    ref_mol: string, Reference Mol File
    ali_mol: string, Random Mol File
    '''
    mol1 = AllChem.MolFromMolFile(ref_mol, removeHs=False)
    mol2 = AllChem.MolFromMolFile(ali_mol, removeHs=False)
    out = rdMolAlign.GetO3A(mol2, mol1)
    rmsd, trans_matrix = out.Trans()
    AllChem.TransformMol(mol2, trans_matrix)
    return mol2

def MCSAlignMolecules(ref_mol, ali_mol):
    from rdkit.Chem import MCS
    mol1 = AllChem.MolFromMolFile(ref_mol, removeHs=False)
    mol2 = AllChem.MolFromMolFile(ali_mol, removeHs=False)
    mcs = MCS.FindMCS([mol1,mol2],ringMatchesRingOnly=True)
    core = AllChem.MolFromSmarts(mcs.smarts)
    smi_core = AllChem.MolFromSmiles(AllChem.MolToSmiles(core))
    match1 = mol1.GetSubstructMatch(smi_core)
    match2 = mol2.GetSubstructMatch(smi_core)
    AllChem.AlignMol(mol2,mol1,atomMap=list(zip(match2,match1)))   # <- m2 is aligned to m1
    return mol2

def PDB2MOL(pdb):
    '''
    Warning: 
    Need to have babel executable in the path
    pdb: string pdb file name with extension
    '''
    os.system('babel -ipdb %s -omol %s.mol' % (pdb, pdb[:-4]))
    return (pdb[:-4] + '.mol')

def ZMAT2MOL(pdb):
    '''
    Warning: 
    Need to have babel executable in the path
    pdb: string pdb file name with extension
    '''
    assert ('BOSSdir' in os.environ) and os.path.isfile((os.environ[
        'BOSSdir'] + '/scripts/xZmol')), 'Please Make sure $BOSSdir is defined \n xZmol and related files are in scripts directory of BOSS'
    execfile = os.environ['BOSSdir'] + '/scripts/xZmol > /tmp/olog'
    coma = execfile + ' ' + pdb[:-4] 
    os.system(coma)
    return (pdb[:-4] + '.mol')



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='RelativeFEP_LigAlign.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
        This code is created as a tool to facilitate Relative FEP's with NAMD
        Initial and Final states of ligand need to aligned to do Relative FEP
        This code aligns the ligand in which ever requested. 
        Initial and Final states of Ligands can be generated from LigParGen Server
        http://zarbi.chem.yale.edu/ligpargen/
         
        @author: Leela Sriram Dodda, Jorgensel Lab @ Yale
        @email:  leela.dodda@yale.edu
        
        MIT License
        Usage: 
        python RelativeFEP_LigAlign.py -r BNZ/UNK_1C46F7.pdb -a TOL/UNK_D65EA5.pdb -o test.pdb -i TOL 

        Copyright (c) 2017 Leela S. Dodda
	REQUIREMENTS:
	Preferably Anaconda python with following modules
	Biopython 
	argparse
	RDKit for using with SMILES code
        openbabel/babel executable
        BOSS (Optional but can give you better alignment)
	"""
    )
    parser.add_argument(
        "-r", "--ref",  help="reference molecule pdb file",required=True)
    parser.add_argument(
        "-a", "--ali",  help="molecule to be aligned in pdb file",required=True)
    parser.add_argument(
        "-o", "--output", help="output file name for the aligned molecule",required=True)
    parser.add_argument(
        "-i", "--id_res", type=str,help="residue name in aligned_molecule \n Default is UNK")
    parser.add_argument(
        "-b", "--use_boss", help="Use BOSS to create mol file instead of Babel",action="store_true")
    args = parser.parse_args()
    if args:
       assert (spawn.find_executable("babel") is not None),'Make sure you have openbabel executable'
       try:
          if args.use_boss:
             ref_mol = ZMAT2MOL(args.ref)
             ali_mol = ZMAT2MOL(args.ali)
          else:
             ref_mol = PDB2MOL(args.ref)
             ali_mol = PDB2MOL(args.ali)
          try:
              mol2 = MCSAlignMolecules(ref_mol, ali_mol)
          except (AttributeError): 
              print('aligning with Max Common Substructure failed \n aligning to minimize RMSD')
              mol2 = AlignMolecules(ref_mol, ali_mol)
          AllChem.MolToMolFile(mol2, filename='Aligned.mol')
          atom_names, connects = GetAtomNames(args.ali)
          coos = ReadMolFile(molfile='Aligned.mol')
          os.system('/bin/rm  Aligned.mol plt.pdb sum out optzmat log')
          pdb_prep(atom_names, coos, connects, resid=args.id_res, opdb=args.output)
          print('Completed Successfully')
       except (RuntimeError, TypeError, NameError):
           print('Error Occured')

