import os
import sys
import argparse
import pandas as pd
import numpy as np
import warnings


def ReadMolFile(mollines):
    [nats, nbonds] = map(int, (mollines[3][0:3], mollines[3][3:6]))
    cooslines = mollines[4:4 + nats]
    coos = {}
    atypes = {}
    for i in range(nats):
        els = cooslines[i].split()
        coos[i + 1] = [float(e) for e in els[0:3]]
        atypes[i + 1] = els[3]
    bondlines = mollines[4 + nats:4 + nats + nbonds]
    bonds = {'BI': [], 'BJ': [], 'RIJ': [], 'UID': []}
    for line in bondlines:
        [bi, bj] = map(int, [line[0:3], line[3:6]])
        bonds['BI'].append(bi)
        bonds['BJ'].append(bj)
    return (coos, atypes, bonds)


def InfoFromPDB(pdb_file):
    pdb_lines = open(pdb_file, 'r').readlines()
    atom_lines = []
    conn_lines = []
    for l in pdb_lines:
        if 'ATOM' in l:
            atom_lines.append(l[0:29])
        elif 'CONECT' in l:
            conn_lines.append(l)
    return(atom_lines, conn_lines)


def WriteToPDB(atom_lines, conn_lines, new_coos):
    opdb = open('LIG_final.pdb', 'w+')
    opdb.write('REMARK LIGPARGEN GENERATED PDB FILE\n')
    num = 0
    for i in range(len(atom_lines)):
        jc = new_coos[i+1]
        opdb.write('%s %8.3f%8.3f%8.3f\n' %
                   (atom_lines[i], jc[0], jc[1], jc[2]))
    opdb.write('TER \n')
    for i in conn_lines:
        opdb.write('%s\n' % i.strip())
    opdb.write('END\n')
    opdb.close()
    return None


def PDB2MOL(pdbfile):
    os.system('babel -ipdb %s -omol2 %s.mol2' % (pdbfile, pdbfile[:-4]))
    return ('%s.mol2' % pdbfile[:-4])


def MCSAlignMolecules(ref_mol, ali_mol):
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    from rdkit.Chem import rdFMCS
    from rdkit.Chem.rdFMCS import FindMCS, AtomCompare, BondCompare
    '''
    Do not sanitize the molecules, RDKit will freak out and give errors
    And All we want is to do MCSS, we dont care much about health of molecule
    '''
    mol1 = Chem.MolFromMol2File(ref_mol, removeHs=False, sanitize=False)
    mol2 = Chem.MolFromMol2File(ali_mol, removeHs=False, sanitize=False)
    _fmcs_params = dict(maximizeBonds=False, threshold=1.0, timeout=60,
                        verbose=False, matchValences=True,
                        ringMatchesRingOnly=True, completeRingsOnly=True,
                        atomCompare=AtomCompare.CompareAny,
                        bondCompare=BondCompare.CompareAny)
    try:
        mcs = rdFMCS.FindMCS([mol1, mol2], **_fmcs_params)
    except ValueError:
        print('\n Max Common Substructure calculation \n failed for this molecule!! \n Please be judicious ')
        sys.exit()
    core = Chem.MolFromSmarts(mcs.smartsString)
    match1 = mol1.GetSubstructMatch(core)
    match2 = mol2.GetSubstructMatch(core)
    from rdkit.Chem import AllChem
    AllChem.AlignMol(mol2, mol1, atomMap=list(zip(match2, match1)))
    Chem.MolToMolFile(mol2, 'aligned.mol', kekulize=False)
    return mol2

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
         
        @author: Leela Sriram Dodda, Jorgensen Lab @ Yale
        @email:  leela.dodda@yale.edu
        
        MIT License
        Usage: 
        python alignpdb.py -f LIG.pdb -t ligand_wH.pdb  

        Copyright (c) 2017 Leela S. Dodda
	REQUIREMENTS:
	Preferably Anaconda python with following modules
	argparse
	RDKit for MCS overlap  
        openbabel/babel executable
	"""
    )
    parser.add_argument(
        "-f", "--fromm",  help="molecule from LPG generated Position",required=True)
    parser.add_argument(
        "-t", "--to",  help="molecule in Bound Position ",required=True)
    args = parser.parse_args()
    if args:
       try:
           Fr_pdb_name = args.fromm
           To_pdb_name = args.to
           Fr_mol_name = PDB2MOL(Fr_pdb_name)
           To_mol_name = PDB2MOL(To_pdb_name)
           aligned_mol = MCSAlignMolecules(ref_mol=To_mol_name, ali_mol=Fr_mol_name)
           final_coos, final_atypes, final_bonds = ReadMolFile(
               mollines=open('aligned.mol', 'r').readlines())
           ats, cons = InfoFromPDB(pdb_file=Fr_pdb_name)
           assert len(final_coos.keys()) == len(ats), 'SOMETHINGS NOT RIGHT'
           os.system('/bin/rm aligned.mol %s %s'%(Fr_mol_name,To_mol_name))
           WriteToPDB(ats, cons, final_coos)
           print('Completed Successfully')
       except (RuntimeError, TypeError, NameError):
           print('Error Occured')
