import sys
from rdkit import Chem
from rdkit.Chem.rdchem import EditableMol

def BondTypeToInt(bondtype):
    if(bondtype == Chem.rdchem.BondType.SINGLE):
       return 1   
    if(bondtype == Chem.rdchem.BondType.DOUBLE):
       return 2
    if(bondtype == Chem.rdchem.BondType.TRIPLE):
       return 3
    if(bondtype == Chem.rdchem.BondType.AROMATIC):
       return 4
    print()
    print('================Please Update Bond Type Definition===============')
    print(bondtype)
    print()
    sys.exit(0)


def Assemble(smi):
    mol = Chem.MolFromSmiles(smi,sanitize = False)
    mol.UpdatePropertyCache(strict=False)

    core_anchors = []
    r_attachment = []
    core_maps = []
    r_maps = []
    cc = {}
    core = True

    for atom in mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0:
           if map_num in cc:
              core = False
           else:
              cc[map_num] = ''
           if core:
              core_anchors.append(atom.GetIdx())
           else:
              r_attachment.append(atom.GetIdx())

    for idx in core_anchors:
        atom = mol.GetAtomWithIdx(idx)
        map_num = atom.GetAtomMapNum()
        nbrs = atom.GetNeighbors()
        for nbr in nbrs:
            idx_n = nbr.GetIdx()
            bondtype = mol.GetBondBetweenAtoms(idx,idx_n).GetBondType()
            bt = BondTypeToInt(bondtype)   
            core_maps.append([idx_n,map_num,bt])

    rr = {}
    for idx in r_attachment:
        atom = mol.GetAtomWithIdx(idx)
        map_num = atom.GetAtomMapNum()
        nbrs = atom.GetNeighbors()
        for nbr in nbrs:
            idx_n = nbr.GetIdx()
            bondtype= mol.GetBondBetweenAtoms(idx,idx_n).GetBondType()
            bt = BondTypeToInt(bondtype)   
            if idx_n not in rr:
               r_maps.append([idx_n,map_num,bt])
               rr[idx_n] = ''

    #print (core_maps)
    #print (r_maps)
    
    bond_join = []
    for ii in core_maps:
        for jj in r_maps:
            if ii[1] == jj[1]:
               if ii[2] == 1 and jj[2] == 4:
                  bt = 4
               else:
                  bt = ii[2]
               bond_join.append([ii[0],jj[0],bt])

    # Make an editable molecule and add bonds between atoms with correspoing AtomMapNum
    em = EditableMol(mol)
    for bond in bond_join:
        start_atm = bond[0]
        end_atm = bond[1]
        bondtype = bond[2]
        if bondtype == 2:
           em.AddBond(start_atm, end_atm, order=Chem.rdchem.BondType.DOUBLE)
        elif bondtype == 3:
           em.AddBond(start_atm, end_atm, order=Chem.rdchem.BondType.TRIPLE)
        elif bondtype == 4:
           em.AddBond(start_atm, end_atm, order=Chem.rdchem.BondType.AROMATIC)
        else:
           em.AddBond(start_atm, end_atm, order=Chem.rdchem.BondType.SINGLE)
    final_mol = em.GetMol()

    # Nuke all of the dummy atoms
    final_mol = Chem.DeleteSubstructs(final_mol, Chem.MolFromSmarts('[#0]'))

    # remove the AtomMapNum values
    for atom in final_mol.GetAtoms():
        atom.SetAtomMapNum(0)

    # sanity check
    sanitFail = Chem.SanitizeMol(final_mol, catchErrors=True)
    if sanitFail: return None

    final_mol = Chem.RemoveHs(final_mol)

    return final_mol


def BatchAssemble(core_rgs):
    mols = []
    for cpd in core_rgs:
        uniq = {}
        rgs = cpd.split('|')
        for rr in rgs[1:]:
            uniq[rr] = ''
        smi = rgs[0] 
        for ff in uniq: smi = smi + '.' + ff 
        mol = Assemble(smi)
        if mol:
           mols.append(Chem.MolToSmiles(mol,isomericSmiles = True))
        else:
           mols.append(None) 
    return mols


if __name__ == "__main__":
    # First fragment must be the core
    # !!!         
    smi = 'c(c(c1ccccc1)[*:2])[*:1].[H]C([H])([H])[*:1].[H]c(o:[*:2])c([H]):[*:1]'
    #smi = 'C1=[*:1]CCCC1.N([*:1])=[*:1]'    
    #smi = 'C1=[*:2]C([*:1])CCC1.[H]C([H])(C([H])([H])[*:1])[*:1].[H]C([*:2])=[*:2]'

    welded_mol = Assemble(smi)
    print(Chem.MolToSmiles(welded_mol))



