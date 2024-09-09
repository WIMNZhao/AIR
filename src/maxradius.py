import rdkit
import numpy as np
from rdkit import Chem

def MaxRadius(smi,map_root):
    # read r-group
    mol = Chem.MolFromSmiles(smi)
    if mol:
       mol = Chem.RemoveHs(mol)
    else:
       mol = Chem.MolFromSmiles(smi,sanitize=False)
       mol = Chem.RemoveHs(mol,sanitize=False)

    # topolical distance
    roots = []
    for atom in mol.GetAtoms():
        atom.SetIntProp('distance',999)
        map_num = atom.GetAtomMapNum()
        if map_num == map_root:
           roots.append(atom)

    radius = []

    for root in roots:
        if root.GetIntProp('distance') != 999: continue
        root.SetIntProp('distance',0)
        atoms_visited = {}
        atoms_visited[root.GetIdx()] = ''
        tree = []
        tree.append(root)
        for search in tree: 
            for atom in search.GetNeighbors():
                if atom.GetIdx() not in atoms_visited:
                   tree.append(atom)
                   atoms_visited[atom.GetIdx()] = ''
                   atom.SetIntProp('distance',search.GetIntProp('distance')+1)
        radius.append(tree[-1].GetIntProp('distance'))
    return max(radius)


def MaxRadiusList(r_grps):
    nr = len(r_grps)
    radius = [ [] for i in range(nr) ]
    for i, grps in enumerate(r_grps):
        for r in grps:
            radius[i].append(MaxRadius(r,i+1))
   
    r_radius = []
    for rs in radius:
        r_radius.append(max(rs))

    return r_radius
            


if __name__ == "__main__":
     
     smi  = '[H][*:1]'
     xx = MaxRadius(smi,1)
     print(xx)




