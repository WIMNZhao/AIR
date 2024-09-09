import rdkit
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def FingerPrint(smi,map_root,radius=8):
    # features are element, valence and pharmacohore based [F,#1] ---> F
    # hydrophobic feature [#6] ---> [#6,#16,F,Cl,Br,I]
    feature_smt = {'[#6,#16,F,Cl,Br,I]':0,'[#7]':1,'[#8]':2,'F':3,'Cl':4,'Br':5,'I':6,'[F,Cl,Br,I]':7,'[#16]':8,
                   '[CR]':9,'[c,$(C=*),$([#7X3H0+0][c,$(C=O)])]':10,'S':11,'s':12,
                   '[$(O=[C,c])]':13,'[$(O=S)]':14,'[$(O=N)]':15,'[$([OH]C=O)]':16,'[OH;!$([OH]C=O)]':17,
                   '[nX2]':18,'[nH,$([NH][C,S]=O),$([NH]c)]':19,'[$(N#C)]':20,
                   '[$([NH2][CX4]),$([NH]([CX4])[CX4]),$(N([CX4])([CX4])[CX4])]':21,
                   '[OH,nH,NH,NH2,NH3]':22,'[O,o,nX2,NX2]':23,
                   '[CX4H0]':24,'[CX4H1]':25,'[$([R](!@[*])@[R]!@[*])]':26,'[$([R](@[R]!@[*])@[R]!@[*])]':27,
                   '[$([#6R](@[#7,#8,#16])@[#7,#8,#16])]':28,'[$([#6R](@[#7,#8,#16])@[#6])]':29,
                   '[$([R](!@[*])(@[#7,#8,#16])@[#7,#8,#16])]':30,'[$([R](!@[*])(@[#7,#8,#16])@[#6])]':31,
                   '[$(C#*)]':32,                   

               }
    features = {}
    for smt in feature_smt:
        ff = Chem.MolFromSmarts(smt)
        features[ff] = feature_smt[smt]

    # read r-group
    mol = Chem.MolFromSmiles(smi)
    if mol:
       mol = Chem.RemoveHs(mol)
    else:
       mol = Chem.MolFromSmiles(smi,sanitize=False)
       mol = Chem.RemoveHs(mol,sanitize=False)
       mol.UpdatePropertyCache(strict=False)
       Chem.GetSymmSSSR(mol)

    # featurize
    feat = [[] for atom in mol.GetAtoms()]
    for smt in features:
        matches = mol.GetSubstructMatches(smt)
        for match in matches:
            for idx in match:
                feat[idx].append(features[smt])
    
    # chiral features
    f_chiral = {'R':len(features),'S':len(features)+1}
    nchiral = AllChem.FindMolChiralCenters(mol)
    for ch in nchiral:
        feat[ch[0]].append(f_chiral[ch[1]])

    # atom topological enviroment for N & O
    ate_per_atom = 2
    ate = 4
    topo_cutoff = ate_per_atom + 1
    f_ate = [[0]*ate for atom in mol.GetAtoms()]
    for ii, root in enumerate(mol.GetAtoms()):
        am = root.GetAtomicNum()
        if am == 8 or am == 7:
           root.SetIntProp('dd',0)
           atoms_visited = {}
           atoms_visited[root.GetIdx()] = ''
           tree = []
           tree.append(root)
           flag = False
           for search in tree:
               for atom in search.GetNeighbors():
                   if atom.GetIdx() not in atoms_visited:
                      dd = search.GetIntProp('dd')+1
                      if dd > topo_cutoff:
                         flag = True
                         break
                      tree.append(atom)
                      atoms_visited[atom.GetIdx()] = ''
                      atom.SetIntProp('dd',dd)
               if flag: break
           for atom in tree[1:]:
               dd = atom.GetIntProp('dd')
               if dd > 1:
                  f_ate[ii][(am-7)*ate_per_atom+dd-2] += 1

    # topolical distance
    roots = []
    for atom in mol.GetAtoms():
        atom.SetIntProp('distance',999)
        map_num = atom.GetAtomMapNum()
        if map_num == map_root:
           roots.append(atom)

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

    # Coding
    nf_no_ate = len(features) + len(f_chiral) 
    nfeatures = nf_no_ate + ate
    fp = np.zeros(radius*nfeatures,dtype=np.float32)
    for atom in mol.GetAtoms():
        dist = atom.GetIntProp('distance')
        if dist > 0 and dist <= radius:
           idx = atom.GetIdx()
           index = (dist-1)*nfeatures 
           for ii in feat[idx]:
               fp[index + ii] += 1
           index += nf_no_ate 
           for jj in range(ate):
               fp[index+jj] += f_ate[idx][jj]

    return fp
   
    
    #------------------------------------------------------------------------------
    # check
    for atom in mol.GetAtoms():
        print(atom.GetSmarts(),feat[atom.GetIdx()],atom.GetIntProp('distance'))
    for i in range(len(fp)):
        if i%nfeatures == 0:
           print()
        print(int(fp[i]),end=' ')
    print()      
    

if __name__ == "__main__":
     
     smi  = 'OCCCCN[*:1]'
     xx = FingerPrint(smi,1)





