import rdkit
from rdkit import Chem

def CoreFromSmarts(smt):
    with open(smt) as f:
         for line in f:
              if line.startswith('!'): continue
              if line.strip():
                 core = Chem.MolFromSmarts(line.split()[0])
                 if core: return core

def MolMatchSubsFromSmi(smifile,core):
    mols = []
    prop = []
    with open(smifile) as f:
         for line in f:
             if line.startswith('!'): continue
             if line.strip():
                cells = line.split()
                mol = Chem.MolFromSmiles(cells[0])
                if mol:
                   xmol = Chem.AddHs(mol)
                   if xmol.HasSubstructMatch(core,useChirality=True):
                      mols.append(mol)
                      prop.append(' '.join(cells))
    return mols, prop

def ValidCheck(smi):  
    rr = {}
    rs = '.'.join(smi.split('|')[1:]).split('.')
    for ss in rs:
        cells = ss.split('[*:')
        if len(cells) > 2:
           map_num = {}
           for xx in cells[1:]:
               map_num[xx.split(']')[0]] = ''
           if len(map_num) > 1:
              if ss in rr: rr[ss] += 1
              else: rr[ss] = 1
    for r in rr:
        if rr[r] == 1:
           return False
    return True


#===========================================================================







