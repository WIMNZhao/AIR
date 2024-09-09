from rdkit import Chem

out = open('r_gen_new.smi','w')

made = {}

with open('ref.smi') as f:
     for line in f.readlines():
         cells = line.split()
         mol = Chem.MolFromSmiles(cells[0])
         if mol is not None:
            smi = Chem.MolToSmiles(mol) #, isomericSmiles=False)
            made[smi] = ''

with open('test.smi') as f:
     for line in f.readlines():
         cells = line.split()
         mol = Chem.MolFromSmiles(cells[0])
         if mol is not None:
            smi = Chem.MolToSmiles(mol) #, isomericSmiles=False)
            made[smi] = ''

with open('r_gen.smi') as f:
     for line in f.readlines():
         cells = line.split()
         mol = Chem.MolFromSmiles(cells[0])
         if mol is not None:
            smi = Chem.MolToSmiles(mol) #, isomericSmiles=False)
            if smi not in made:
               out.write(line)


