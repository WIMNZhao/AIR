import math,random

ref = open("ref.smi",'w')
test = open("test.smi",'w')

cpd = {}

with open("all.smi") as f:
      for line in f:
          if line.startswith('SMILES'): continue
          if len(line.split())>=2:
             cells = line.split()
             smi = cells[0]
             pki = float(cells[1])
             cpd[smi] = pki

l = list(cpd.items())
random.shuffle(l)

ncpd = len(l)
nref = 0.8 * ncpd

for i in range(ncpd):
    if i < nref:
       ref.write(l[i][0]+' %.2f'%l[i][1]+'\n')
    else:
       test.write(l[i][0]+' %.2f'%l[i][1]+'\n')
     
          
