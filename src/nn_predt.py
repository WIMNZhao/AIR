import torch
import torch.utils.data as Data
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition as rdRGD
from .assemble import BatchAssemble
from .coding import BatchCodingPred
from .nn import NNRG


def nn_predt(test, mappedcore, path, r_radius, ensemble, ipt_dim, mean_r, std_r, h_size, cutoff):
    ff = open('r_predt.smi','w')
    fr = open('r_raw.smi','w')
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    params = rdRGD.RGroupDecompositionParameters()
    params.onlyMatchAtRGroups = True
    core = Chem.MolFromSmarts(mappedcore)
    ngrps = 0
    for atom in core.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > ngrps:
           ngrps = map_num

    ncpd = 0
    core_rgs  = []
    mols = []

    with open(test) as f:
         for line in f:
             if line.startswith('!'): continue
             if line.strip():
                cells = line.split()
                mol = Chem.MolFromSmiles(cells[0])
                if mol:
                   xmol = Chem.AddHs(mol)
                   if xmol.HasSubstructMatch(core,useChirality=True):
                      res,unmatched = rdRGD.RGroupDecompose([core],[mol],asSmiles=True,options=params)
                      for cpd in res:
                          add = cpd['Core']
                          for jj in range(ngrps):
                              key = 'R' + str(jj+1)
                              if key not in cpd:
                                 cpd[key]  = '[H][*:' + str(jj+1) + ']'
                              add = add + '|' + cpd[key] 
                          core_rgs.append(add)
                          mols.append(' '.join(cells))
                          ncpd += 1

                          if ncpd % 10000 == 0:
                             ipts = BatchCodingPred(core_rgs,r_radius)
                             plist = np.empty(shape=(ensemble,len(ipts)))
                             for ii in range(ensemble):
                                 lenn = NNRG(ipt_dim,h_size,1)
                                 lenn.to(device)
                                 lenn.load_state_dict(torch.load(path+'/model_'+'%02d'%(ii+1)+'.ckpt'))
                                 pred = lenn(ipts.to(device)).cpu().data.numpy()
                                 plist[ii] = pred[:,0]*std_r + mean_r
                             for ii in range(len(mols)): 
                                 mean = np.mean(plist[:,ii])
                                 if mean >= cutoff:
                                    stdev = np.std(plist[:,ii])
                                    ff.write(mols[ii] + ' %.2f'%mean + ' %.2f'%stdev + '\r\n')
                                    fr.write(mols[ii])
                                    for jj in range(ensemble):
                                        fr.write(' %.2f'%plist[jj][ii])
                                    fr.write('\r\n')
                             core_rgs = []
                             mols = []

         if core_rgs:
            ipts = BatchCodingPred(core_rgs,r_radius)
            plist = np.empty(shape=(ensemble,len(ipts)))
            for ii in range(ensemble):
                lenn = NNRG(ipt_dim,h_size,1)
                lenn.to(device)
                lenn.load_state_dict(torch.load(path+'/model_'+'%02d'%(ii+1)+'.ckpt'))
                pred = lenn(ipts.to(device)).cpu().data.numpy()
                plist[ii] = pred[:,0]*std_r + mean_r
            for ii in range(len(mols)): 
                mean = np.mean(plist[:,ii])
                if mean >= cutoff:
                   stdev = np.std(plist[:,ii])
                   ff.write(mols[ii] + ' %.2f'%mean + ' %.2f'%stdev + '\r\n')
                   fr.write(mols[ii])
                   for jj in range(ensemble):
                       fr.write(' %.2f'%plist[jj][ii])
                   fr.write('\r\n')





