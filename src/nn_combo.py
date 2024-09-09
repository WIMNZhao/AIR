import torch,sys
import torch.utils.data as Data
import numpy as np
from .assemble import BatchAssemble
from .coding import CodingFromDict
from .utils import ValidCheck
from .nn import NNRG


def nn_combo(mappedcore, r_grps, coded_rgs, path, r_radius, ensemble, ipt_dim, h_size, cutoff, number, mean_r, std_r):
    ff = open('r_gen.smi','w')
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    ncpd = 0
    core_rgs  = []

    #Enumerate
    rs = len(r_grps)
    grps = [[r for r in rs] for rs in r_grps]

    idx = [0]*rs
    flag = True
    while flag:
          smi = mappedcore
          for ii in range(rs):
              smi = smi + '|' + grps[ii][idx[ii]]
          idx[-1] += 1
          for jj in range(rs-1,-1,-1):
              if idx[jj] == len(grps[jj]):
                 idx[jj] = 0
                 idx[jj-1] += 1
                 if idx[0] == len(grps[0]):
                    flag = False
                    break
          if ValidCheck(smi):
             core_rgs.append(smi)
             ncpd += 1

             if ncpd % 10000 == 0:
                ipts = CodingFromDict(core_rgs,coded_rgs)
                plist = np.empty(shape=(ensemble,len(ipts)))
                for ii in range(ensemble):
                    lenn = NNRG(ipt_dim,h_size,1)
                    lenn.to(device)
                    lenn.load_state_dict(torch.load(path+'/model_'+'%02d'%(ii+1)+'.ckpt'))
                    pred = lenn(ipts.to(device)).cpu().data.numpy()
                    plist[ii] = pred[:,0]*std_r + mean_r
                for ii in range(len(core_rgs)): 
                    mean = np.mean(plist[:,ii])
                    if mean >= cutoff:
                       smi = BatchAssemble([core_rgs[ii]])[0] 
                       if smi:
                          stdev = np.std(plist[:,ii])
                          ff.write(smi + ' %.2f'%mean + ' %.2f'%stdev + '\r\n')
                core_rgs = []
                sys.stdout.write('      Progress:   '+'%.1f%%\r'%(100.*(ncpd)/number))
                sys.stdout.flush()


    if core_rgs:
       ipts = CodingFromDict(core_rgs,coded_rgs)
       plist = np.empty(shape=(ensemble,len(ipts)))
       for ii in range(ensemble):
           lenn = NNRG(ipt_dim,h_size,1)
           lenn.to(device)
           lenn.load_state_dict(torch.load(path+'/model_'+'%02d'%(ii+1)+'.ckpt'))
           pred = lenn(ipts.to(device)).cpu().data.numpy()
           plist[ii] = pred[:,0]*std_r + mean_r
       for ii in range(len(core_rgs)): 
           mean = np.mean(plist[:,ii])
           if mean >= cutoff:
              smi = BatchAssemble([core_rgs[ii]])[0] 
              if smi:
                 stdev = np.std(plist[:,ii])
                 ff.write(smi + ' %.2f'%mean + ' %.2f'%stdev + '\r\n')
       sys.stdout.write('      Progress:   100%\r')
       sys.stdout.flush()







