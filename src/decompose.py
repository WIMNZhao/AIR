import sys
from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition as rdRGD
from .utils import CoreFromSmarts, MolMatchSubsFromSmi


def Decomp(corefile,cpdfile):
    core = CoreFromSmarts(corefile)
    mols,prop = MolMatchSubsFromSmi(cpdfile,core)

    groups = []
    core_rgs = []
    uniq = {}

    #first alignment with unlabelled core
    res,unmatched = rdRGD.RGroupDecompose([core],mols,asSmiles=True)

    r_grps = [{} for ii in range(len(res[0])-1)]    

    for ii in range(len(res)):
        add = res[ii]['Core']
        for jj in range(len(r_grps)):
            key = 'R' + str(jj+1)
            hh  = '[H][*:' + str(jj+1) + ']'
            res[ii][key] = res[ii][key].replace('.'+hh,'').replace(hh+'.','')
            add = add + '|' + res[ii][key] 
        uniq[add] = ''
        core_rgs.append((add,prop[ii].split()[0:2]))

    groups.extend(res)
    mapped_core = res[0]['Core']

    #second alignment for symmetrical molecules
    core = Chem.MolFromSmarts(mapped_core)
    params = rdRGD.RGroupDecompositionParameters()
    params.onlyMatchAtRGroups = True
    for ii in range(len(mols)):
        res,unmatched = rdRGD.RGroupDecompose([core],[mols[ii]],asSmiles=True,options=params)
        for cpd in res:
           add = mapped_core
           for jj in range(len(r_grps)):
               key = 'R' + str(jj+1)
               if key not in cpd:
                  cpd[key]  = '[H][*:' + str(jj+1) + ']'
               add = add + '|' + cpd[key] 
           if add not in uniq:
              uniq[add] = ''
              core_rgs.append((add,prop[ii].split()[0:2]))
              groups.append(cpd)

    for cpd in groups:
        for ii in range(len(r_grps)):
            r_grps[ii][cpd['R'+str(ii+1)]] = 'R' + str(ii+1)

    rgs = open('rgs.smi','w')
    for cpd in core_rgs:
        rgs.write(' '.join(cpd[0].split('|')) + ' ' + ' '.join(cpd[1]) + '\r\n')


    return core_rgs, mapped_core, r_grps


