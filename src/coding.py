import torch
import numpy as np
from .fingerprint import FingerPrint

def Coding(core_rgps,radius):
    cells = core_rgps.split('|')[1:]
    fp = []
    for ii in range(len(cells)):
        fp.extend(FingerPrint(cells[ii],ii+1,radius[ii]))

    return fp

def BatchCodingTrain(core_rgps,radius):
    ipt = []
    pic50 = []
    for cpd in core_rgps:
        ipt.append(Coding(cpd[0],radius))
        pic50.append(float(cpd[1][1]))
    return torch.from_numpy(np.array(ipt)), np.array(pic50,dtype=np.float32) 
    #return torch.from_numpy(np.array(ipt)), torch.from_numpy(np.array(pic50,dtype=np.float32).reshape(-1,1)) 

def BatchCodingPred(core_rgps,radius):
    ipt = []
    for cpd in core_rgps:
        ipt.append(Coding(cpd,radius))
    return torch.from_numpy(np.array(ipt))

def CodingRgroupsDict(r_grps,r_radius):
    coded = {}
    for ii, rs in enumerate(r_grps):
        for rr in rs:
            coded[rr+'|'+rs[rr]] = FingerPrint(rr,ii+1,r_radius[ii])
    return coded

def CodingFromDict(core_rgps,coded):
    ipt = []
    for cpd in core_rgps:
        cells = cpd.split('|')[1:]
        fp = []
        for ii in range(len(cells)):
            fp.extend(coded[cells[ii]+'|R'+str(ii+1)]) 
        ipt.append(fp)
    return torch.from_numpy(np.array(ipt))
    







