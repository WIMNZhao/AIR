import torch,sys
import torch.utils.data as Data
from torch.utils.data import TensorDataset
import numpy as np

from .nn import NNRG
from .coding import BatchCodingTrain


def nn_train(core_rgs, r_radius, ensemble=9, epochs=1000, lrate=0.001, h_size=16):
    inps,tgts = BatchCodingTrain(core_rgs,r_radius)

    mean = np.mean(tgts)
    std = np.std(tgts)
    tgts = (tgts - mean)/std
    tgts = torch.from_numpy(tgts).reshape(-1,1)

    ps = int((len(inps)+19)/20) * 100
    train_ds  = TensorDataset(inps,tgts)
 
    fit = [[] for ii in range(ensemble)]

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    ipt_dim = len(inps[0])

    for ii in range(ensemble):
        lenn = NNRG(ipt_dim,h_size,1)
        lenn.to(device)

        optimizer = torch.optim.Adam(lenn.parameters(), lr = lrate) 
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer,patience=ps,factor=0.9,min_lr=0.0001)
        loss_func = torch.nn.MSELoss()  

        for epoch in range(epochs):
            loader = Data.DataLoader(dataset=train_ds,batch_size=20,shuffle=True,drop_last=False)
            for x, y in loader:
                optimizer.zero_grad()   
                pred = lenn(x.to(device))     
                loss = loss_func(pred, y.to(device))    
                loss.backward() 
                optimizer.step()
                scheduler.step(loss)        

        sys.stdout.write('      Training in progress:   '+'%.1f%%\r'%(100.*(ii+1)/ensemble))
        sys.stdout.flush()
        torch.save(lenn.state_dict(),'r_models/model_'+'%02d'%(ii+1)+'.ckpt')
        pred = lenn(inps.to(device)).cpu().data.numpy()
        for p in pred:
            fit[ii].append(' %.3f'%(p[0]*std+mean))

    ff = open('r_train.smi','w')
    for ii in range(len(inps)):
        ff.write(' '.join(core_rgs[ii][1]))
        for jj in range(ensemble):
            ff.write(fit[jj][ii])
        ff.write('\r\n')

    return ipt_dim, mean, std








