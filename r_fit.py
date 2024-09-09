#!/usr/bin/env python

import argparse,textwrap,os
from src.decompose import Decomp
from src.maxradius import MaxRadiusList
from src.nn_train import nn_train
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def rgen_train(core,ref,ensemble=9,epochs=1000,lr=0.001,hsize=16):
    print()
    print('                                                                        AIR')
    print('===========================================================================')
    print()

    print('      R-group decomposition of the training set')
    core_rgs, core, r_grps = Decomp(core,ref)
    r_radius = MaxRadiusList(r_grps)

    # save 
    save_dir = 'r_models'
    if not os.path.isdir(save_dir):
       os.makedirs(save_dir)

    rgs = open(save_dir+'/rgs.str','w')
    for cpd in core_rgs:
        rgs.write(cpd[0] + '\r\n')
    #

    print()
    print('      Neural network training with number models: ' + '%d'%ensemble) 
    print('                                   number epochs: ' + '%g'%epochs) 
    print('                                   hidden neuron: ' + '%d'%hsize) 
    print('                                   learning rate: ' + '%g'%lr) 
    print()

    ipt_dim,mean,std = nn_train(core_rgs, r_radius, ensemble, epochs, lr, hsize)   

    prm = open(save_dir+'/r.prm','w')
    prm.write('input_dim: '+'%d'%ipt_dim+'\r\n')
    prm.write('hidden_size: '+'%d'%hsize+'\r\n')
    prm.write('epochs: '+'%d'%epochs+'\r\n')
    prm.write('lr: '+'%g'%lr+'\r\n')
    prm.write('mean: '+'%f'%mean+'\r\n')
    prm.write('std: '+'%f'%std+'\r\n')
    prm.write('ensemble: '+'%d'%ensemble+'\r\n')
    prm.write('core: ' + core +'\r\n')
    prm.write('r_radius: '+' '.join(str(radius) for radius in r_radius)+'\r\n')

    print()
    print()
    print('      Fitted data written to r_train.smi')
    print('      R-group decomposition written to rgs.smi')
    print('      Models and parameters in r_models')
    print()
    print('===========================================================================')
    print()
    print('                              ------End------')
    print()


if __name__ == "__main__":
   parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                    description=textwrap.dedent('''\
                                    AI-augmented R-group exploration in medicinal chemistry
                                    Hongtao Zhao'''))     

   parser.add_argument('-core', action='store', dest='core', default='core.smi',
                        help=textwrap.dedent('''path to the SMARTS defined core\ndefault: core.smi'''))

   parser.add_argument('-ref', action='store', dest='ref', default='ref.smi',
                        help=textwrap.dedent('''path to the reference smiles and activity\ndefault: ref.smi'''))

   parser.add_argument('-epochs', action='store', dest='epochs', type=int,
                        default = 200,
                        help=textwrap.dedent('''training epochs\ndefault: 200'''))

   parser.add_argument('-lr', action='store', dest='lr', type=float,
                        default = 0.001,
                        help=textwrap.dedent('''learning rate\ndefault: 0.001'''))

   parser.add_argument('-hsize', action='store', dest='hsize', type=int,
                        default = 16,
                        help=textwrap.dedent('''hidden size\ndefault: 16'''))

   parser.add_argument('-ensemble', action='store', dest='ensemble', type=int,
                        default = 9,
                        help=textwrap.dedent('number of models to evaluate uncertainty\ndefault: 9'))

   arg_dict = vars(parser.parse_args())

   rgen_train(**arg_dict)








   








   


