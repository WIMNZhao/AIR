#!/usr/bin/env python

import argparse,textwrap
from src.nn_predt import nn_predt
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def r_gen(test,cutoff,path):
    print()
    print('                                                                        AIR')
    print('===========================================================================')
    print()

    prm = {'test':test,'cutoff':cutoff,'path':path,'ensemble':None,'ipt_dim':None,
           'mean_r':None, 'std_r':None, 'h_size':None,'mappedcore':None,'r_radius':[]}

    with open(path+'/r.prm') as f:
         for line in f:
             if line.startswith('!'): continue
             if line.strip():
                cells = line.strip().split(':')
                if cells[0] == 'ensemble':
                   prm['ensemble'] = int(cells[1])
                elif cells[0] == 'input_dim':
                   prm['ipt_dim'] = int(cells[1])
                elif cells[0] == 'mean':
                   prm['mean_r'] = float(cells[1])
                elif cells[0] == 'std':
                   prm['std_r'] = float(cells[1])
                elif cells[0] == 'hidden_size':
                   prm['h_size'] = int(cells[1])
                elif cells[0] == 'core':
                   prm['mappedcore'] = line.strip().split('core: ')[1]
                elif cells[0] == 'r_radius':
                   rr = cells[1].split()
                   for r in rr:
                       prm['r_radius'].append(int(r))   

    print('      Prediction with cutoff: ' + '%g'%cutoff) 
    nn_predt(**prm)

    print()
    print('      Prediction written to r_predt.smi')
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

   parser.add_argument('-f', action='store', dest='test', default='test.smi',
                        help=textwrap.dedent('''path to the test smiles\ndefault: test.smi'''))

   parser.add_argument('-cutoff', action='store', dest='cutoff', type=float,
                        default = -999,
                        help=textwrap.dedent('''cutoff to keep compounds\ndefault: -999'''))

   parser.add_argument('-path', action='store', dest='path', default='./r_models',
                        help=textwrap.dedent('''path to the models\ndefault: ./r_models'''))

   arg_dict = vars(parser.parse_args())

   r_gen(**arg_dict)








   








   


