#!/usr/bin/env python

import argparse,textwrap
from src.nn_combo import nn_combo
from src.coding import CodingRgroupsDict
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def r_gen(cutoff, path):
    print()
    print('                                                                        AIR')
    print('===========================================================================')
    print()

    prm = {'cutoff':cutoff,'path':path,'ensemble':None,'ipt_dim':None,'coded_rgs':None,
           'mappedcore':None,'r_grps':None,'h_size':None,'r_radius':[]}

    with open(path+'/r.prm') as f:
         for line in f:
             if line.startswith('!'): continue
             if line.strip():
                cells = line.strip().split(':')
                if cells[0] == 'ensemble':
                   prm['ensemble'] = int(cells[1])
                elif cells[0] == 'input_dim':
                   prm['ipt_dim'] = int(cells[1])
                elif cells[0] == 'hidden_size':
                   prm['h_size'] = int(cells[1])
                elif cells[0] == 'mean':
                   prm['mean_r'] = float(cells[1])
                elif cells[0] == 'std':
                   prm['std_r'] = float(cells[1])
                elif cells[0] == 'core':
                   prm['mappedcore'] = line.strip().split('core: ')[1]
                elif cells[0] == 'r_radius':
                   rr = cells[1].split()
                   for r in rr:
                       prm['r_radius'].append(int(r)) 

    prm['r_grps'] = [{} for ii in range(len(prm['r_radius']))]
    with open(path+'/rgs.str') as f:
         for line in f:
             if line.startswith('!'): continue
             if line.strip():
                cells = line.strip().split('|')[1:]
                for ii in range(len(prm['r_radius'])):
                    prm['r_grps'][ii][cells[ii]] = 'R' + str(ii+1)
 
    num = 1
    for ii in range(len(prm['r_radius'])):
        num = num * len(prm['r_grps'][ii])
    prm['number'] = num

    prm['coded_rgs'] = CodingRgroupsDict(prm['r_grps'],prm['r_radius'])
 
    print('      Rough estimate of enumerated compounds: '+'{:,}'.format(num))
    print()
    print('      Enumeration with cutoff: ' + '%g'%cutoff) 
    print()

    nn_combo(**prm)

    print()
    print()
    print('      Enumeration written to r_gen.smi')
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

   parser.add_argument('-cutoff', action='store', dest='cutoff', type=float,
                        default = 7.0,
                        help=textwrap.dedent('cutoff to keep compounds\ndefault: 7.0'))

   parser.add_argument('-path', action='store', dest='path', default='./r_models',
                        help=textwrap.dedent('''path to the models\ndefault: ./r_models'''))

   arg_dict = vars(parser.parse_args())

   r_gen(**arg_dict)








   


