import iotbx
import numpy as np
import matplotlib.pyplot as plt
import random
import sys
sys.path.append("..")
from AEVclass import AEV, radial_aev_class
#compare CCC and generate a matrix

CCCC_pdb = '''REMARK  99 electronic Ligand Builder and Optimisation Workbench (eLBOW)
REMARK  99   - a module of PHENIX version dev-svn-
REMARK  99   - file written: Sat Apr 13 03:55:40 2019
REMARK  99
REMARK  99   Random seed: 3628800
REMARK  99   SMILES string: CCCC
FORMUL   1  LIG    C4 H10 
HETATM    1  C01 LIG A   1      -1.817  -0.000  -0.707  1.00 20.00      A    C  
HETATM    2  C02 LIG A   1      -0.290  -0.000  -0.707  1.00 20.00      A    C  
HETATM    3  C03 LIG A   1       0.290  -0.000   0.707  1.00 20.00      A    C  
HETATM    4  C04 LIG A   1       1.817   0.000   0.707  1.00 20.00      A    C  
HETATM    5 H011 LIG A   1      -2.178   0.884  -0.196  1.00 20.00      A    H  
HETATM    6 H012 LIG A   1      -2.178  -0.000  -1.728  1.00 20.00      A    H  
HETATM    7 H013 LIG A   1      -2.178  -0.885  -0.196  1.00 20.00      A    H  
HETATM    8 H021 LIG A   1       0.061  -0.882  -1.229  1.00 20.00      A    H  
HETATM    9 H022 LIG A   1       0.061   0.881  -1.229  1.00 20.00      A    H  
HETATM   10 H031 LIG A   1      -0.061   0.881   1.229  1.00 20.00      A    H  
HETATM   11 H032 LIG A   1      -0.061  -0.882   1.229  1.00 20.00      A    H  
HETATM   12 H041 LIG A   1       2.178  -0.884   0.196  1.00 20.00      A    H  
HETATM   13 H042 LIG A   1       2.178   0.000   1.728  1.00 20.00      A    H  
HETATM   14 H043 LIG A   1       2.178   0.885   0.196  1.00 20.00      A    H  
CONECT    1    2    5    6    7
CONECT    2    1    3    8    9
CONECT    3    2    4   10   11
CONECT    4    3   12   13   14
CONECT    5    1
CONECT    6    1
CONECT    7    1
CONECT    8    2
CONECT    9    2
CONECT   10    3
CONECT   11    3
CONECT   12    4
CONECT   13    4
CONECT   14    4
END'''

CCCS_pdb = '''REMARK  99 electronic Ligand Builder and Optimisation Workbench (eLBOW)
REMARK  99   - a module of PHENIX version dev-svn-
REMARK  99   - file written: Sat Apr 13 03:50:10 2019
REMARK  99
REMARK  99   Random seed: 3628800
REMARK  99   SMILES string: CCCS
FORMUL   1  LIG    C3 H8 S 
HETATM    1  C01 LIG A   1      -1.515   0.000  -0.678  1.00 20.00      A    C  
HETATM    2  C02 LIG A   1       0.012   0.000  -0.678  1.00 20.00      A    C  
HETATM    3  C03 LIG A   1       0.591   0.000   0.735  1.00 20.00      A    C  
HETATM    4  S04 LIG A   1       2.405  -0.000   0.652  1.00 20.00      A    S  
HETATM    5 H011 LIG A   1      -1.877   0.885  -0.168  1.00 20.00      A    H  
HETATM    6 H012 LIG A   1      -1.877   0.000  -1.700  1.00 20.00      A    H  
HETATM    7 H013 LIG A   1      -1.877  -0.885  -0.168  1.00 20.00      A    H  
HETATM    8 H021 LIG A   1       0.362   0.881  -1.201  1.00 20.00      A    H  
HETATM    9 H022 LIG A   1       0.362  -0.881  -1.201  1.00 20.00      A    H  
HETATM   10 H031 LIG A   1       0.254   0.884   1.262  1.00 20.00      A    H  
HETATM   11 H032 LIG A   1       0.254  -0.884   1.262  1.00 20.00      A    H  
HETATM   12 H041 LIG A   1       2.904  -0.000   1.882  1.00 20.00      A    H  
CONECT    1    2    5    6    7
CONECT    2    1    3    8    9
CONECT    3    2    4   10   11
CONECT    4    3   12
CONECT    5    1
CONECT    6    1
CONECT    7    1
CONECT    8    2
CONECT    9    2
CONECT   10    3
CONECT   11    3
CONECT   12    4
END'''


def _sort(k1, k2):
  if len(k1)==1:
    return 1
  else:
    return -1

#test 2 different element pdb file

def merge(a, b):
  for key, item in b.items():
    if key in a:
      for key1, item1 in item.items():
        if key1 not in a[key].keys():
          a[key].setdefault(key1, [])
          a[key][key1] = item1
    else:
      a.setdefault(key, {})
      a[key] = item
  return a

def randomcolor():
    colorArr = ['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
    color = ""
    for i in range(6):
        color += colorArr[random.randint(0, 14)]
    return "#"+color

def compare(pdb_file1=None, pdb_file2=None, record1=None, record2=None):
  if pdb_file1 and pdb_file2:
    aev1 = merge(AEV(pdb_file_name=pdb_file1).get_AEVS(), AEV(pdb_file_name=pdb_file2).get_items())
    aev2 = merge(AEV(pdb_file_name=pdb_file2).get_AEVS(), AEV(pdb_file_name=pdb_file1).get_items())
  else:
    aev1 = merge(AEV(raw_records=record1).get_AEVS(), AEV(raw_records=record2).get_items())
    aev2 = merge(AEV(raw_records=record2).get_AEVS(), AEV(raw_records=record1).get_items())
  list = ['C0', 'C1', 'C2']
  diff = {}
  for element in list:
    diff.setdefault(element, [])
    all1 = []
    for r_or_a, value in aev1[element].items():
      for v1 in value:
        all1.append(v1)
    for element2 in list:
      all2 = []
      for r_or_a2 in aev1[element].keys():
        value2 = aev2[element2][r_or_a2]
        for v2 in value2:
          all2.append(v2)
      covalue = np.corrcoef(all1, all2).tolist()
      diff[element].append(covalue[1][0])
  print diff
  return diff

def main(pdb_file_name1=None, pdb_file_name2=None):
  if pdb_file_name1 and pdb_file_name2:
    compare(pdb_file_name1, pdb_file_name2)
  else:
    compare(record1=CCCC_pdb, record2=CCCS_pdb)

if __name__ == '__main__':
    import os, sys
    main(*tuple(sys.argv[1:]))

