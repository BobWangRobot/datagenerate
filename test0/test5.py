import iotbx
import numpy as np
import matplotlib.pyplot as plt
import random
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



#test 2 different element pdb file



# def main(pdb_file_name1=None, pdb_file_name2=None, elelist=None):
#   if pdb_file_name1 and pdb_file_name2:
#     aev1 = AEV(pdb_file_name=pdb_file_name1)
#     aev2 = AEV(pdb_file_name=pdb_file_name2)
#   else:
#     aev1 = AEV(raw_records=CCCC_pdb)
#     aev2 = AEV(raw_records=CCCS_pdb)
#   print(aev1.Rpart(), aev2.Apart())
#   #print(aev1.compare(aev2,elelist))
def main(pdb_file_name):
  aev = AEV(pdb_file_name).Rpart()
  print(aev)
if __name__ == '__main__':
    import os, sys
    main(*tuple(sys.argv[1:]))

