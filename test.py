import iotbx
from AEVclass import AEV
import numpy as np

def compare(a,b):
  diff = {}
  for atom1, list1 in a.items():
    diff.setdefault(atom1, [])
    for atom2, list2 in b.items():
      if atom1 == atom2:
        diff[atom1].append(np.corrcoef(list1, list2))
  return diff

def main(pdb_file_name1,pdb_file_name2):
    a = AEV(pdb_file_name1)
    b = AEV(pdb_file_name2)
    print(type(a,b))
    f = compare(a, b)
    print(f)
    #print(a.Rpart())
    #print(a.Apart())
    print(a.get_AEVS())




if __name__ == '__main__':
    import os, sys
    main(*tuple(sys.argv[1:]))
#
# #a = AEV('acid060.pdb')
# #print(a.Atome_classify())
# #b = Rpart('acid060.pdb')
# #print(b.R_AEV())
# #c = Apart('acid060.pdb')
# #print(c.Aeq())
#
# a = AEV(pdb_file_name='acid060.pdb')
# b = AEV(pdb_file_name='acid300.pdb')
# def compare(dic1, dic2):

