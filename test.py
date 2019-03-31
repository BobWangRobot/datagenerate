import iotbx
from AEVclass import AEV
import numpy as np

# def compare(a,b):
#   diff = {}
#   for atom1, list1 in a.items():
#     diff.setdefault(atom1, [])
#     for atom2, list2 in b.items():
#       if atom1 == atom2:
#         diff[atom1].append(np.corrcoef(list1, list2))
#   return diff

def main(pdb_file_name1):
  a = AEV(pdb_file_name1)
  print(a.Rpart())
  print(a.Apart())
  print(a.get_AEVS())




if __name__ == '__main__':
    import os, sys
    main(*tuple(sys.argv[1:]))

