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
    a = a.get_AEVS()
    b = b.get_AEVS()
    # print(a.get_AEVS(), b.get_AEVS())
    f = compare(a, b)
    print(f)


if __name__ == '__main__':
    import os, sys
    main(*tuple(sys.argv[1:]))
