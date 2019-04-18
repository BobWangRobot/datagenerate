import iotbx
import numpy as np
import matplotlib.pyplot as plt
import random
import sys
sys.path.append("..")
from AEVclass import AEV, radial_aev_class
#compare CCC and generate a matrix
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
        color += colorArr[random.randint(0,14)]
    return "#"+color

def compare(pdb_file1, pdb_file2):
  aev1 = merge(AEV(pdb_file1).get_AEVS(), AEV(pdb_file2).get_items())
  aev2 = merge(AEV(pdb_file2).get_AEVS(), AEV(pdb_file1).get_items())
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
      for r_or_a2, value2 in aev2[element2].items():
        for v2 in value2:
          all2.append(v2)
      covalue = np.corrcoef(all1, all2).tolist()
      diff[element].append(covalue[1][0])
  print diff
  return diff




def main(pdb_file_name1, pdb_file_name2):
  compare(pdb_file_name1, pdb_file_name2)

if __name__ == '__main__':
    import os, sys
    main(*tuple(sys.argv[1:]))

