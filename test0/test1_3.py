import iotbx
import numpy as np
import matplotlib.pyplot as plt
import random
from AEVclass import *

def _sort(k1, k2):
  if len(k1)==1:
    return 1
  else:
    return -1

#test 2 different pdb-files and plot


def randomcolor():
    colorArr = ['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
    color = ""
    for i in range(6):
        color += colorArr[random.randint(0,14)]
    return "#"+color

def plotvalue(pdb1, pdb2, elename):
  y1 = []
  y2 = []
  name = []
  aev1 = AEV(pdb1).merge(AEV(pdb2).get_items())
  aev2 = AEV(pdb2).merge(AEV(pdb1).get_items())
  print(aev1, aev2)
  for ele, values in aev1.items():
    if elename in ele:
      for r_or_a, value in values.items():
        for list1 in value:
          y1.append(list1)
        name.append(r_or_a)
  for element, values in aev2.items():
    if elename in element:
      for r_or_a in name:
        value = values[r_or_a]
        for list2 in value:
          y2.append(list2)
  x = range(len(y1))
  plt.title("AEV of C element between %s and %s" % (pdb1.replace('.pdb',''), pdb2.replace('.pdb','')))
  plt.xlabel("radial or angular of %s atom" % elename)
  plt.ylabel("AEV value")
  plt.plot(x, y1, 'r', label='%s' % pdb1.replace('.pdb', ''))
  plt.plot(x, y2, 'g', label='%s' % pdb2.replace('.pdb', ''))
  plt.xticks(x[::8], name)
  plt.legend()
  # plt.savefig('./difference/%s.jpg' % pdb2.replace('.pdb','%s'%elename))
  plt.show()

def plotcompare(diff,ele,pdb_filename):
  name = []
  y1 = []
  for element, values in diff.items():
    if ele in element:
      for r_or_a, value in values.items():
        name.append(r_or_a)
        y1.append(value[0])
  x = range(len(y1))
  plt.title("Correlation coefficient of %s atom"%ele)
  plt.xlabel("radial or angular to %s"%ele)
  plt.ylabel("value")
  plt.plot(x, y1, 'ob')
  plt.xticks(x, name)
  plt.savefig('./corrcoee/%s.jpg' % pdb_filename.replace('.pdb','%s'%ele))
  plt.show()

def main(pdb_file_name1, pdb_file_name2, elename):
  # diff = compare(pdb_file_name1, pdb_file_name2)
  # print(diff)
  # plotcompare(diff, elename, pdb_file_name1)
  plotvalue(pdb_file_name1, pdb_file_name2, elename)



if __name__ == '__main__':
    import os, sys
    main(*tuple(sys.argv[1:]))
import sys

