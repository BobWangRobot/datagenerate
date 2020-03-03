from iotbx import pdb
import numpy as np
import iotbx
import math
import collections
import os, sys
import time
import mmtbx
sys.path.append("..")
from AEVclass2 import *


def main(filename1=None, filename2=None):
  diff = {}
  all = 0
  if filename1 and filename2:
    a = AEV(pdb_file_name=filename1)
    b = AEV(pdb_file_name=filename2)
    for a.five in a.generate_ca():
      a.Rpart()
    for b.five in b.generate_ca():
      b.Rpart()
  for aev1, aev2 in zip(a.Rpart().items(), b.Rpart().items()):
    ele1 = aev1[0]
    ele2 = aev2[0]
    list1 = aev1[1].values()
    list2 = aev2[1].values()
    covalue = np.corrcoef(list1, list2).tolist()
    a.diffs.setdefault(ele1 + ele2, covalue[1][0])
    all += covalue[1][0]
    # if covalue[1][0] > limit:
    a.diffs.setdefault(ele1 + ele2, covalue[1][0])
  a.diffs.setdefault('all', all / 5)
  print (a.diffs)
  return a.diffs

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))