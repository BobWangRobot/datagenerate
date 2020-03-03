from iotbx import pdb
import numpy as np
import matplotlib.pyplot as plt
import random
import iotbx
import math
import collections
import os, sys
import time
import mmtbx

sys.path.append("..")
from AEVclass2 import *


def main(filename1=None, filename2=None):
  if filename1 and filename2:
    y1 = []
    y2 = []
    name = []
    a = AEV(pdb_file_name=filename1)
    b = AEV(pdb_file_name=filename2)
    a.five = next(a.generate_ca())
    b.five = next(b.generate_ca())
    a.Rpart()
    b.Rpart()
    aev1 = a.AEVs
    aev2 = b.AEVs
    for ele, values in aev1.items():
      for r_or_a, value in values.items():
        for list1 in value:
          y1.append(list1)
    for ele, values in aev2.items():
      for r_or_a, value in values.items():
        for list2 in value:
          y2.append(list2)
    for i in range(5):
      y11 = y1[0:8]
      y22 = y2[0:8]
      x = range(len(y11))
      for v in a.rdistance['C1']['C']:
        name.append('%0.3f' % v)
      plt.title("C%s AEV values of %s and %s  " % (i+1, filename1.replace('.pdb', ''), (filename2.replace('.pdb', ''))))
      plt.xlabel("Rs_values")
      plt.ylabel("value")
      plt.plot(x, y11, 'r', label='%s' % filename1.replace('.pdb', ''))
      plt.plot(x, y22, 'g', label='%s' % filename2.replace('.pdb', ''))
      plt.xticks(x[::1], a.rs_values)
      plt.legend()
      plt.savefig('./difference/C%s %s.jpg' % (i+1, filename2.replace('.pdb', '%s' % ele)))
      plt.show()
      del y1[0:8]
      del y2[0:8]
    else:
      a = AEV(pdb_file_name=filename1)
      print(a.find_function())



if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))