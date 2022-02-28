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
    #a.Rpart()
    #b.Rpart()
    aev1 = a.Rpart('all')
    aev2 = b.Rpart('all')
    print(aev1,aev2)
    for ele, values in aev1.items():
      for r_or_a, value in values.items():
        for list1 in value:
          y1.append(list1)
      name.append(ele)
    for ele, values in aev2.items():
      for r_or_a, value in values.items():
        for list2 in value:
          y2.append(list2)
    x = range(len(y1))
    plt.title("%s and %s radial part AEV values " % (filename1.replace('.pdb', ''), (filename2.replace('.pdb', ''))))
    plt.xlabel("atom number")
    plt.ylabel("value")
    plt.plot(x, y1, 'r', label='%s' % filename1.replace('.pdb', ''))
    plt.plot(x, y2, 'g', label='%s' % filename2.replace('.pdb', ''))
    plt.xticks(x[::8], name)
    plt.legend()
    # plt.savefig('./difference/%s.jpg' % filename2.replace('.pdb', '%s' % ele))
    plt.show()
  else:
    a = AEV(pdb_file_name=filename1)
    print(a.find_function())
    
if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))
