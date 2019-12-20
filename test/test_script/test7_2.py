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


def main(filename1):
  if filename1 :
    y1 = []
    name = []
    a = AEV(pdb_file_name=filename1)
    a.five = next(a.generate_ca())
    aev1 = a.Rpart()
    for ele, values in aev1.items():
      for r_or_a, value in values.items():
        for list1 in value:
          y1.append(list1)
      name.append(ele)
    x = range(len(y1))
    plt.title("%s radial part AEV values " % (filename1.replace('.pdb', '')))
    plt.xlabel("atom-number")
    plt.ylabel("value")
    plt.plot(x, y1, 'r', label='%s' % filename1.replace('.pdb', ''))
    plt.xticks(x[::8], name)
    plt.legend()
    plt.savefig('./difference/%s.jpg' % filename1.replace('.pdb', '%s' % ele))
    plt.show()
  else:
    pass


if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))
