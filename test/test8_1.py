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
#plot all atomes compare result
sys.path.append("..")
from AEVclass3 import *


def main(filename1=None, filename2=None):
  if filename1 and filename2:
    a = AEV(pdb_file_name=filename1)
    b = AEV(pdb_file_name=filename2)
    for a.five in a.generate_ca():
      a.Rpart()
    for b.five in b.generate_ca():
      b.Rpart()
    y1 = a.AEVs
    y2 = b.AEVs
    x = range(len(y1))
    plt.title("%s and %s radial part AEV values " % (filename1.replace('.pdb', ''), (filename2.replace('.pdb', ''))))
    plt.xlabel("atom number")
    plt.ylabel("value")
    plt.plot(x, y1, 'r', label='%s' % filename1.replace('.pdb', ''))
    plt.plot(x, y2, 'g', label='%s' % filename2.replace('.pdb', ''))
    plt.xticks(x[::8], (a.e_name+b.e_name))
    plt.legend()
    plt.savefig('./difference/%s.jpg' % (filename1.replace('.pdb','')+filename2.replace('.pdb','')))
    plt.show()
  else:
    a = AEV(pdb_file_name=filename1)
    print(a.find_function())


if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))