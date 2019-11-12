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
from AEVclass3 import *

def plot(AEV1, AEV2=None):
  if AEV1 and AEV2:
    list1 = []
    list2 = []
    name1 = []
    name2 = []
    for element1, value1 in AEV1.items():
      name1.append(element1)
      for item1 in value1.values():
        list1 += item1
    for element2, value2 in AEV2.items():
      name2.append(element2)
      for item2 in value2.values():
        list2 += item2
    x = range(len(list1))
    plt.title("helix1 and helix2 AEV values ")
    plt.xlabel("atom number")
    plt.ylabel("AEV values")
    plt.plot(x, list1, 'orange', label='helix1' )
    plt.plot(x, list2, 'green', label='helix2',linestyle=':')
    plt.xticks(x[::16], (a_name + b_name for a_name, b_name in zip(name1, name2)))
    plt.legend()
    plt.show()
    #plt.savefig('./difference/1.jpg)
  elif AEV1:
    list1 = []
    name1 = []
    for element1, value1 in AEV1.items():
      name1.append(element1)
      for item1 in value1.values():
        list1 += item1
    x = range(len(list1))
    plt.title("helix1 AEV values ")
    plt.xlabel("atom number")
    plt.ylabel("AEV values")
    plt.plot(x, list1, 'green', label='perfect_helix')
    plt.xticks(x[::16], (a_name for a_name in name1))
    plt.legend()
    plt.show()

def main(direction, scope, AEV1=None, AEV2=None):
  if AEV1 and AEV2:
    a = AEV(direction,scope, pdb_file_name=AEV1)
    b = AEV(direction,scope, pdb_file_name=AEV2)
    for a.five in a.generate_ca():
      a.get_AEVs()
    for b.five in b.generate_ca():
      b.get_AEVs()
    print(a.AEVs, b.AEVs)
    plot(a.AEVs, b.AEVs)
  elif AEV1:
    a = AEV(direction, scope, pdb_file_name=AEV1)
    for a.five in a.generate_ca():
      a.get_AEVs()
    plot(a.AEVs)

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))

