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
import xlwt
sys.path.append("..")
from AEVclass3 import *

def plot(AEV1, AEV2=None):
  font1 = {'family': 'Times New Roman',
           'weight': 'normal',
           'size': 20,
           }
  font2 = {'family': 'Times New Roman',
           'weight': 'normal',
           'size': 18,
           }
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
    plt.title("AEV values figure",font1)
    plt.xlabel("atom index",font1)
    plt.ylabel("AEV values",font1)
    plt.plot(x, list1, 'orange', label='reference_structure' )
    plt.plot(x, list2, 'green', label='real_structure',linestyle=':')
    plt.xticks(x[::16], (range(1, 30)))
    #plt.xticks(x[::16], (a_name + b_name for a_name, b_name in zip(name1, name2)))
    plt.legend(prop=font2,loc='best')
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
    plt.title("reference structure AEV values ", font1)
    plt.xlabel("atom index",font1)
    plt.ylabel("AEV values",font1)
    plt.plot(x, list1, 'green', label='reference structure')
    # plt.plot(x, list1, 'green', label='real structure')
    plt.legend(prop=font2,loc='best')
    plt.xticks(x[::16], (range(1,23)))
#small plot
    # plt.axes([0.17, 0.42, 0.1, 0.2])
    # y1 = list1[:8]
    # x1 = range(len(y1))
    # plt.plot(x1, y1, 'green')
    # plt.xticks(x1,('Rs%r'%i for i in range(1,9)))
    # plt.axes([0.25, 0.66, 0.1, 0.2])
    # y2 = list1[8:16]
    # x2 = range(len(y2))
    # plt.plot(x2, y2, 'green')
    # plt.xticks(x2, ('Rs%r' % i for i in range(1, 9)))
    # plt.xticks(x[::16], (a_name for a_name in name1))
    plt.show()

def data_save(datas):
  f = xlwt.Workbook()
  sheet1 = f.add_sheet(u'sheet1', cell_overwrite_ok=True)
  i = -1
  j = 1
  for key1,dict1 in datas.items():
    sheet1.write(i+2, 0, key1)
    i += 1
    sheet1.write(0, j+1, 'Rs%s'%j)
    i += 1
    j += 1
  i = 1
  for dict2 in datas.values():
    for list1 in dict2.values():
      j = 2
      for value in list1:
        value = float('%.2f'%value)
        sheet1.write(i, j, value)
        j = j + 1
      i = i + 1
  f.save('helix_sheet.xls')

def main(direction, scope, AEV1=None, AEV2=None):
  if AEV1 and AEV2:
    a = AEV(direction,scope, pdb_file_name=AEV1)
    b = AEV(direction,scope, pdb_file_name=AEV2)
    for a.five in a.generate_ca():
      a.get_AEVs()
      # a.Rpart()
    for b.five in b.generate_ca():
      b.get_AEVs()
      # b.Rpart()
    print(a.AEVs, b.AEVs)
    plot(a.AEVs, b.AEVs)
    # data_save(a.AEVs)
  elif AEV1:
    a = AEV(direction, scope, pdb_file_name=AEV1)
    for a.five in a.generate_ca():
      a.get_AEVs()
    plot(a.AEVs)
    data_save(a.AEVs)


if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))

