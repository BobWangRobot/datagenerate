from iotbx import pdb
import numpy as np
import matplotlib.pyplot as plt
import iotbx
import math
import collections
import os, sys
import time
import mmtbx
import xlwt
from AEVclass4 import *

perfect_helix = {
  'B1':[0.0000, 0.5487, 0.4894, 0.4189, 0.1576, 0.0031, 0.0000, 0.0000, 7.4866, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
  'B2':[0.0000, 0.5487, 0.4893, 0.4190, 0.1577, 0.0031, 0.0000, 0.0000, 7.4751, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
  'M':[0.0000, 1.0957, 0.4034, 0.4808, 0.0441, 0.0000, 0.0000, 0.0000, 16.3600, 4.8532, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
  'E1':[0.0000, 0.5488, 0.4894, 0.4191, 0.1577, 0.0031, 0.0000, 0.0000, 3.8169, 19.6605, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
  'E2':[0.0000, 0.5487, 0.4895, 0.4190, 0.1576, 0.0031, 0.0000, 0.0000, 3.8187, 19.6749, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000]
}


def sort_dic():
  perfect_helix_12 = collections.OrderedDict()
  for i in ["B1", "B2", "M", "E2", "E1"]:
    perfect_helix_12[i] = perfect_helix[i]
  return perfect_helix_12

# def compare(AEV):
#   diffs = diff_class()
#   perfect_helix_12 = sort_dic()
#   for ele1, item1 in AEV.items():
#       diffs.setdefault(ele1, OrderedDict())
#       for ele2, item2 in perfect_helix_12.items():
#         com_list1 = []
#         com_list2 = []
#         for list1, list2 in zip(item1.values(), item2.values()):
#           com_list1.extend(list1)
#           com_list2.extend(list2)
#         covalue = np.corrcoef(com_list1, com_list2).tolist()
#         diffs[ele1].setdefault(ele2, covalue[1][0])
#   return diffs
def compare(data):
  result = diff_class()
  perfect_helix = sort_dic()
  for key,value in perfect_helix.items():
    if key == 'B1':
      for key1, value1 in data.BAEVs.items():
        covalue = np.corrcoef(value, value1).tolist()
        result.setdefault(key1, OrderedDict())
        result[key1].setdefault(key, covalue[1][0])
    elif key == 'B2':
      for key1, value1 in data.BAEVs.items():
        covalue = np.corrcoef(value, value1).tolist()
        print(covalue[1][0])
        result[key1].setdefault(key, covalue[1][0])
    elif key == 'M':
      for key1, value1 in data.MAEVs.items():
        covalue = np.corrcoef(value, value1).tolist()
        try:
          result[key1].setdefault(key, covalue[1][0])
        except KeyError:
          result.setdefault(key1, OrderedDict())
          result[key1].setdefault(key, covalue[1][0])
    elif key == 'E2':
      for key1, value1 in data.EAEVs.items():
        covalue = np.corrcoef(value, value1).tolist()
        try:
          result[key1].setdefault(key, covalue[1][0])
        except KeyError:
          result.setdefault(key1, OrderedDict())
          result[key1].setdefault(key, covalue[1][0])
    elif key == 'E1':
      for key1, value1 in data.EAEVs.items():
        covalue = np.corrcoef(value, value1).tolist()
        try:
          result[key1].setdefault(key, covalue[1][0])
        except KeyError:
          result.setdefault(key1, OrderedDict())
          result[key1].setdefault(key, covalue[1][0])
  return result

def data_save(datas,name):
  f = xlwt.Workbook()
  sheet1 = f.add_sheet(u'sheet1', cell_overwrite_ok=True)
  line_list = []
  row_list = datas.keys()
  #write header
  for value in datas['%s'%row_list[0]].keys():
    line_list.append(value)
  for i in range(len(row_list)):
    sheet1.write(i+1, 0, row_list[i])
  for j in range(len(line_list)):
    sheet1.write(0, j+1, line_list[j])
  #write values
  i = 1
  all_num = 0
  for dict2 in datas.values():
    j = 1
    max_num = 0
    for value in dict2.values():
      if value > max_num:
        max_num = value
      value = float('%.4f' % value)
      sheet1.write(i, j, value)
      j = j + 1
    all_num += max_num
    # print(max_num)
    i = i + 1
  average = all_num / (i - 1)
  print(i)
  f.save('%s(%0.6f).xls' % (name, average))

def main(filename):
  t0 = time.time()
  a = AEV(pdb_file_name=filename)
  a.generate_AEV()
  # print(a.BAEVs, a.MAEVs, a.EAEVs)
  print(compare(a))
  print('time', time.time()-t0)
  # data_save(com_result,name=filename.replace('.pdb',''))



if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))


