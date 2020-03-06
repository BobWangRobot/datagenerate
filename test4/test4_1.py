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
from AEVclass4 import *

perfect_helix = {
  'B1': {'CA': [1.2350877996290163e-06, 0.5487007234202762, 0.4893913937038461, 0.418942571898836, 0.15761683314261377,
                0.0030591217473060262, 3.341899561282178e-13, 3.2577683449187475e-29],
         'CACA': [0.05848912440308649, 0.18385850960221692, 0.06695068795914315, 0.0015606036396054806,
                  0.6146490985751465, 0.7199060452933045, 0.08600518255430749, 0.0007865537569315071,
                  0.4672260473848402, 0.4925547568981271, 0.04695075906892053, 0.00022117411047351224,
                  0.0599557279392109, 0.0785639219508055, 0.009300905815509148, 4.343254929231361e-05]},
  'B2': {'CA': [2.468386914195106e-06, 1.0964924447336064, 0.48948905102281054, 0.419009313808839, 0.15770496142721793,
                0.003067782562540799, 3.368584375149749e-13, 3.298566244146484e-29],
         'CACA': [0.1425640639458447, 1.2100152577462586, 1.0794338630583522, 0.0759922334544158, 0.668743132321037,
                  0.9527091974037247, 0.3888026731123253, 0.09077633474373654, 0.4762487599258398, 0.5378655945988641,
                  0.1375704818038189, 0.030731534894061725, 0.06005858429856836, 0.07876374838171123,
                  0.009728469863022782, 0.000196620432088233]},
  'M': {'CA': [2.4763915687013514e-06, 1.0974629285678776, 0.9789480882731942, 0.8382561256541073, 0.31526093727861015,
               0.006118653038005154, 6.685142725857298e-13, 6.51776314019727e-29],
        'CACA': [0.22807578336727607, 1.564283596277074, 1.2567864976416725, 0.08471677967053275, 1.3986393894577918,
                 2.3649373013530255, 1.880130400156958, 0.9106602548484867, 1.0238600715441326, 1.5511137419048662,
                 1.3783657046202034, 0.7405290062464193, 0.12659794010545605, 0.2138684775596562, 0.20721852219687614,
                 0.16595733851456881]},

'E1': {'CA': [2.4817862217911787e-06, 1.096645326043909, 0.4897576948693963, 0.41900810958387413, 0.15759246799832402,
                0.003059304919665151, 3.3425613764519415e-13, 3.2587779959841476e-29],
         'CACA': [0.14267331692834523, 1.2105590864961961, 1.0805470641589316, 0.07614789578699133, 0.6687070017322899,
                  0.952357352408722, 0.38855964811273935, 0.09086094756334955, 0.4761761280327183, 0.5378333667316271,
                  0.13742980195741467, 0.030727359496837303, 0.059959948945917954, 0.0786678130700521,
                  0.009719950682316056, 0.00019620765952365482]},
  'E2': {
    'CA': [1.2447365820713111e-06, 0.5488012464794341, 0.48941376612736126, 0.4191044892239558, 0.15770932332511725,
           0.00307227728922915, 3.382645769684755e-13, 3.320139206216676e-29],
    'CACA': [0.05852159362395412, 0.18374541933360714, 0.06683507231438032, 0.0015558587127007689, 0.6141657695061328,
             0.7191866402720278, 0.08590175957667005, 0.0007851974366386065, 0.46713569769991814, 0.49222127637650337,
             0.04689088136132048, 0.0002207246047364507, 0.060088854556646136, 0.07867669339194452,
             0.009307195462884547, 4.34233482508004e-05]}}

def sort_dic():
  perfect_helix_12 = collections.OrderedDict()
  for i in ["B1", "B2", "M", "E1", "E2"]:
    perfect_helix_12[i] = perfect_helix[i]
  return perfect_helix_12

def compare(AEV):
  diffs = diff_class()
  perfect_helix_12 = sort_dic()
  for ele1, item1 in AEV.items():
      diffs.setdefault(ele1, OrderedDict())
      for ele2, item2 in perfect_helix_12.items():
        com_list1 = []
        com_list2 = []
        for list1, list2 in zip(item1.values(), item2.values()):
          com_list1.extend(list1)
          com_list2.extend(list2)
        covalue = np.corrcoef(com_list1, com_list2).tolist()
        diffs[ele1].setdefault(ele2, covalue[1][0])
  return diffs

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

def main(direction, scope, filename):
  a = AEV(direction,scope,pdb_file_name=filename)
  for a.five in a.generate_ca():
    # a.get_AEVs()
    a.Rpart()
  com_result = compare(a.AEVs)
  print(com_result)
  # data_save(com_result,name=filename.replace('.pdb',''))



if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))


