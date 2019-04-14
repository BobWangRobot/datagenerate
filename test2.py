import iotbx
from AEVclass import AEV, radial_aev_class
import numpy as np
import matplotlib.pyplot as plt
import random

def _sort(k1, k2):
  if len(k1)==1:
    return 1
  else:
    return -1



def randomcolor():
    colorArr = ['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
    color = ""
    for i in range(6):
        color += colorArr[random.randint(0,14)]
    return "#"+color

def compare(aev1, aev2):
  diff = aev2.get_AEVS()
  for element, values in aev1.get_AEVS().items():
    for r_or_a, value in values.items():
      if 0 in value:
        diff[element][r_or_a] = [1]
      else:
        covalue = np.corrcoef(value, diff[element][r_or_a]).tolist()
        diff[element][r_or_a] = [covalue[1][0]]
  return diff

def compare_all(origin, march, number):
  a = AEV(origin)
  b = AEV(march)
  total = compare(a, b)
  for i in range(2, int(number) + 1):
    march = march.replace('1.%s.pdb' % (i-1), '1.%s.pdb' % i)
    b = AEV(march)
    diff = compare(a, b)
    for element, values in diff.items():
        for r_or_a, value in values.items():
          total[element][r_or_a].append(value[0])
    print(march)
  return total

def plotvalue(aev1, aev2, elename,pdb_filename):
  y1 = []
  y2 = []
  name = []
  d = aev1.get_AEVS()
  keys = d.keys()
  keys.sort(_sort)
  for key in keys:
    values = d[key]
    if elename in key:
      for r_or_a, value in values.items():
        for list1 in value:
          y1.append(list1)
        name.append(r_or_a)
  for element, values in aev2.get_AEVS().items():
    if elename in element:
      for value in values.values():
        for list2 in value:
          y2.append(list2)
  x = range(len(y1))
  plt.title("AEV difference with two %s atom" % elename)
  plt.xlabel("radial or angular of %s atom" % elename)
  plt.ylabel("value")
  plt.plot(x, y1, 'r')
  plt.plot(x, y2, 'g')
  plt.xticks(x[::8], name)
  plt.savefig('./difference/%s.jpg' % pdb_filename.replace('.pdb',''))
  plt.show()

def plotcompare(diff,ele,pdb_filename):
  name = []
  i = 0
  for element, values in diff.items():
    if ele in element:
      for r_or_a, value in values.items():
        c = randomcolor()
        plt.plot(range(len(value)), value, label = '%r'%r_or_a, color = c)
        plt.legend()
  plt.title("Correlation coefficient of %s atom"%ele)
  plt.xlabel("all files")
  plt.ylabel("value")
  plt.savefig('./corrcoee/%s.jpg' % pdb_filename.replace('.pdb',''))
  plt.show()

def main(pdb_file_name1, pdb_file_name2, number, elename):
  a = AEV(pdb_file_name1)
  b = AEV(pdb_file_name2)
  total = compare_all(pdb_file_name1, pdb_file_name2, number)
  print(total)
  plotvalue(a, b, elename, 'ethane.pdb')
  plotcompare(total, elename, pdb_file_name2)



if __name__ == '__main__':
    import os, sys
    main(*tuple(sys.argv[1:]))

