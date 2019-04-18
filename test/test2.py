import iotbx
import numpy as np
import matplotlib.pyplot as plt
import random
import sys
sys.path.append("..")
from AEVclass import AEV, radial_aev_class

#test 9 same element pdb file

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
    diff[element].setdefault(element, [])
    all1 = []
    all2 = []
    for r_or_a, value in values.items():
      covalue = np.corrcoef(value, diff[element][r_or_a]).tolist()
      diff[element][r_or_a] = [covalue[1][0]]
      for v1 in value:
        all1.append(v1)
      for v2 in aev2.get_AEVS()[element][r_or_a]:
        all2.append(v2)
    covalue = np.corrcoef(all1, all2).tolist()
    diff[element][element] = [covalue[1][0]]
  return diff

def compare_all(origin, march, number):
  a = AEV(origin)
  b = AEV(march)
  total = compare(a, b)
  for i in range(0, int(number) + 1):
    march = march.replace('1.%s.pdb' % (i-1), '1.%s.pdb' % i)
    b = AEV(march)
    diff = compare(a, b)
    for element, values in diff.items():
        for r_or_a, value in values.items():
          total[element][r_or_a].append(value[0])
    print(march)
  return total

def plotvalue(aev1, aev2, elename,number):
  a = AEV(aev1)
  y1 = []
  name = []
  d = a.get_AEVS()
  keys = d.keys()
  keys.sort(_sort)
  for key in keys:
    values = d[key]
    if elename in key:
      for r_or_a, value in values.items():
        for list1 in value:
          y1.append(list1)
        name.append(r_or_a)
  for i in range(0, int(number) + 1):
    y2 = []
    aev2 = aev2.replace('1.%s.pdb' % (i-1), '1.%s.pdb' % i)
    b = AEV(aev2)
    for element, values in b.get_AEVS().items():
      if elename in element:
        for value in values.values():
          for list2 in value:
            y2.append(list2)
    x = range(len(y1))
    plt.title("%s AEV difference with two %s atom" %(aev2.replace('.pdb',''), elename))
    plt.xlabel("radial or angular of %s atom" % elename)
    plt.ylabel("value")
    plt.plot(x, y1, 'r', label='%s'%aev1.replace('.pdb',''))
    plt.plot(x, y2, 'g', label='%s'%aev2.replace('.pdb',''))
    plt.xticks(x[::8], name)
    plt.legend()
    plt.savefig('./difference/%s.jpg' % aev2.replace('.pdb',''))
    plt.show()

def plotcompare(diff,ele,pdb_filename):
  display = radial_aev_class()
  for element, values in diff.items():
    if ele in element:
      display.setdefault(ele, {})
      for r_or_a, value in values.items():
        c = randomcolor()
        plt.plot(range(len(value)), value, label = '%r'%r_or_a, color = c)
        plt.legend()
        display[ele].update(values)
        print (display)
  plt.title("Correlation coefficient of %s atom"%ele)
  plt.xlabel("all files")
  plt.ylabel("value")
  plt.savefig('./corrcoee/%s.jpg' % pdb_filename.replace('.pdb',''))
  plt.show()

def main(pdb_file_name1,pdb_file_name2, elename, number):
  total = compare_all(pdb_file_name1, pdb_file_name2, number)
  #plotvalue(pdb_file_name1, pdb_file_name2, elename, number)
  plotcompare(total, elename, pdb_file_name2)



if __name__ == '__main__':
    import os, sys
    main(*tuple(sys.argv[1:]))

