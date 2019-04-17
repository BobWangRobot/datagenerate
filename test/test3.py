import iotbx
import numpy as np
import matplotlib.pyplot as plt
import random
import sys
sys.path.append("..")
from AEVclass import AEV, radial_aev_class

#test 2 different element pdb file

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

def compare(pdb_file1, pdb_file2):
  aev1 = AEV(pdb_file1)
  aev2 = AEV(pdb_file2)
  diff = aev1.get_items()
  for element, values in aev2.get_AEVS().items():
    if element in aev1.get_AEVS().keys():
      print(element)
      for r_or_a, value in values.items():
        if r_or_a in aev1.get_AEVS()[element]:
          if 0 in value:
            diff[element][r_or_a] = [1]
          else:
            covalue = np.corrcoef(value, aev1.get_AEVS()[element][r_or_a]).tolist()
            diff[element][r_or_a] = [covalue[1][0]]
  #print(diff)
  return diff

def plotvalue(pdb1, pdb2, elename):
  y1 = []
  y2 = []
  name = []
  aev1 = AEV(pdb1)
  aev2 = AEV(pdb2)
  for ele, values in aev1.get_AEVS().items():
    if elename in ele:
      for r_or_a, value in values.items():
        for list1 in value:
          y1.append(list1)
        name.append(r_or_a)
  for element, values in aev2.get_AEVS().items():
    if elename in element:
      for value in values.values():
        for list2 in value:
          y2.append(list2)
  x = range(len(y2))
  plt.title("%s AEV difference with two %s atom" %(pdb2.replace('.pdb',''), elename))
  plt.xlabel("radial or angular of %s atom" % elename)
  plt.ylabel("value")
  plt.plot(x, y1, 'r', label='%s' % pdb1.replace('.pdb', ''))
  plt.plot(x, y2, 'g', label='%s' % pdb2.replace('.pdb', ''))
  plt.xticks(x[::8], name)
  plt.legend()
  plt.savefig('./difference/%s.jpg' % pdb2.replace('.pdb',''))
  plt.show()

def plotcompare(diff,ele,pdb_filename):
  name = []
  y1 = []
  y2 = []
  for element, values in diff.items():
    if ele in element:
      for r_or_a, value in values.items():
        name.append(r_or_a)
        y1.append(value[0])
  x = range(len(y1))
  plt.title("Correlation coefficient of %s atom"%ele)
  plt.xlabel("radial or angular to %s"%ele)
  plt.ylabel("value")
  plt.plot(x, y1, 'ob')
  plt.xticks(x, name)
  plt.savefig('./corrcoee/%s.jpg' % pdb_filename.replace('.pdb','%s'%ele))
  plt.show()

def main(pdb_file_name1, pdb_file_name2, elename):
  #print(AEV(pdb_file_name1).get_AEVS(), AEV(pdb_file_name2).get_AEVS())
  diff = compare(pdb_file_name1, pdb_file_name2)
  plotcompare(diff, elename, pdb_file_name1)
  # plotvalue(pdb_file_name1, pdb_file_name2, elename)
  # plotcompare(total, elename, pdb_file_name2)



if __name__ == '__main__':
    import os, sys
    main(*tuple(sys.argv[1:]))

