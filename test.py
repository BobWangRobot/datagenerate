import iotbx
from AEVclass import AEV
import numpy as np
import matplotlib.pyplot as plt

def compare(aev1, aev2):
  diff = aev2.get_AEVS()
  for element, values in aev1.get_AEVS().items():
    for r_or_a, value in values.items():
      if 0 in value:
        diff[element][r_or_a] = (1, 1)
      else:
        covalue = list(np.corrcoef(value, diff[element][r_or_a]))
        diff[element][r_or_a] = covalue[1]
  return diff

def plotvalue(aev1, aev2, elename):
  y1 = []
  y2 = []
  name = []
  for element, values in aev1.get_AEVS().items():
    if elename in element:
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
  plt.show()

def plotcompare(diff,ele):
  name = []
  y1 = []
  y2 = []
  for element, values in diff.items():
    if ele in element:
      for r_or_a, value in values.items():
        name.append(r_or_a)
        y1.append(value[0])
        y2.append(value[1])
  x = range(len(y1))
  plt.title("Correlation coefficient of %s atom"%ele)
  plt.xlabel("radial or angular to %s"%ele)
  plt.ylabel("value")
  plt.plot(x, y1, 'ob')
  plt.plot(x, y2, 'r')
  plt.xticks(x, name)
  plt.show()

def main(pdb_file_name1, pdb_file_name2, element):
  a = AEV(pdb_file_name1)
  b = AEV(pdb_file_name2)
  print(a.get_AEVS(), b.get_AEVS())
  plotcompare(compare(a, b), element)
  plotvalue(a, b, element)



if __name__ == '__main__':
    import os, sys
    main(*tuple(sys.argv[1:]))

