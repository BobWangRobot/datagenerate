import iotbx
from AEVclass import AEV
import numpy as np
import matplotlib.pyplot as plt

def compare(a,b):
  y1 = []
  y2 = []
  name = []
  diff = {}
  for atoms, values in a.get_AEVS().items():
    for atom, value in values.items():
      name.append(atom)
      y1.append(value)
  for atoms, values in b.get_AEVS().items():
    for atom, value in values.items():
      y2.append(value)
  for i in range(len(y1)):
    print(i,'\n',y1[i],y2[i])
    diff.setdefault(name[i], [])
    #print(name[i], np.corrcoef(y1[i], y2[i]))
    a = np.corrcoef(y1[i], y2[i])
    print(a[0])
    diff[name[i]].append(a[0])
  return diff

def plotvalue(a, b, atomname):
  y1 = []
  y2 = []
  name = []
  for atoms, values in a.get_AEVS().items():
    if atomname in atoms:
      for atom, value in values.items():
        for a in value:
          y1.append(a)
        name.append(atom)
  for atoms, values in b.get_AEVS().items():
    if atomname in atoms:
      for atom, value in values.items():
        for a in value:
          y2.append(a)
  x = range(len(y1))
  plt.plot(x,y1,'r')
  plt.plot(x,y2,'g')
  plt.xticks(x[::8], name)
  plt.show()

def main(pdb_file_name1, pdb_file_name2):
  a = AEV(pdb_file_name1)
  b = AEV(pdb_file_name2)
  print(compare(a,b))
  plotvalue(a, b, 'C')
  print(a.get_AEVS())

def plot(xvelue,yvelue):
  x = []


if __name__ == '__main__':
    import os, sys
    main(*tuple(sys.argv[1:]))

