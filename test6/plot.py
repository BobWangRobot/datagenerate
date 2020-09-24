from __future__ import absolute_import, division, print_function
import sys
import time
import mmtbx
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
import __init__ as aev
import matplotlib.pyplot as plt
import xlwt

def plot(data1,data2):
  label = []
  result1 = []
  result2 = []
  for key, value in data1.items():
    label.append(key)
    result1.extend(value[7:])
  for value in data2.values():
    result2.extend(value[7:])
  x = range(len(result1))
  plt.title("angular figure")
  plt.xlabel("residue")
  plt.ylabel("value")
  plt.plot(x, result1, 'r', label='CA')
  plt.plot(x, result2, 'b', label='perfect_helix', linestyle=':')
  plt.xticks(x[::16], label)
  plt.legend()
  plt.savefig('angular.png')
  plt.show()

def main(filename1, filename2):
  t0 = time.time()
  pdb_inp1 = iotbx.pdb.input(file_name = filename1)
  model1 = mmtbx.model.manager(
    model_input   = pdb_inp1,
    build_grm     = True,
    log           = null_out())
  pdb_inp2 = iotbx.pdb.input(file_name=filename2)
  model2 = mmtbx.model.manager(
    model_input=pdb_inp2,
    build_grm=True,
    log=null_out())
  a = aev.AEV(model = model1)
  b = aev.AEV(model = model2)
  # c = aev.compare(a,1)
  d = aev.compare(b, 1)
  print(d)
  # data1 = a.MAEVs
  # data2 = b.MAEVs
  # plot(data1, data2)
  print('time', time.time()-t0)

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))