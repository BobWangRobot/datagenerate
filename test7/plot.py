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

def main(precesion):
  data = []
  path = '/home/bob/Desktop/datagenerate/test6/per_structure/per1.0/'
  of = open("".join([path, "INDEX"]), "r")
  files = ["".join(f).strip() for f in of.readlines()]
  h_list = []
  i = 0
  non_h = []
  for f in files:
    pdb_inp = iotbx.pdb.input(file_name=f)
    model = mmtbx.model.manager(
      model_input=pdb_inp,
      build_grm=True,
      log=null_out())
    a = aev.AEV(model = model)
    b = aev.compare(a)
    result = aev.format_HELIX_records_from_AEV(b, precesion)
    if result!=[]:
      h_list.append(result)
    else:
      non_h.append(i)
    i += 1
    # sys.stdout.write("process:%s" %i + '%' + '\r')
    # sys.stdout.flush()
  # print(len(h_list), non_h)
  len_number = len(h_list)
  data.append(len_number)
  data.append(non_h)
  data = str(data)
  return data

if __name__ == '__main__':
  t0 = time.time()
  # main(*tuple(sys.argv[1:]))
  file = open("log1.log", 'a')
  for i in range(205, 295, 2):
    j = (1-(280 - i)/70) * 100
    sys.stdout.write("process:%s" % j + '%' + '\r')
    sys.stdout.flush()
    i = i/100
    a = main(i)
    file.write(str(i) + a + '\n')
  file.close()
  print('time', time.time() - t0)
