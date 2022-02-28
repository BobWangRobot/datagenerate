from __future__ import absolute_import, division, print_function
import mmtbx
import iotbx.pdb
import mmtbx.model
from libtbx import easy_run
from libtbx.utils import null_out
import __init__ as aev
import os
import sys
import matplotlib.pyplot as plt


pdb_str = """
ATOM      1  N   GLY A   1      -5.606  -2.251 -12.878  1.00  0.00           N
ATOM      2  CA  GLY A   1      -5.850  -1.194 -13.852  1.00  0.00           C
ATOM      3  C   GLY A   1      -5.186  -1.524 -15.184  1.00  0.00           C
ATOM      4  O   GLY A   1      -5.744  -1.260 -16.249  1.00  0.00           O
ATOM      1  N   GLY A   2      -3.992  -2.102 -15.115  1.00  0.00           N
ATOM      2  CA  GLY A   2      -3.261  -2.499 -16.313  1.00  0.00           C
ATOM      3  C   GLY A   2      -3.961  -3.660 -17.011  1.00  0.00           C
ATOM      4  O   GLY A   2      -4.016  -3.716 -18.240  1.00  0.00           O
ATOM      1  N   GLY A   3      -4.492  -4.585 -16.219  1.00  0.00           N
ATOM      2  CA  GLY A   3      -5.216  -5.731 -16.755  1.00  0.00           C
ATOM      3  C   GLY A   3      -6.531  -5.289 -17.389  1.00  0.00           C
ATOM      4  O   GLY A   3      -6.939  -5.814 -18.425  1.00  0.00           O
ATOM      1  N   GLY A   4      -7.189  -4.323 -16.758  1.00  0.00           N
ATOM      2  CA  GLY A   4      -8.442  -3.785 -17.273  1.00  0.00           C
ATOM      3  C   GLY A   4      -8.205  -3.003 -18.561  1.00  0.00           C
ATOM      4  O   GLY A   4      -9.007  -3.065 -19.492  1.00  0.00           O
ATOM      1  N   GLY A   5      -7.099  -2.269 -18.604  1.00  0.00           N
ATOM      2  CA  GLY A   5      -6.735  -1.498 -19.787  1.00  0.00           C
ATOM      3  C   GLY A   5      -6.358  -2.423 -20.939  1.00  0.00           C
ATOM      4  O   GLY A   5      -6.687  -2.157 -22.094  1.00  0.00           O
ATOM      1  N   GLY A   6      -5.665  -3.509 -20.614  1.00  0.00           N
ATOM      2  CA  GLY A   6      -5.268  -4.493 -21.614  1.00  0.00           C
ATOM      3  C   GLY A   6      -6.485  -5.236 -22.153  1.00  0.00           C
ATOM      4  O   GLY A   6      -6.565  -5.533 -23.345  1.00  0.00           O
ATOM      1  N   GLY A   7      -7.430  -5.532 -21.267  1.00  0.00           N
ATOM      2  CA  GLY A   7      -8.660  -6.212 -21.655  1.00  0.00           C
ATOM      3  C   GLY A   7      -9.529  -5.303 -22.518  1.00  0.00           C
ATOM      4  O   GLY A   7     -10.158  -5.756 -23.474  1.00  0.00           O
ATOM      1  N   GLY A   8      -9.559  -4.021 -22.172  1.00  0.00           N
ATOM      2  CA  GLY A   8     -10.324  -3.039 -22.930  1.00  0.00           C
ATOM      3  C   GLY A   8      -9.706  -2.819 -24.306  1.00  0.00           C
ATOM      4  O   GLY A   8     -10.416  -2.660 -25.299  1.00  0.00           O
ATOM      1  N   GLY A   9      -8.378  -2.810 -24.356  1.00  0.00           N
ATOM      2  CA  GLY A   9      -7.658  -2.641 -25.613  1.00  0.00           C
ATOM      3  C   GLY A   9      -7.843  -3.861 -26.508  1.00  0.00           C
ATOM      4  O   GLY A   9      -7.980  -3.734 -27.725  1.00  0.00           O
ATOM      1  N   GLY A  10      -7.846  -5.040 -25.897  1.00  0.00           N
ATOM      2  CA  GLY A  10      -8.046  -6.284 -26.631  1.00  0.00           C
ATOM      3  C   GLY A  10      -9.473  -6.375 -27.160  1.00  0.00           C
ATOM      4  O   GLY A  10      -9.704  -6.850 -28.272  1.00  0.00           O
"""
def compute(data):
  Min = 0
  Max = 0
  i = 0
  B_mean = 0
  E_mean = 0
  M_mean = 0
  for number in data.values():
    i += 1
    for key, value in number.items():
      if value == None:
        value = 0
      if Min > value:
        Min = value
      if Max < value:
        Max = value
      if key=='B':
        B_mean += value
      if key=='M':
        M_mean += value
      if key=='E':
        E_mean += value
  B_mean = B_mean / i
  E_mean = E_mean / i
  M_mean = M_mean / i
  # print("mean of forward:", B_mean)
  # print("mean of backward:", E_mean)
  # print("mean of middle:", M_mean)
  # print("Min value:%r  Max value:%r"%(Min, Max))
  return B_mean, E_mean, M_mean, Min, Max

def plot(B,E,M,Min,Max,RMSD):
  x = range(len(B))
  plt.title("RMSD=%r figure"%float(RMSD))
  plt.xlabel("Perture structrue number")
  plt.ylabel("CC value")
  plt.plot(x, B, 'r', label='backward_means')
  plt.plot(x, E, 'g', label='forward_means')
  plt.plot(x, M, 'b', label='middle_means')
  plt.plot(x, Min, 'y', label='Min value')
  plt.plot(x, Max, 'm', label='max value')
  plt.legend(bbox_to_anchor=(0.95, 0.1), loc=4)
  plt.savefig("New_RMSD=%r.jpg"%float(RMSD))
  plt.clf()
  # plt.show()


def run(num):
  path = '/home/bob/datagenerate/test6/per_structure/per%s/'%num
  of = open("".join([path, "INDEX"]), "r")
  files = ["".join(f).strip() for f in of.readlines()]
  RMSD = files[0][46:49]
  B_mean = []
  E_mean = []
  M_mean = []
  Min = []
  Max = []
  j = 0
  for f in files:
    j += 1
    sys.stdout.write("process:%s"%(j) +'%' + '\r')
    sys.stdout.flush()
    pdb_inp = iotbx.pdb.input(file_name=f)
    model = mmtbx.model.manager(
      model_input=pdb_inp,
      build_grm=True,
      log=null_out())
    a = aev.AEV(model = model)
    b = aev.compare(a)
    list = compute(b)
    B_mean.append(list[0])
    E_mean.append(list[1])
    M_mean.append(list[2])
    Min.append(list[3])
    Max.append(list[4])
    if not os.path.exists(f):
      files.remove(f)
  print('\n')
  of.close()
  plot(B_mean, E_mean, M_mean, Min, Max, RMSD)

if __name__ == '__main__':
  for i in range(0, 51, 5):
    if i != 0:
      run(i/10)
