import sys, os
import numpy as np
import iotbx
import math
import collections
import time
import mmtbx
from iotbx import pdb
from AEVclass4 import *

perfect_helix_12 = """
ATOM      1  N   ALA A   1      -5.606  -2.251 -12.878  1.00  0.00           N
ATOM      2  CA  ALA A   1      -5.850  -1.194 -13.852  1.00  0.00           C
ATOM      3  C   ALA A   1      -5.186  -1.524 -15.184  1.00  0.00           C
ATOM      4  O   ALA A   1      -5.744  -1.260 -16.249  1.00  0.00           O
ATOM      5  CB  ALA A   1      -5.339   0.137 -13.324  1.00  0.00           C
ATOM      6  N   CYS A   2      -3.992  -2.102 -15.115  1.00  0.00           N
ATOM      7  CA  CYS A   2      -3.261  -2.499 -16.313  1.00  0.00           C
ATOM      8  C   CYS A   2      -3.961  -3.660 -17.011  1.00  0.00           C
ATOM      9  O   CYS A   2      -4.016  -3.716 -18.240  1.00  0.00           O
ATOM     10  CB  CYS A   2      -1.828  -2.894 -15.955  1.00  0.00           C
ATOM     11  SG  CYS A   2      -0.819  -1.533 -15.323  1.00  0.00           S
ATOM     12  N   GLU A   3      -4.492  -4.585 -16.219  1.00  0.00           N
ATOM     13  CA  GLU A   3      -5.216  -5.731 -16.755  1.00  0.00           C
ATOM     14  C   GLU A   3      -6.531  -5.289 -17.389  1.00  0.00           C
ATOM     15  O   GLU A   3      -6.939  -5.814 -18.425  1.00  0.00           O
ATOM     16  CB  GLU A   3      -5.488  -6.753 -15.651  1.00  0.00           C
ATOM     17  CG  GLU A   3      -4.238  -7.425 -15.107  1.00  0.00           C
ATOM     18  CD  GLU A   3      -4.542  -8.417 -14.002  1.00  0.00           C
ATOM     19  OE1 GLU A   3      -5.712  -8.490 -13.572  1.00  0.00           O
ATOM     20  OE2 GLU A   3      -3.610  -9.124 -13.561  1.00  0.00           O
ATOM     21  N   ASP A   4      -7.189  -4.323 -16.758  1.00  0.00           N
ATOM     22  CA  ASP A   4      -8.442  -3.785 -17.273  1.00  0.00           C
ATOM     23  C   ASP A   4      -8.205  -3.003 -18.561  1.00  0.00           C
ATOM     24  O   ASP A   4      -9.007  -3.065 -19.492  1.00  0.00           O
ATOM     25  CB  ASP A   4      -9.101  -2.881 -16.230  1.00  0.00           C
ATOM     26  CG  ASP A   4     -10.476  -2.406 -16.657  1.00  0.00           C
ATOM     27  OD1 ASP A   4     -11.425  -3.219 -16.622  1.00  0.00           O
ATOM     28  OD2 ASP A   4     -10.610  -1.220 -17.027  1.00  0.00           O
ATOM     29  N   GLY A   5      -7.099  -2.269 -18.604  1.00  0.00           N
ATOM     30  CA  GLY A   5      -6.735  -1.498 -19.787  1.00  0.00           C
ATOM     31  C   GLY A   5      -6.358  -2.423 -20.939  1.00  0.00           C
ATOM     32  O   GLY A   5      -6.687  -2.157 -22.094  1.00  0.00           O
ATOM     33  N   PHE A   6      -5.665  -3.509 -20.614  1.00  0.00           N
ATOM     34  CA  PHE A   6      -5.268  -4.493 -21.614  1.00  0.00           C
ATOM     35  C   PHE A   6      -6.485  -5.236 -22.153  1.00  0.00           C
ATOM     36  O   PHE A   6      -6.565  -5.533 -23.345  1.00  0.00           O
ATOM     37  CB  PHE A   6      -4.274  -5.490 -21.014  1.00  0.00           C
ATOM     38  CG  PHE A   6      -2.949  -4.883 -20.653  1.00  0.00           C
ATOM     39  CD1 PHE A   6      -1.984  -4.668 -21.623  1.00  0.00           C
ATOM     40  CD2 PHE A   6      -2.667  -4.527 -19.344  1.00  0.00           C
ATOM     41  CE1 PHE A   6      -0.762  -4.110 -21.294  1.00  0.00           C
ATOM     42  CE2 PHE A   6      -1.448  -3.968 -19.009  1.00  0.00           C
ATOM     43  CZ  PHE A   6      -0.495  -3.759 -19.986  1.00  0.00           C
ATOM     44  N   ILE A   7      -7.430  -5.532 -21.267  1.00  0.00           N
ATOM     45  CA  ILE A   7      -8.660  -6.212 -21.655  1.00  0.00           C
ATOM     46  C   ILE A   7      -9.529  -5.303 -22.518  1.00  0.00           C
ATOM     47  O   ILE A   7     -10.158  -5.756 -23.474  1.00  0.00           O
ATOM     48  CB  ILE A   7      -9.465  -6.668 -20.424  1.00  0.00           C
ATOM     49  CG1 ILE A   7      -9.887  -5.460 -19.580  1.00  0.00           C
ATOM     50  CG2 ILE A   7      -8.641  -7.639 -19.590  1.00  0.00           C
ATOM     51  CD1 ILE A   7     -10.877  -5.789 -18.477  1.00  0.00           C
ATOM     52  N   HIS A   8      -9.559  -4.021 -22.172  1.00  0.00           N
ATOM     53  CA  HIS A   8     -10.324  -3.039 -22.930  1.00  0.00           C
ATOM     54  C   HIS A   8      -9.706  -2.819 -24.306  1.00  0.00           C
ATOM     55  O   HIS A   8     -10.416  -2.660 -25.299  1.00  0.00           O
ATOM     56  CB  HIS A   8     -10.390  -1.714 -22.170  1.00  0.00           C
ATOM     57  CG  HIS A   8     -11.130  -1.798 -20.872  1.00  0.00           C
ATOM     58  ND1 HIS A   8     -12.504  -1.872 -20.802  1.00  0.00           N
ATOM     59  CD2 HIS A   8     -10.687  -1.820 -19.593  1.00  0.00           C
ATOM     60  CE1 HIS A   8     -12.876  -1.935 -19.537  1.00  0.00           C
ATOM     61  NE2 HIS A   8     -11.792  -1.905 -18.782  1.00  0.00           N
ATOM     62  N   LYS A   9      -8.378  -2.810 -24.356  1.00  0.00           N
ATOM     63  CA  LYS A   9      -7.658  -2.641 -25.613  1.00  0.00           C
ATOM     64  C   LYS A   9      -7.843  -3.861 -26.508  1.00  0.00           C
ATOM     65  O   LYS A   9      -7.980  -3.734 -27.725  1.00  0.00           O
ATOM     66  CB  LYS A   9      -6.170  -2.411 -25.347  1.00  0.00           C
ATOM     67  CG  LYS A   9      -5.859  -1.086 -24.669  1.00  0.00           C
ATOM     68  CD  LYS A   9      -4.365  -0.906 -24.459  1.00  0.00           C
ATOM     69  CE  LYS A   9      -4.055   0.409 -23.761  1.00  0.00           C
ATOM     70  NZ  LYS A   9      -2.595   0.594 -23.539  1.00  0.00           N
ATOM     71  N   MET A  10      -7.846  -5.040 -25.897  1.00  0.00           N
ATOM     72  CA  MET A  10      -8.046  -6.284 -26.631  1.00  0.00           C
ATOM     73  C   MET A  10      -9.473  -6.375 -27.160  1.00  0.00           C
ATOM     74  O   MET A  10      -9.704  -6.850 -28.272  1.00  0.00           O
ATOM     75  CB  MET A  10      -7.751  -7.485 -25.730  1.00  0.00           C
ATOM     76  CG  MET A  10      -6.281  -7.646 -25.376  1.00  0.00           C
ATOM     77  SD  MET A  10      -5.715  -6.430 -24.171  1.00  0.00           S
ATOM     78  CE  MET A  10      -4.130  -7.123 -23.708  1.00  0.00           C
ATOM     79  N   LEU A  11     -10.426  -5.917 -26.355  1.00  0.00           N
ATOM     80  CA  LEU A  11     -11.829  -5.915 -26.751  1.00  0.00           C
ATOM     81  C   LEU A  11     -12.072  -4.912 -27.873  1.00  0.00           C
ATOM     82  O   LEU A  11     -12.848  -5.170 -28.793  1.00  0.00           O
ATOM     83  CB  LEU A  11     -12.720  -5.579 -25.552  1.00  0.00           C
ATOM     84  CG  LEU A  11     -12.745  -6.598 -24.407  1.00  0.00           C
ATOM     85  CD1 LEU A  11     -13.557  -6.060 -23.238  1.00  0.00           C
ATOM     86  CD2 LEU A  11     -13.294  -7.941 -24.870  1.00  0.00           C
ATOM     87  N   ASN A  12     -11.403  -3.767 -27.788  1.00  0.00           N
ATOM     88  CA  ASN A  12     -11.519  -2.734 -28.811  1.00  0.00           C
ATOM     89  C   ASN A  12     -10.882  -3.194 -30.117  1.00  0.00           C
ATOM     90  O   ASN A  12     -11.397  -2.918 -31.201  1.00  0.00           O
ATOM     91  CB  ASN A  12     -10.860  -1.441 -28.338  1.00  0.00           C
ATOM     92  CG  ASN A  12     -11.640  -0.766 -27.227  1.00  0.00           C
ATOM     93  OD1 ASN A  12     -12.868  -0.703 -27.263  1.00  0.00           O
ATOM     94  ND2 ASN A  12     -10.926  -0.255 -26.230  1.00  0.00           N
TER
END"""

# It is format calss of corelation coefficient values.
class diff_class(OrderedDict):
  def __repr__(self):
    outl = '...\n'
    for key, item in self.items():
      outl += '  %s :' % (key)
      for key1,value in item.items():
        outl += ' %s: '%key1
        outl += '%0.2f, ' % value
      outl += '\n'
    return outl

# The AEV values of the perfect helix.
def generate_perfect_helix():
  perfect_helix = OrderedDict()
  a = AEV(raw_records=perfect_helix_12)
  a.generate_AEV()
  perfect_helix['B'] = a.BAEVs.values()[0]
  perfect_helix['M'] = a.MAEVs.values()[5]
  perfect_helix['E'] = a.EAEVs.values()[-1]
  return perfect_helix

# Comparing perfect helix with a target structure and getting correlation coefficient values.
# The result include 3 direction AEVs'Co-Co values. If the c-alpha doesn't have BAEVs or EAEVs
# the co-co values would be 0.
def compare(data):
  result = diff_class()
  perfect_helix = generate_perfect_helix()
  for key,value in perfect_helix.items():
    if key == 'B':
      for key1, value1 in data.BAEVs.items():
        if value1 != []:
          covalue = np.corrcoef(value, value1).tolist()
          result.setdefault(key1, OrderedDict())
          result[key1].setdefault(key, covalue[1][0])
        else:
          result.setdefault(key1, OrderedDict())
          result[key1].setdefault(key, 1)
    elif key == 'M':
      for key1, value1 in data.MAEVs.items():
        if value1 != []:
          covalue = np.corrcoef(value, value1).tolist()
          result.setdefault(key1, OrderedDict())
          result[key1].setdefault(key, covalue[1][0])
        else:
          result.setdefault(key1, OrderedDict())
          result[key1].setdefault(key, 1)
    elif key == 'E':
      for key1, value1 in data.EAEVs.items():
        if value1 != []:
          covalue = np.corrcoef(value, value1).tolist()
          result.setdefault(key1, OrderedDict())
          result[key1].setdefault(key, covalue[1][0])
        else:
          result.setdefault(key1, OrderedDict())
          result[key1].setdefault(key, 1)
  return result

# Geting heliecs of target structure.
# According to the threshold(percision), select helices of target structure.
# If the co-co values of atoms > threshold, the atom would be considered a part of a helix.
def HELIX_record(data, precision):
  start = []
  end = []
  M = 0
  i = 0
  precision = float(precision)
  for key,value in data.items():
    if start==[] and value['B'] > precision and value['M'] > precision - 0.1:
      start = key
      length = 1
    elif start and value['M'] > precision:
      M = 1
      length += 1
      if value['E'] > precision:
        end = key
    elif start and M == 1 and value['E'] + value['M'] > 2 * precision:
      end = key
      length += 1
    else:
      if start and end and M ==1 and length > 3:
        i += 1
        print("HELIX   {0:>2}  {0:>2} {1:>}  {2:>}   {3:>36}".format(i, start, end, length))
      start = []
      end = []
      M = 0
  return 0
# Writing results in a file.
def write_file(data,filename):
  name = str(filename + '.log')
  f = open(name, 'w')
  f.write(str(data))
  f.close()

def main(filename, precision=0.9, display_AEV=False, display_coco=False):
  t0 = time.time()
  a = AEV(pdb_file_name=filename)
  a.generate_AEV()
  if display_AEV:
    print('forward-only', a.BAEVs)
    print('all', a.MAEVs) # all AEV
    print('backward-only', a.EAEVs)
  b = compare(a)
  if display_coco:
    print(b)
  HELIX_record(b, precision)
  print('time', time.time()-t0)




if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))


