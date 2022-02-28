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
sys.path.append("..")
from AEVclass2 import *

#plot two parallel helices(same length)

#helix_5_12 from 1an0
helix1='''
ATOM      1  N   LEU A 165      50.506  10.993  25.761  1.00  2.00           N
ATOM      2  CA  LEU A 165      50.928  12.143  24.961  1.00  2.00           C
ATOM      3  C   LEU A 165      49.803  13.115  24.596  1.00  2.00           C
ATOM      4  O   LEU A 165      49.523  14.055  25.330  1.00 20.53           O
ATOM      5  CB  LEU A 165      51.627  11.639  23.699  1.00 20.53           C
ATOM      6  CG  LEU A 165      52.441  12.551  22.777  1.00 20.53           C
ATOM      7  CD1 LEU A 165      53.737  12.941  23.433  1.00 20.53           C
ATOM      8  CD2 LEU A 165      52.727  11.798  21.491  1.00 20.53           C
ATOM      9  N   LYS A 166      49.118  12.829  23.497  1.00 21.44           N
ATOM     10  CA  LYS A 166      48.008  13.654  23.010  1.00 21.44           C
ATOM     11  C   LYS A 166      47.093  14.097  24.142  1.00 21.44           C
ATOM     12  O   LYS A 166      46.610  15.230  24.168  1.00 35.87           O
ATOM     13  CB  LYS A 166      47.194  12.859  21.980  1.00 35.87           C
ATOM     14  CG  LYS A 166      46.036  13.609  21.354  1.00 35.87           C
ATOM     15  CD  LYS A 166      45.274  12.733  20.375  1.00 35.87           C
ATOM     16  CE  LYS A 166      44.151  13.508  19.717  1.00 35.87           C
ATOM     17  NZ  LYS A 166      43.261  14.153  20.718  1.00 35.87           N
ATOM     18  N   ASN A 167      46.884  13.200  25.094  1.00 23.47           N
ATOM     19  CA  ASN A 167      46.032  13.480  26.234  1.00 23.47           C
ATOM     20  C   ASN A 167      46.596  14.663  26.986  1.00 23.47           C
ATOM     21  O   ASN A 167      45.874  15.613  27.309  1.00 12.55           O
ATOM     22  CB  ASN A 167      45.971  12.274  27.157  1.00 12.55           C
ATOM     23  CG  ASN A 167      45.209  12.559  28.425  1.00 12.55           C
ATOM     24  OD1 ASN A 167      45.792  12.627  29.500  1.00 12.55           O
ATOM     25  ND2 ASN A 167      43.901  12.740  28.305  1.00 12.55           N
ATOM     26  N   VAL A 168      47.905  14.621  27.200  1.00  2.00           N
ATOM     27  CA  VAL A 168      48.596  15.686  27.900  1.00  2.00           C
ATOM     28  C   VAL A 168      48.161  16.998  27.272  1.00  2.00           C
ATOM     29  O   VAL A 168      47.716  17.914  27.966  1.00 20.74           O
ATOM     30  CB  VAL A 168      50.120  15.570  27.742  1.00 20.74           C
ATOM     31  CG1 VAL A 168      50.805  16.641  28.556  1.00 20.74           C
ATOM     32  CG2 VAL A 168      50.601  14.185  28.146  1.00 20.74           C
ATOM     33  N   PHE A 169      48.184  17.041  25.944  1.00 18.32           N
ATOM     34  CA  PHE A 169      47.820  18.245  25.227  1.00 18.32           C
ATOM     35  C   PHE A 169      46.332  18.581  25.184  1.00 18.32           C
ATOM     36  O   PHE A 169      45.954  19.750  25.019  1.00 14.82           O
ATOM     37  CB  PHE A 169      48.473  18.238  23.856  1.00 14.82           C
ATOM     38  CG  PHE A 169      49.956  18.437  23.917  1.00 14.82           C
ATOM     39  CD1 PHE A 169      50.481  19.689  24.215  1.00 14.82           C
ATOM     40  CD2 PHE A 169      50.828  17.377  23.721  1.00 14.82           C
ATOM     41  CE1 PHE A 169      51.855  19.882  24.320  1.00 14.82           C
ATOM     42  CE2 PHE A 169      52.212  17.562  23.826  1.00 14.82           C
ATOM     43  CZ  PHE A 169      52.722  18.815  24.125  1.00 14.82           C
ATOM     44  N   ASP A 170      45.477  17.581  25.366  1.00 21.62           N
ATOM     45  CA  ASP A 170      44.042  17.844  25.397  1.00 21.62           C
ATOM     46  C   ASP A 170      43.783  18.490  26.738  1.00 21.62           C
ATOM     47  O   ASP A 170      43.070  19.499  26.844  1.00 28.00           O
ATOM     48  CB  ASP A 170      43.255  16.552  25.294  1.00 28.00           C
ATOM     49  CG  ASP A 170      43.147  16.070  23.883  1.00 28.00           C
ATOM     50  OD1 ASP A 170      42.236  16.553  23.180  1.00 28.00           O
ATOM     51  OD2 ASP A 170      43.977  15.238  23.466  1.00 28.00           O
ATOM     52  N   GLU A 171      44.417  17.915  27.755  1.00 27.58           N
ATOM     53  CA  GLU A 171      44.325  18.416  29.113  1.00 27.58           C
ATOM     54  C   GLU A 171      44.859  19.846  29.140  1.00 27.58           C
ATOM     55  O   GLU A 171      44.329  20.704  29.844  1.00 32.00           O
ATOM     56  CB  GLU A 171      45.113  17.517  30.062  1.00 32.00           C
ATOM     57  CG  GLU A 171      44.464  16.162  30.311  1.00 32.00           C
ATOM     58  CD  GLU A 171      43.111  16.265  31.006  1.00 32.00           C
ATOM     59  OE1 GLU A 171      42.707  17.375  31.421  1.00 32.00           O
ATOM     60  OE2 GLU A 171      42.446  15.223  31.149  1.00 32.00           O
ATOM     61  N   ALA A 172      45.888  20.091  28.335  1.00 24.54           N
ATOM     62  CA  ALA A 172      46.497  21.416  28.201  1.00 24.54           C
ATOM     63  C   ALA A 172      45.524  22.373  27.510  1.00 24.54           C
ATOM     64  O   ALA A 172      45.417  23.544  27.882  1.00 19.75           O
ATOM     65  CB  ALA A 172      47.788  21.324  27.396  1.00 19.75           C
ATOM     66  N   ILE A 173      44.823  21.879  26.497  1.00 14.40           N
ATOM     67  CA  ILE A 173      43.875  22.696  25.767  1.00 14.40           C
ATOM     68  C   ILE A 173      42.713  23.044  26.672  1.00 14.40           C
ATOM     69  O   ILE A 173      42.362  24.215  26.812  1.00  9.13           O
ATOM     70  CB  ILE A 173      43.350  21.961  24.530  1.00  9.13           C
ATOM     71  CG1 ILE A 173      44.464  21.848  23.496  1.00  9.13           C
ATOM     72  CG2 ILE A 173      42.137  22.663  23.954  1.00  9.13           C
ATOM     73  CD1 ILE A 173      44.061  21.080  22.295  1.00  9.13           C
ATOM     74  N   LEU A 174      42.123  22.033  27.305  1.00  2.00           N
ATOM     75  CA  LEU A 174      40.984  22.258  28.196  1.00  2.00           C
ATOM     76  C   LEU A 174      41.343  23.143  29.386  1.00  2.00           C
ATOM     77  O   LEU A 174      40.552  23.987  29.782  1.00 24.72           O
ATOM     78  CB  LEU A 174      40.394  20.927  28.660  1.00 24.72           C
ATOM     79  CG  LEU A 174      39.846  20.061  27.514  1.00 24.72           C
ATOM     80  CD1 LEU A 174      39.875  18.579  27.865  1.00 24.72           C
ATOM     81  CD2 LEU A 174      38.449  20.512  27.137  1.00 24.72           C
ATOM     82  N   ALA A 175      42.553  23.008  29.914  1.00 39.89           N
ATOM     83  CA  ALA A 175      42.973  23.836  31.047  1.00 39.89           C
ATOM     84  C   ALA A 175      43.243  25.288  30.653  1.00 39.89           C
ATOM     85  O   ALA A 175      43.217  26.179  31.505  1.00 17.77           O
ATOM     86  CB  ALA A 175      44.216  23.246  31.708  1.00 17.77           C
ATOM     87  N   ALA A 176      43.490  25.515  29.363  1.00  2.00           N
ATOM     88  CA  ALA A 176      43.809  26.843  28.851  1.00  2.00           C
ATOM     89  C   ALA A 176      42.632  27.634  28.355  1.00  2.00           C
ATOM     90  O   ALA A 176      42.750  28.833  28.119  1.00 13.66           O
ATOM     91  CB  ALA A 176      44.836  26.734  27.752  1.00 13.66           C
'''
#helix_5_12 from 10gs
helix2='''
ATOM      1  N   GLY B  12      15.107  25.510   7.614  1.00 17.29           N
ATOM      2  CA  GLY B  12      14.919  24.170   8.147  1.00 19.03           C
ATOM      3  C   GLY B  12      16.043  23.580   8.979  1.00 21.68           C
ATOM      4  O   GLY B  12      16.525  24.201   9.930  1.00 19.53           O
ATOM      5  N   ARG B  13      16.494  22.393   8.587  1.00 18.60           N
ATOM      6  CA  ARG B  13      17.557  21.700   9.303  1.00 16.78           C
ATOM      7  C   ARG B  13      18.984  22.143   8.997  1.00 17.85           C
ATOM      8  O   ARG B  13      19.941  21.511   9.453  1.00 14.77           O
ATOM      9  CB  ARG B  13      17.419  20.189   9.123  1.00 17.91           C
ATOM     10  CG  ARG B  13      16.250  19.595   9.896  1.00 16.36           C
ATOM     11  CD  ARG B  13      16.233  18.086   9.786  1.00 17.76           C
ATOM     12  NE  ARG B  13      15.255  17.479  10.688  1.00 16.09           N
ATOM     13  CZ  ARG B  13      15.511  17.133  11.947  1.00 14.76           C
ATOM     14  NH1 ARG B  13      16.715  17.344  12.474  1.00  9.54           N
ATOM     15  NH2 ARG B  13      14.558  16.583  12.686  1.00 12.19           N
ATOM     16  N   CYS B  14      19.137  23.223   8.234  1.00 13.59           N
ATOM     17  CA  CYS B  14      20.473  23.722   7.915  1.00 12.78           C
ATOM     18  C   CYS B  14      20.777  25.056   8.583  1.00 14.39           C
ATOM     19  O   CYS B  14      21.923  25.518   8.566  1.00 17.05           O
ATOM     20  CB  CYS B  14      20.666  23.837   6.404  1.00 11.20           C
ATOM     21  SG  CYS B  14      20.802  22.248   5.584  1.00 18.79           S
ATOM     22  N   ALA B  15      19.754  25.658   9.187  1.00 12.33           N
ATOM     23  CA  ALA B  15      19.894  26.944   9.860  1.00 11.57           C
ATOM     24  C   ALA B  15      21.024  26.948  10.892  1.00 15.07           C
ATOM     25  O   ALA B  15      21.909  27.802  10.851  1.00 19.37           O
ATOM     26  CB  ALA B  15      18.577  27.330  10.509  1.00 12.32           C
ATOM     27  N   ALA B  16      21.022  25.956  11.777  1.00 18.57           N
ATOM     28  CA  ALA B  16      22.038  25.856  12.816  1.00 12.94           C
ATOM     29  C   ALA B  16      23.447  25.658  12.272  1.00 14.92           C
ATOM     30  O   ALA B  16      24.385  26.332  12.714  1.00 20.76           O
ATOM     31  CB  ALA B  16      21.683  24.747  13.791  1.00 14.79           C
ATOM     32  N   LEU B  17      23.600  24.755  11.305  1.00 14.18           N
ATOM     33  CA  LEU B  17      24.920  24.494  10.730  1.00 16.72           C
ATOM     34  C   LEU B  17      25.433  25.696   9.934  1.00 12.67           C
ATOM     35  O   LEU B  17      26.639  25.927   9.860  1.00 16.54           O
ATOM     36  CB  LEU B  17      24.928  23.195   9.900  1.00 18.29           C
ATOM     37  CG  LEU B  17      24.156  23.036   8.588  1.00 21.35           C
ATOM     38  CD1 LEU B  17      25.082  23.325   7.424  1.00 16.63           C
ATOM     39  CD2 LEU B  17      23.615  21.612   8.473  1.00 20.25           C
ATOM     40  N   ARG B  18      24.515  26.476   9.367  1.00 12.74           N
ATOM     41  CA  ARG B  18      24.891  27.674   8.618  1.00 18.88           C
ATOM     42  C   ARG B  18      25.375  28.751   9.588  1.00 16.51           C
ATOM     43  O   ARG B  18      26.380  29.416   9.331  1.00 13.17           O
ATOM     44  CB  ARG B  18      23.710  28.195   7.804  1.00 15.29           C
ATOM     45  CG  ARG B  18      23.418  27.381   6.566  1.00  8.72           C
ATOM     46  CD  ARG B  18      22.061  27.743   6.005  1.00 15.71           C
ATOM     47  NE  ARG B  18      21.810  27.073   4.734  1.00 15.01           N
ATOM     48  CZ  ARG B  18      20.604  26.797   4.254  1.00 12.98           C
ATOM     49  NH1 ARG B  18      19.517  27.128   4.940  1.00 13.28           N
ATOM     50  NH2 ARG B  18      20.486  26.195   3.078  1.00 20.09           N
ATOM     51  N   MET B  19      24.661  28.905  10.706  1.00 15.95           N
ATOM     52  CA  MET B  19      25.034  29.882  11.728  1.00 15.28           C
ATOM     53  C   MET B  19      26.402  29.532  12.294  1.00 17.96           C
ATOM     54  O   MET B  19      27.211  30.411  12.574  1.00 19.34           O
ATOM     55  CB  MET B  19      24.021  29.901  12.870  1.00 18.95           C
ATOM     56  CG  MET B  19      22.626  30.337  12.477  1.00 24.87           C
ATOM     57  SD  MET B  19      21.554  30.545  13.921  1.00 34.59           S
ATOM     58  CE  MET B  19      21.340  28.900  14.412  1.00 33.63           C
ATOM     59  N   LEU B  20      26.648  28.236  12.459  1.00 18.78           N
ATOM     60  CA  LEU B  20      27.916  27.744  12.985  1.00 17.18           C
ATOM     61  C   LEU B  20      29.064  28.113  12.046  1.00 20.02           C
ATOM     62  O   LEU B  20      30.083  28.655  12.485  1.00 16.87           O
ATOM     63  CB  LEU B  20      27.842  26.223  13.179  1.00 13.92           C
ATOM     64  CG  LEU B  20      28.971  25.482  13.907  1.00 19.86           C
ATOM     65  CD1 LEU B  20      28.451  24.135  14.389  1.00 18.43           C
ATOM     66  CD2 LEU B  20      30.198  25.303  13.017  1.00 13.99           C
ATOM     67  N   LEU B  21      28.897  27.801  10.760  1.00 19.50           N
ATOM     68  CA  LEU B  21      29.910  28.096   9.752  1.00 19.18           C
ATOM     69  C   LEU B  21      30.167  29.596   9.656  1.00 17.53           C
ATOM     70  O   LEU B  21      31.316  30.031   9.662  1.00 15.83           O
ATOM     71  CB  LEU B  21      29.484  27.549   8.389  1.00 15.49           C
ATOM     72  CG  LEU B  21      29.547  26.030   8.229  1.00 15.85           C
ATOM     73  CD1 LEU B  21      28.894  25.612   6.926  1.00 13.39           C
ATOM     74  CD2 LEU B  21      30.993  25.564   8.281  1.00 12.67           C
ATOM     75  N   ALA B  22      29.092  30.378   9.604  1.00 15.69           N
ATOM     76  CA  ALA B  22      29.191  31.831   9.515  1.00 16.05           C
ATOM     77  C   ALA B  22      29.941  32.414  10.711  1.00 19.44           C
ATOM     78  O   ALA B  22      30.927  33.135  10.548  1.00 23.53           O
ATOM     79  CB  ALA B  22      27.809  32.436   9.426  1.00 13.27           C
ATOM     80  N   ASP B  23      29.500  32.051  11.912  1.00 19.70           N
ATOM     81  CA  ASP B  23      30.104  32.536  13.145  1.00 16.42           C
ATOM     82  C   ASP B  23      31.563  32.119  13.318  1.00 20.57           C
ATOM     83  O   ASP B  23      32.342  32.832  13.950  1.00 27.27           O
ATOM     84  CB  ASP B  23      29.277  32.073  14.346  1.00 26.04           C
ATOM     85  CG  ASP B  23      29.571  32.867  15.601  1.00 27.38           C
ATOM     86  OD1 ASP B  23      29.380  34.102  15.590  1.00 28.68           O
ATOM     87  OD2 ASP B  23      29.966  32.245  16.603  1.00 29.14           O
'''

def main(filename1=None, filename2=None):
  y1 = []
  y2 = []
  name1 = []
  name2 = []
  if filename1 and filename2:
    a = AEV(pdb_file_name=filename1)
    b = AEV(pdb_file_name=filename2)
    for a.five in a.generate_ca():
      a.get_AEVs('all')
    for b.five in b.generate_ca():
      b.get_AEVs('all')
    for key1,value1 in a.AEVs.items():
      name1.append(key1)
      for v1 in value1.values():
        y1.extend(v1)
    for key2,value2 in b.AEVs.items():
      name2.append(key2)
      for v2 in value2.values():
        y2.extend(v2)
    print(a.AEVs, b.AEVs)
    x1 = range(len(y1))
    x2 = range(len(y2))
    plt.title("%s and %s AEV values " % (filename1.replace('.pdb', ''),
                                                     filename2.replace('.pdb', '')))
    plt.xlabel("atom number")
    plt.ylabel("value")
    plt.plot(x1, y1, 'r', label='%s' % filename1.replace('.pdb', ''))
    plt.plot(x2, y2, 'g', label='%s' % filename2.replace('.pdb', ''))
    plt.xticks(x1[::16], (a_name+b_name for a_name, b_name in zip(name1, name2)))
    plt.legend()
    plt.savefig('%s all.jpg' % (filename1.replace('.pdb','')+filename2.replace('.pdb','')))
    plt.show()
  else:
    a = AEV(raw_records=helix1)
    b = AEV(raw_records=helix2)
    for a.five in a.generate_ca():
      a.Rpart('all')
    for b.five in b.generate_ca():
      b.Rpart('all')
    for key1,value1 in a.AEVs.items():
      name1.append(key1)
      for v1 in value1.values():
        y1.extend(v1)
    for key2,value2 in b.AEVs.items():
      name2.append(key2)
      for v2 in value2.values():
        y2.extend(v2)
    x1 = range(len(y1))
    x2 = range(len(y2))
    plt.title("1an0-helix and 10gs-helix radial part AEV values" )
    plt.xlabel("atom number")
    plt.ylabel("value")
    plt.plot(x1, y1, 'r', label='1an0')
    plt.plot(x2, y2, 'g', label='10gs')
    plt.xticks(x1[::8], name1)
    plt.xticks(x2[::8], name2)
    plt.legend()
    plt.savefig('1an0-10gs.jpg' )
    plt.show()



if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))
