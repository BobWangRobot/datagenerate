import os, sys
import time
import mmtbx

from iotbx import pdb

perfect_helix = '''ATOM      1  N   ALA A   1      -5.606  -2.251 -12.878  1.00  0.00           N
ATOM      2  CA  ALA A   1      -5.850  -1.194 -13.852  1.00  0.00           C
ATOM      3  C   ALA A   1      -5.186  -1.524 -15.184  1.00  0.00           C
ATOM      4  O   ALA A   1      -5.744  -1.260 -16.249  1.00  0.00           O
ATOM      5  CB  ALA A   1      -5.339   0.137 -13.324  1.00  0.00           C
ATOM      1  N   CYS A   2      -3.992  -2.102 -15.115  1.00  0.00           N
ATOM      2  CA  CYS A   2      -3.261  -2.499 -16.313  1.00  0.00           C
ATOM      3  C   CYS A   2      -3.961  -3.660 -17.011  1.00  0.00           C
ATOM      4  O   CYS A   2      -4.016  -3.716 -18.240  1.00  0.00           O
ATOM      5  CB  CYS A   2      -1.828  -2.894 -15.955  1.00  0.00           C
ATOM      6  SG  CYS A   2      -0.819  -1.533 -15.323  1.00  0.00           S
ATOM      1  N   GLU A   3      -4.492  -4.585 -16.219  1.00  0.00           N
ATOM      2  CA  GLU A   3      -5.216  -5.731 -16.755  1.00  0.00           C
ATOM      3  C   GLU A   3      -6.531  -5.289 -17.389  1.00  0.00           C
ATOM      4  O   GLU A   3      -6.939  -5.814 -18.425  1.00  0.00           O
ATOM      5  CB  GLU A   3      -5.488  -6.753 -15.651  1.00  0.00           C
ATOM      6  CG  GLU A   3      -4.238  -7.425 -15.107  1.00  0.00           C
ATOM      7  CD  GLU A   3      -4.542  -8.417 -14.002  1.00  0.00           C
ATOM      8  OE1 GLU A   3      -5.712  -8.490 -13.572  1.00  0.00           O
ATOM      9  OE2 GLU A   3      -3.610  -9.124 -13.561  1.00  0.00           O
ATOM      1  N   ASP A   4      -7.189  -4.323 -16.758  1.00  0.00           N
ATOM      2  CA  ASP A   4      -8.442  -3.785 -17.273  1.00  0.00           C
ATOM      3  C   ASP A   4      -8.205  -3.003 -18.561  1.00  0.00           C
ATOM      4  O   ASP A   4      -9.007  -3.065 -19.492  1.00  0.00           O
ATOM      5  CB  ASP A   4      -9.101  -2.881 -16.230  1.00  0.00           C
ATOM      6  CG  ASP A   4     -10.476  -2.406 -16.657  1.00  0.00           C
ATOM      7  OD1 ASP A   4     -11.425  -3.219 -16.622  1.00  0.00           O
ATOM      8  OD2 ASP A   4     -10.610  -1.220 -17.027  1.00  0.00           O
ATOM      1  N   GLY A   5      -7.099  -2.269 -18.604  1.00  0.00           N
ATOM      2  CA  GLY A   5      -6.735  -1.498 -19.787  1.00  0.00           C
ATOM      3  C   GLY A   5      -6.358  -2.423 -20.939  1.00  0.00           C
ATOM      4  O   GLY A   5      -6.687  -2.157 -22.094  1.00  0.00           O
ATOM      1  N   PHE A   6      -5.665  -3.509 -20.614  1.00  0.00           N
ATOM      2  CA  PHE A   6      -5.268  -4.493 -21.614  1.00  0.00           C
ATOM      3  C   PHE A   6      -6.485  -5.236 -22.153  1.00  0.00           C
ATOM      4  O   PHE A   6      -6.565  -5.533 -23.345  1.00  0.00           O
ATOM      5  CB  PHE A   6      -4.274  -5.490 -21.014  1.00  0.00           C
ATOM      6  CG  PHE A   6      -2.949  -4.883 -20.653  1.00  0.00           C
ATOM      7  CD1 PHE A   6      -1.984  -4.668 -21.623  1.00  0.00           C
ATOM      8  CD2 PHE A   6      -2.667  -4.527 -19.344  1.00  0.00           C
ATOM      9  CE1 PHE A   6      -0.762  -4.110 -21.294  1.00  0.00           C
ATOM     10  CE2 PHE A   6      -1.448  -3.968 -19.009  1.00  0.00           C
ATOM     11  CZ  PHE A   6      -0.495  -3.759 -19.986  1.00  0.00           C
ATOM      1  N   ILE A   7      -7.430  -5.532 -21.267  1.00  0.00           N
ATOM      2  CA  ILE A   7      -8.660  -6.212 -21.655  1.00  0.00           C
ATOM      3  C   ILE A   7      -9.529  -5.303 -22.518  1.00  0.00           C
ATOM      4  O   ILE A   7     -10.158  -5.756 -23.474  1.00  0.00           O
ATOM      5  CB  ILE A   7      -9.465  -6.668 -20.424  1.00  0.00           C
ATOM      6  CG1 ILE A   7      -9.887  -5.460 -19.580  1.00  0.00           C
ATOM      7  CG2 ILE A   7      -8.641  -7.639 -19.590  1.00  0.00           C
ATOM      8  CD1 ILE A   7     -10.877  -5.789 -18.477  1.00  0.00           C
ATOM      1  N   HIS A   8      -9.559  -4.021 -22.172  1.00  0.00           N
ATOM      2  CA  HIS A   8     -10.324  -3.039 -22.930  1.00  0.00           C
ATOM      3  C   HIS A   8      -9.706  -2.819 -24.306  1.00  0.00           C
ATOM      4  O   HIS A   8     -10.416  -2.660 -25.299  1.00  0.00           O
ATOM      5  CB  HIS A   8     -10.390  -1.714 -22.170  1.00  0.00           C
ATOM      6  CG  HIS A   8     -11.130  -1.798 -20.872  1.00  0.00           C
ATOM      7  ND1 HIS A   8     -12.504  -1.872 -20.802  1.00  0.00           N
ATOM      8  CD2 HIS A   8     -10.687  -1.820 -19.593  1.00  0.00           C
ATOM      9  CE1 HIS A   8     -12.876  -1.935 -19.537  1.00  0.00           C
ATOM     10  NE2 HIS A   8     -11.792  -1.905 -18.782  1.00  0.00           N
ATOM      1  N   LYS A   9      -8.378  -2.810 -24.356  1.00  0.00           N
ATOM      2  CA  LYS A   9      -7.658  -2.641 -25.613  1.00  0.00           C
ATOM      3  C   LYS A   9      -7.843  -3.861 -26.508  1.00  0.00           C
ATOM      4  O   LYS A   9      -7.980  -3.734 -27.725  1.00  0.00           O
ATOM      5  CB  LYS A   9      -6.170  -2.411 -25.347  1.00  0.00           C
ATOM      6  CG  LYS A   9      -5.859  -1.086 -24.669  1.00  0.00           C
ATOM      7  CD  LYS A   9      -4.365  -0.906 -24.459  1.00  0.00           C
ATOM      8  CE  LYS A   9      -4.055   0.409 -23.761  1.00  0.00           C
ATOM      9  NZ  LYS A   9      -2.595   0.594 -23.539  1.00  0.00           N
ATOM      1  N   MET A  10      -7.846  -5.040 -25.897  1.00  0.00           N
ATOM      2  CA  MET A  10      -8.046  -6.284 -26.631  1.00  0.00           C
ATOM      3  C   MET A  10      -9.473  -6.375 -27.160  1.00  0.00           C
ATOM      4  O   MET A  10      -9.704  -6.850 -28.272  1.00  0.00           O
ATOM      5  CB  MET A  10      -7.751  -7.485 -25.730  1.00  0.00           C
ATOM      6  CG  MET A  10      -6.281  -7.646 -25.376  1.00  0.00           C
ATOM      7  SD  MET A  10      -5.715  -6.430 -24.171  1.00  0.00           S
ATOM      8  CE  MET A  10      -4.130  -7.123 -23.708  1.00  0.00           C
ATOM      1  N   LEU A  11     -10.426  -5.917 -26.355  1.00  0.00           N
ATOM      2  CA  LEU A  11     -11.829  -5.915 -26.751  1.00  0.00           C
ATOM      3  C   LEU A  11     -12.072  -4.912 -27.873  1.00  0.00           C
ATOM      4  O   LEU A  11     -12.848  -5.170 -28.793  1.00  0.00           O
ATOM      5  CB  LEU A  11     -12.720  -5.579 -25.552  1.00  0.00           C
ATOM      6  CG  LEU A  11     -12.745  -6.598 -24.407  1.00  0.00           C
ATOM      7  CD1 LEU A  11     -13.557  -6.060 -23.238  1.00  0.00           C
ATOM      8  CD2 LEU A  11     -13.294  -7.941 -24.870  1.00  0.00           C
ATOM      1  N   ASN A  12     -11.403  -3.767 -27.788  1.00  0.00           N
ATOM      2  CA  ASN A  12     -11.519  -2.734 -28.811  1.00  0.00           C
ATOM      3  C   ASN A  12     -10.882  -3.194 -30.117  1.00  0.00           C
ATOM      4  O   ASN A  12     -11.397  -2.918 -31.201  1.00  0.00           O
ATOM      5  CB  ASN A  12     -10.860  -1.441 -28.338  1.00  0.00           C
ATOM      6  CG  ASN A  12     -11.640  -0.766 -27.227  1.00  0.00           C
ATOM      7  OD1 ASN A  12     -12.868  -0.703 -27.263  1.00  0.00           O
ATOM      8  ND2 ASN A  12     -10.926  -0.255 -26.230  1.00  0.00           N
ATOM      1  N   GLN A  13      -9.760  -3.897 -30.006  1.00  0.00           N
ATOM      2  CA  GLN A  13      -9.067  -4.426 -31.174  1.00  0.00           C
ATOM      3  C   GLN A  13      -9.884  -5.532 -31.832  1.00  0.00           C
ATOM      4  O   GLN A  13      -9.933  -5.636 -33.057  1.00  0.00           O
ATOM      5  CB  GLN A  13      -7.688  -4.959 -30.777  1.00  0.00           C
ATOM      6  CG  GLN A  13      -6.848  -5.450 -31.945  1.00  0.00           C
ATOM      7  CD  GLN A  13      -5.470  -5.917 -31.518  1.00  0.00           C
ATOM      8  OE1 GLN A  13      -5.141  -5.911 -30.332  1.00  0.00           O
ATOM      9  NE2 GLN A  13      -4.656  -6.324 -32.485  1.00  0.00           N
ATOM      1  N   PRO A  14     -10.522  -6.357 -31.008  1.00  0.00           N
ATOM      2  CA  PRO A  14     -11.364  -7.438 -31.505  1.00  0.00           C
ATOM      3  C   PRO A  14     -12.615  -6.883 -32.177  1.00  0.00           C
ATOM      4  O   PRO A  14     -13.068  -7.405 -33.195  1.00  0.00           O
ATOM      5  CB  PRO A  14     -11.735  -8.208 -30.242  1.00  0.00           C
ATOM      6  CG  PRO A  14     -10.627  -7.921 -29.292  1.00  0.00           C
ATOM      7  CD  PRO A  14     -10.227  -6.498 -29.550  1.00  0.00           C
ATOM      1  N   SER A  15     -13.168  -5.822 -31.598  1.00  0.00           N
ATOM      2  CA  SER A  15     -14.348  -5.173 -32.155  1.00  0.00           C
ATOM      3  C   SER A  15     -14.013  -4.480 -33.471  1.00  0.00           C
ATOM      4  O   SER A  15     -14.808  -4.494 -34.411  1.00  0.00           O
ATOM      5  CB  SER A  15     -14.915  -4.157 -31.162  1.00  0.00           C
ATOM      6  OG  SER A  15     -16.045  -3.493 -31.701  1.00  0.00           O
ATOM      1  N   ARG A  16     -12.832  -3.875 -33.530  1.00  0.00           N
ATOM      2  CA  ARG A  16     -12.373  -3.203 -34.740  1.00  0.00           C
ATOM      3  C   ARG A  16     -12.089  -4.214 -35.844  1.00  0.00           C
ATOM      4  O   ARG A  16     -12.376  -3.965 -37.015  1.00  0.00           O
ATOM      5  CB  ARG A  16     -11.112  -2.383 -34.449  1.00 10.00           C
ATOM      6  CG  ARG A  16     -11.306  -1.268 -33.426  1.00 10.00           C
ATOM      7  CD  ARG A  16     -12.240  -0.167 -33.923  1.00 10.00           C
ATOM      8  NE  ARG A  16     -11.698   0.544 -35.083  1.00 10.00           N
ATOM      9  CZ  ARG A  16     -12.389   1.321 -35.921  1.00 10.00           C
ATOM     10  NH1 ARG A  16     -13.697   1.536 -35.778  1.00 10.00           N
ATOM     11  NH2 ARG A  16     -11.755   1.900 -36.932  1.00 10.00           N
ATOM      1  N   THR A  17     -11.525  -5.355 -35.463  1.00  0.00           N
ATOM      2  CA  THR A  17     -11.229  -6.420 -36.413  1.00  0.00           C
ATOM      3  C   THR A  17     -12.516  -7.047 -36.938  1.00  0.00           C
ATOM      4  O   THR A  17     -12.617  -7.387 -38.117  1.00  0.00           O
ATOM      5  CB  THR A  17     -10.365  -7.521 -35.771  1.00  0.00           C
ATOM      6  OG1 THR A  17     -11.073  -8.112 -34.675  1.00  0.00           O
ATOM      7  CG2 THR A  17      -9.048  -6.946 -35.272  1.00  0.00           C
ATOM      1  N   TRP A  18     -13.497  -7.197 -36.054  1.00  0.00           N
ATOM      2  CA  TRP A  18     -14.791  -7.754 -36.431  1.00  0.00           C
ATOM      3  C   TRP A  18     -15.545  -6.794 -37.344  1.00  0.00           C
ATOM      4  O   TRP A  18     -16.211  -7.216 -38.290  1.00  0.00           O
ATOM      5  CB  TRP A  18     -15.626  -8.050 -35.184  1.00  0.00           C
ATOM      6  CG  TRP A  18     -15.072  -9.158 -34.343  1.00  0.00           C
ATOM      7  CD1 TRP A  18     -14.290  -9.033 -33.231  1.00  0.00           C
ATOM      8  CD2 TRP A  18     -15.258 -10.565 -34.544  1.00  0.00           C
ATOM      9  NE1 TRP A  18     -13.977 -10.273 -32.729  1.00  0.00           N
ATOM     10  CE2 TRP A  18     -14.560 -11.230 -33.517  1.00  0.00           C
ATOM     11  CE3 TRP A  18     -15.947 -11.326 -35.494  1.00  0.00           C
ATOM     12  CZ2 TRP A  18     -14.530 -12.619 -33.412  1.00  0.00           C
ATOM     13  CZ3 TRP A  18     -15.917 -12.705 -35.388  1.00  0.00           C
ATOM     14  CH2 TRP A  18     -15.214 -13.337 -34.355  1.00  0.00           C
ATOM      1  N   VAL A  19     -15.436  -5.502 -37.054  1.00  0.00           N
ATOM      2  CA  VAL A  19     -16.080  -4.477 -37.866  1.00  0.00           C
ATOM      3  C   VAL A  19     -15.427  -4.388 -39.241  1.00  0.00           C
ATOM      4  O   VAL A  19     -16.105  -4.196 -40.250  1.00  0.00           O
ATOM      5  CB  VAL A  19     -16.005  -3.094 -37.188  1.00  0.00           C
ATOM      6  CG1 VAL A  19     -16.591  -2.020 -38.095  1.00  0.00           C
ATOM      7  CG2 VAL A  19     -16.724  -3.118 -35.845  1.00  0.00           C
ATOM      1  N   TYR A  20     -14.106  -4.527 -39.271  1.00  0.00           N
ATOM      2  CA  TYR A  20     -13.360  -4.496 -40.524  1.00  0.00           C
ATOM      3  C   TYR A  20     -13.670  -5.726 -41.368  1.00  0.00           C
ATOM      4  O   TYR A  20     -13.779  -5.640 -42.591  1.00  0.00           O
ATOM      5  CB  TYR A  20     -11.857  -4.420 -40.246  1.00  0.00           C
ATOM      6  CG  TYR A  20     -11.413  -3.111 -39.633  1.00  0.00           C
ATOM      7  CD1 TYR A  20     -12.247  -1.999 -39.641  1.00  0.00           C
ATOM      8  CD2 TYR A  20     -10.161  -2.986 -39.045  1.00  0.00           C
ATOM      9  CE1 TYR A  20     -11.846  -0.801 -39.082  1.00  0.00           C
ATOM     10  CE2 TYR A  20      -9.751  -1.791 -38.482  1.00  0.00           C
ATOM     11  CZ  TYR A  20     -10.598  -0.703 -38.504  1.00  0.00           C
ATOM     12  OH  TYR A  20     -10.194   0.489 -37.945  1.00  0.00           O
'''
def get_geometry_restraints_manager(pdb_filename=None, raw_records=None):
  t0=time.time()
  from mmtbx.monomer_library import server
  from mmtbx.monomer_library import pdb_interpretation
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  processed_pdb = pdb_interpretation.process(mon_lib_srv, ener_lib, raw_records=raw_records,file_name=pdb_filename)
  geometry_restraints_manager = processed_pdb.geometry_restraints_manager()
  print 'time',time.time()-t0
  return geometry_restraints_manager

def generate_ca(filename=None, raw_records=None):
  from mmtbx.conformation_dependent_library import generate_protein_fragments
  if filename:
    pdb_inp = pdb.input(filename)
  else:
    pdb_inp = pdb.input(lines=raw_records, source_info='perfect_helix')
  hierarchy = pdb_inp.construct_hierarchy()
  hierarchy.reset_atom_i_seqs()
  geometry_restraints_manager = get_geometry_restraints_manager(pdb_filename=filename, raw_records=raw_records)
  hierarchy.reset_i_seq_if_necessary()
  for five in generate_protein_fragments(hierarchy,geometry_restraints_manager,length=5):
    rc = []
    for atom in five.atoms():
      if atom.name==' CA ':
        rc.append(atom)
    if len(rc)==5:
      yield rc


def main(filename=None):
  if filename:
    for five in generate_ca(filename):
      print five
  else:
    for five in generate_ca(raw_records=perfect_helix):
      print five

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))
