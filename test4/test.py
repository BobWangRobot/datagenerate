import iotbx
from iotbx import pdb
import sys
from AEVclass3 import *

def name():
  pdb_inp = iotbx.pdb.input('per1.pdb')
  hierarchy = pdb_inp.construct_hierarchy()
  for atom_group in hierarchy.atom_groups():
    print(atom_group.resname)

def main(direction, scope, filename):
  a = AEV(direction,scope,pdb_file_name=filename)
  for a.five in a.generate_ca():
    for atoms in a.five:
      print(dir(atoms))

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))
