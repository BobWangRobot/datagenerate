import iotbx
from iotbx import pdb
import sys
import collections
from AEVclass4 import *

def name():
  pdb_inp = iotbx.pdb.input('per1.pdb')
  hierarchy = pdb_inp.construct_hierarchy()
  for atom_group in hierarchy.atom_groups():
    print(atom_group.resname)

def main(direction, scope, filename):
  a = AEV(direction,scope,pdb_file_name=filename)

  for a.five in a.generate_ca():

    print(a.five)
    # for atoms in a.five:
    #   print(atoms)
def test():
  b = collections.OrderedDict()
  b['a'] = 1
  b['b'] = 2
  print(b[1])

if __name__ == '__main__':
  # main(*tuple(sys.argv[1:]))
  test()