from __future__ import absolute_import, division, print_function
import sys
import mmtbx
import iotbx.pdb
from libtbx.utils import null_out
from scitbx.array_family import flex
import mmtbx.secondary_structure

def get_ss(hierarchy, sec_str_from_pdb_file=None, method="ksdssp",
           use_recs=False):
  list1 = flex.double()
  if(use_recs): params = None
  else:
    params = mmtbx.secondary_structure.manager.get_default_ss_params()
    params.secondary_structure.protein.search_method=method
    params = params.secondary_structure
  ssm = mmtbx.secondary_structure.manager(
    pdb_hierarchy         = hierarchy,
    sec_str_from_pdb_file = sec_str_from_pdb_file,
    params                = params,
    log                   = null_out())
  if sec_str_from_pdb_file != None:
    print("- standard"*10)
    print(ssm.records_for_pdb_file())
    print("-" * 79)
  else:
    print("- %s"%method*10)
    print(ssm.records_for_pdb_file())
    print("-"*79)


def run(args):
  pdb_inp = iotbx.pdb.input(file_name = args[0])
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  v1 = get_ss(
    hierarchy=pdb_hierarchy,
    sec_str_from_pdb_file=pdb_inp.extract_secondary_structure())

if __name__ == '__main__':
  run(args=sys.argv[1:])
