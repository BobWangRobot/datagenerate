from __future__ import absolute_import, division, print_function
import sys
import mmtbx
import iotbx.pdb
from libtbx.utils import null_out
from scitbx.array_family import flex
import mmtbx.secondary_structure
import xlwt

# def get_ss(hierarchy, sec_str_from_pdb_file=None, method="ksdssp",
#            use_recs=False):
#   if(use_recs): params = None
#   else:
#     params = mmtbx.secondary_structure.manager.get_default_ss_params()
#     params.secondary_structure.protein.search_method=method
#     params = params.secondary_structure
#   ssm = mmtbx.secondary_structure.manager(
#     pdb_hierarchy         = hierarchy,
#     sec_str_from_pdb_file = sec_str_from_pdb_file,
#     params                = params,
#     log                   = null_out())
#   print(ssm.records_for_pdb_file())

def run(args):
  pdb_inp = iotbx.pdb.input(file_name = args[0])
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  ssm = mmtbx.secondary_structure.manager(
    pdb_hierarchy=pdb_hierarchy,
    sec_str_from_pdb_file=pdb_inp.extract_secondary_structure(),
    log=null_out())
  helix = ssm.records_for_pdb_file()
  print(helix)



# def write_xls(data):
#   workbook = xwlt.Workbook()
#   worksheet = workbook.add_sheet('PDB')
#   font = xlwt.Font()
#   font.name = 'Times New Roman'
#   for

if __name__ == '__main__':
  run(args=sys.argv[1:])