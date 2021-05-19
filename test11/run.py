from __future__ import absolute_import, division, print_function
import sys
import mmtbx
import iotbx.pdb
from libtbx.utils import null_out
from scitbx.array_family import flex
import mmtbx.secondary_structure

def get_ss(hierarchy, sec_str_from_pdb_file=None, params=None):
  ssm = mmtbx.secondary_structure.manager(
    pdb_hierarchy         = hierarchy,
    sec_str_from_pdb_file = sec_str_from_pdb_file,
    params                = params,
    log                   = null_out())
  print("-"*79)
  print(ssm.records_for_pdb_file())
  print(len(list(ssm.helix_selection())))
  print("-"*79)
  alpha = ssm.helix_selection()
  beta  = ssm.beta_selection()
  assert alpha.size() == beta.size() == hierarchy.atoms().size()
  annotation_vector = flex.double(hierarchy.atoms().size(), 0)
  annotation_vector.set_selected(alpha, 1)
  annotation_vector.set_selected(beta, 2)
  return annotation_vector

def run(args):
  pdb_inp = iotbx.pdb.input(file_name = args[0])
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  # get secodary structure annotation vector from atomic coordinates
  v1 = get_ss(hierarchy = pdb_hierarchy)
  # get secodary structure annotation vector from HELIX/SHEET records (file header)
  v2 = get_ss(
    hierarchy             = pdb_hierarchy,
    sec_str_from_pdb_file = pdb_inp.extract_secondary_structure())
  print(v2)
  # compare them using correlation
  cc = flex.linear_correlation(x = v1, y = v2).coefficient()
  print(cc)

if __name__ == '__main__':
  run(args=sys.argv[1:])
