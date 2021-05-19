from __future__ import absolute_import, division, print_function
import sys
import mmtbx
import iotbx.pdb
from libtbx.utils import null_out
from scitbx.array_family import flex
import mmtbx.secondary_structure

def match_score(x,y):
  assert x.size() == y.size()
  match_cntr = 0
  for x_,y_ in zip(x,y):
    if(x_==y_): match_cntr+=1
  return match_cntr/x.size()

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
  #print("- %s"%method*10)
  #print(ssm.records_for_pdb_file())
  #print("-"*79)
  alpha = ssm.helix_selection()
  CA = ssm.selection_cache.selection('protein')
  for x_, y_ in zip(alpha, CA):
    if x_ == True and y_ == True:
      list1.append(1)
    elif y_ == True:
      list1.append(0)
  print(list(list1))
  return list1

def run(args):
  pdb_inp = iotbx.pdb.input(file_name = args[0])
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  # ss = hierarchy.atom_selection_cache().selection('name CA')
  # pdb_hierarchy = hierarchy.select(ss)
  # get secodary structure annotation vector from HELIX/SHEET records (file header)
  v1 = get_ss(
    hierarchy             = pdb_hierarchy,
    sec_str_from_pdb_file = pdb_inp.extract_secondary_structure())
  v2 = get_ss(hierarchy = pdb_hierarchy, method = "from_ca")
  v3 = get_ss(hierarchy = pdb_hierarchy, method = "ksdssp")
  print("CC:", flex.linear_correlation(x = v1, y = v2).coefficient())
  print("CC:", flex.linear_correlation(x = v1, y = v3).coefficient())
  print("CC:", flex.linear_correlation(x = v3, y = v2).coefficient())
  print("match:", match_score(x = v1, y = v2))
  print("match:", match_score(x = v1, y = v3))
  print("match:", match_score(x = v3, y = v2))
  
if __name__ == '__main__':
  run(args=sys.argv[1:])
