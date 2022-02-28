from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME mmtbx.development.aev
import sys
import time
import mmtbx
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from scitbx.array_family import flex
import __init__ as aev
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
  alpha = ssm.helix_selection()
  CA = ssm.selection_cache.selection('name CA')
  for x_, y_ in zip(alpha, CA):
    if x_ == True and y_ == True:
      list1.append(1)
    elif y_ == True:
      list1.append(0)
  return list1

def get_AEV_result(pdb_inp,cut_num):
  model = mmtbx.model.manager(
    model_input=pdb_inp,
    log=null_out())
  # sel = model.selection(string="protein")
  # model = model.select(selection=sel)
  model.crystal_symmetry()
  a = aev.AEV(model=model)
  CC_value = aev.compare(a)
  recs = aev.format_HELIX_records_from_AEV(CC_value,cut_num)
  return recs

def run(args):
  t0 = time.time()
  pdb_inp = iotbx.pdb.input(file_name = args[0])
  hierarchy = pdb_inp.construct_hierarchy()
  hierarchy.remove_alt_confs(always_keep_one_conformer=True)
  ase = hierarchy.atom_selection_cache().selection('protein')
  pdb_hierarchy = hierarchy.select(ase)
  # v1 = get_AEV_result(pdb_inp, cut_num=[0.9, 4])
  v2 = get_ss(hierarchy=pdb_hierarchy, method="from_ca")
  v3 = get_ss(hierarchy=pdb_hierarchy, method="ksdssp")
  v4 = get_ss(hierarchy=pdb_hierarchy, method="cablam")
  for i in range(73, 100):
    thre = i/100
    v1 = get_AEV_result(pdb_inp, cut_num=[thre, 4])
    print("threshold:", i/100)
    print("CC of from_ca:", flex.linear_correlation(x=v1, y=v2).coefficient())
    # print("match of from_ca:", match_score(x=v1, y=v2))
    print("CC of ksdssp:", flex.linear_correlation(x=v1, y=v3).coefficient())
    # print("match of ksdssp:", match_score(x=v1, y=v3))
    print("CC of cablam:", flex.linear_correlation(x=v1, y=v4).coefficient())
    # print("match of cablam:", match_score(x=v1, y=v4))
    # print(list(v4))
    print("CC of standard", flex.linear_correlation(x=v1, y=v4).coefficient())
  # print("CC of from_ca:", flex.linear_correlation(x=v1, y=v2).coefficient())
  # print("match of from_ca:", match_score(x=v1, y=v2))
  # print("CC of ksdssp:", flex.linear_correlation(x=v2, y=v4).coefficient())
  # print("match of ksdssp:", match_score(x=v2, y=v3))
  print('time', time.time()-t0)

if __name__ == '__main__':
  run(args=sys.argv[1:])
