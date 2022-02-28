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
import matplotlib.pyplot as plt
import os
import numpy as np

old_method = [0.8562308024308593, 0.9510228411791382, 0.07121937990402431, 0.4574000819554095, 0.7248725998408653,
       0.8111345584607443, 0.7172813119270387, 0.2570783734521071, 0.10535485055837265]

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
  recs = aev.format_HELIX_records_from_AEV(CC_value, cut_num)
  return recs

def cal(filename):
  pdb_inp = iotbx.pdb.input(file_name = filename)
  hierarchy = pdb_inp.construct_hierarchy()
  hierarchy.remove_alt_confs(always_keep_one_conformer=True)
  ase = hierarchy.atom_selection_cache().selection('protein')
  pdb_hierarchy = hierarchy.select(ase)
  v1 = get_ss(
    hierarchy=pdb_hierarchy,
    sec_str_from_pdb_file=pdb_inp.extract_secondary_structure())
  v2 = get_ss(hierarchy=pdb_hierarchy, method="from_ca")
  v3 = get_ss(hierarchy=pdb_hierarchy, method="ksdssp")
  v4 = get_ss(hierarchy=pdb_hierarchy, method="cablam")
  v5 = get_AEV_result(pdb_inp, cut_num=[0.98, 6])
  return v1, v2, v3, v4, v5

def run():
  records = flex.double()
  AEV_method = flex.double()
  from_ca = flex.double()
  ksddsp = flex.double()
  cablam = flex.double()
  t0 = time.time()
  path = '/home/bob/datagenerate/test11/modified'
  for root, dirs, files in os.walk(path):
    for file in files:
      print(file)
      filename = str(os.path.join(root, file))
      if 'pdb' in filename:
        value_result = cal(filename)
        records.extend(value_result[0])
        from_ca.extend(value_result[1])
        ksddsp.extend(value_result[2])
        cablam.extend(value_result[3])
        AEV_method.extend(value_result[4])
  from_ca = flex.linear_correlation(x=records, y=from_ca).coefficient()
  ksddsp = flex.linear_correlation(x=records, y=ksddsp).coefficient()
  cablam = flex.linear_correlation(x=records, y=cablam).coefficient()
  AEV_method = flex.linear_correlation(x=records, y=AEV_method).coefficient()
  print("from_ca:", from_ca, "ksddsp:", ksddsp, "cablam: ", cablam, "AEV_method:", AEV_method)
  print('time', time.time() - t0)


if __name__ == '__main__':
  run()
