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
  v5 = get_AEV_result(pdb_inp, cut_num=[0.96, 5])
  value1 = flex.linear_correlation(x=v5, y=v2).coefficient()
  value2 = flex.linear_correlation(x=v5, y=v3).coefficient()
  value3 = flex.linear_correlation(x=v5, y=v4).coefficient()
  value4 = flex.linear_correlation(x=v1, y=v5).coefficient()
  return value1, value2, value3, value4


def run():
  from_ca = []
  cablam = []
  ksdssp=[]
  AEV_method = []
  x = []
  t0 = time.time()
  # path = '/home/bob/datagenerate/test11/tst_pdb'
  path = '/home/bob/datagenerate/test11/oringal'

  for root, dirs, files in os.walk(path):
    for file in files:
      print(file)
      filename = str(os.path.join(root, file))
      if 'pdb' in filename:
        result = cal(filename)
        from_ca.append(result[0])
        ksdssp.append(result[1])
        cablam.append(result[2])
        AEV_method.append(result[3])
        x.append(file.split('/')[-1].split('.')[0])
    fig_title = "0.96 AEV figure"
    plt.title(fig_title)
    plt.xlabel("pdb name")
    plt.ylabel("CC values")
    plt.ylim((-0.2, 1))
    plt.plot(x, from_ca, 'r', label='from_ca')
    plt.plot(x, ksdssp, 'g', label='ksdssp')
    plt.plot(x, cablam, 'b', label='cablam')
    plt.plot(x, AEV_method, 'y', label='Standard')
    plt.legend()
    plt.show()
    # plt.savefig('/home/bob/datagenerate/test11/figure/%s.png' % fig_title)
    plt.clf()
  print('time', time.time()-t0)


if __name__ == '__main__':
  run()
