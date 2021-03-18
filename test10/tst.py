from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME mmtbx.development.aev
import sys
import time
import mmtbx
import iotbx.pdb
import mmtbx.model
import copy
import os
from libtbx.utils import null_out
from scitbx.array_family import flex
import __init__ as aev

def run(filename):
  i = 0
  t0 = time.time()
  pdb_inp = iotbx.pdb.input(file_name = filename)
  model = mmtbx.model.manager(
    model_input   = pdb_inp,
    log           = null_out())
  a = aev.AEV(model=model)
  b = aev.compare(a)
  recs = aev.format_HELIX_records_from_AEV(b)
  for pro_list in recs[1]:
    model1 = copy.copy(model)
    chain, num1, num2 = pro_list
    string1 = "chain " + chain +" and resseq " + str(num1) +":"+str(num2)
    print(string1)
    sel = model1.selection(string=string1)
    model1 = model1.select(selection=sel)
    pdbname = str(os.path.basename(filename))[:-4]+"_"+str(i) + ".pdb"
    model1.get_hierarchy().write_pdb_file(file_name='/home/bob/Project/datagenerate/test9/pdb/'+pdbname)
    i += 1
  print(i)
  print("\n".join(recs[0]))
  print('time', time.time()-t0)

if __name__ == '__main__':
  for arg, res, err_str in easy_mp.multi_core_run(run(*tuple(sys.argv[1:])), args, nproc):
    if err_str:
      print(err_str)

