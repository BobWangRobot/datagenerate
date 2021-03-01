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

def main(filename):
  t0 = time.time()
  pdb_inp = iotbx.pdb.input(file_name = filename)
  model = mmtbx.model.manager(
    model_input   = pdb_inp,
    log           = null_out())
  a = aev.AEV(model=model)
  b = aev.compare(a)
  recs = aev.format_HELIX_records_from_AEV(b)
  num1 = recs[0].split()[-5]
  num2 = recs[0].split()[-2]
  string1 = "resseq " + num1 +":"+num2
  sel = model.selection(string=string1)
  model = model.select(selection=sel)
  hi = model.get_hierarchy()
  # print(dir(hi.write_pdb_file(file_name='1.pdb')))
  model.get_hierarchy().write_pdb_file(file_name='1.pdb')
  print("\n".join(recs))
  print('time', time.time()-t0)

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))
