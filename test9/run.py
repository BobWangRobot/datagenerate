import time
import mmtbx
import iotbx.pdb
import mmtbx.model
import os
import sys
from libtbx.utils import null_out
from scitbx.array_family import flex
import __init__ as aev

def main(pdb_dir):
  t0 = time.time()
  for root, dirs, files in os.walk(pdb_dir):
    for file in files:
      try:
        filename = str(os.path.join(root, file))
        print(filename)
        pdb_inp = iotbx.pdb.input(file_name=filename)
        model = mmtbx.model.manager(
          model_input=pdb_inp,
          log=null_out())
        sel = model.selection(string="protein")
        model = model.select(selection=sel)
        a = aev.AEV(model=model)
        b = aev.compare(a)
        recs = aev.format_HELIX_records_from_AEV(b)
        print("\n".join(recs))
      except:
        pass
      print('time', time.time() - t0)

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))


