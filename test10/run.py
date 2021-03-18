import time
import mmtbx
import iotbx.pdb
import mmtbx.model
import os
import sys
import copy
from libtbx.utils import null_out
from scitbx.array_family import flex
import __init__ as aev

def main(pdb_dir):
  t0 = time.time()
  for root, dirs, files in os.walk(pdb_dir):
    for file in files:
      print(file)
      try:
        filename = str(os.path.join(root, file))
        pdb_inp = iotbx.pdb.input(file_name=filename)
        model = mmtbx.model.manager(
          model_input=pdb_inp,
          log=null_out())
        a = aev.AEV(model=model)
        b = aev.compare(a)
        recs = aev.format_HELIX_records_from_AEV(b)
        for pro_list in recs[1]:
          model1 = copy.copy(model)
          chain, num1, num2 = pro_list
          string1 = "chain " + chain + " and resseq " + str(num1) + ":" + str(num2)
          print(string1)
          sel = model1.selection(string=string1)
          model1 = model1.select(selection=sel)
          pdbname = str(os.path.basename(filename))[:7] + "_" + string1.replace(' ', '_').replace(':', '_') + ".pdb"
          model1.get_hierarchy().write_pdb_file(file_name='/mnt/hgfs/pdb/PRO_helices/' + pdbname)
      except:
        pass
      print('time', time.time() - t0)

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))


