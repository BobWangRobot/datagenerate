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
  mm = mmtbx.model.manager(
    model_input   = pdb_inp,
    log           = null_out())
  sel = mm.selection(string="protein")
  mm = mm.select(selection=sel)
  hierarchy = mm.get_hierarchy()
  for model in hierarchy.models():
    for chain in model.chains():
      for conformer in chain.conformers():
        print(conformer.altloc,
        conformer.atoms().size(),
        conformer.residues_size())
  a = aev.AEV(model = mm)
  print(a.BAEVs)
  b = aev.compare(a)
  print(b)
  recs = aev.format_HELIX_records_from_AEV(b)
  print("\n".join(recs))
  print('time', time.time()-t0)

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))
