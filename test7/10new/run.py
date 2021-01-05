from __future__ import division
import iotbx.pdb
from scitbx.array_family import flex
from libtbx.utils import null_out
import mmtbx.model

def run():
  pdb_inp = iotbx.pdb.input(file_name = "5JPC.pdb")
  mm = mmtbx.model.manager(model_input = pdb_inp, log = null_out)
  sel = mm.selection(string = "protein")
  mm = mm.select(selection = sel)
  hierarchy = mm.get_hierarchy()
  for model in hierarchy.models():
    for chain in model.chains():
      for conformer in chain.conformers():
        print dir(conformer)
        print conformer.altloc, \
              conformer.atoms().size(), \
              conformer.residues_size()

if (__name__ == "__main__"):
  run()


