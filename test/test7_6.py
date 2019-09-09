import iotbx
import time
from iotbx import pdb
from mmtbx import *
from mmtbx.monomer_library import server
from mmtbx.monomer_library import pdb_interpretation
from mmtbx.conformation_dependent_library import generate_protein_fragments

hierarchy = pdb_inp.construct_hierarchy()
mon_lib_srv = server.server()
ener_lib = server.ener_lib()
processed_pdb = pdb_interpretation.process(self.mon_lib_srv, self.ener_lib, file_name=pdb_file_name,
                                                    raw_records=raw_records)
geometry_restraints_manager = processed_pdb.geometry_restraints_manager()

  # generate five ca
def generate_ca(self):
  t0 = time.time()
  self.hierarchy.reset_atom_i_seqs()
  self.hierarchy.reset_i_seq_if_necessary()
  for five in generate_protein_fragments(self.hierarchy,
                                                  self.geometry_restraints_manager, length=5):
    rc = []
    for atom in five.atoms():
      if atom.name == ' CA ':
        rc.append(atom)
    if len(rc) == 5:
      yield rc