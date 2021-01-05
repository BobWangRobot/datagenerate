from __future__ import absolute_import, division, print_function
from iotbx import pdb



pdb_inp = pdb.input(file_name="3nir.pdb")
ph = pdb_inp.construct_hierarchy()
i = 0
for chain in ph.chains():
  for conformer in chain.conformers():
    i += 1
    ph.write_pdb_file(file_name=str(i) + ".pdb")
    print(i,conformer)
    for atom in conformer.atoms():
      print(atom.format_atom_record())
