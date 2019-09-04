import iotbx.pdb

pdb_inp = iotbx.pdb.input('helix2.pdb')
hierarchy = pdb_inp.construct_hierarchy()
hierarchy.reset_atom_i_seqs()
hierarchy.write_pdb_file(file_name='helix3.pdb')