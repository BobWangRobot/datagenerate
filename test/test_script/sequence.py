import iotbx.pdb
import os, sys

def main(filename):
  pdb_inp = iotbx.pdb.pdb_input(filename)
  hir = pdb_inp.construct_hierarchy()
  hir.reset_atom_i_seqs()
  for atom in hir.atoms():
    if ' CA ' in atom.id_str():
      print(atom.id_str())
  #hir.write_pdb_file(file_name="%s"%filename.replace('.pdb','re.pdb'))

if __name__=='__main__':
  main(*tuple(sys.argv[1:]))