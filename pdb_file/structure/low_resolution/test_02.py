import iotbx.pdb

def run():
  pdb_inp = iotbx.pdb.input('1an0.pdb')
  hierarchy = pdb_inp.construct_hierarchy()
  asc = hierarchy.atom_selection_cache()
  sso = pdb_inp.extract_secondary_structure()
  for i, helix in enumerate(sso.helices):
    sel_strs = helix.as_atom_selections()
    if(len(sel_strs) == 1):
      sel = asc.selection(string=sel_strs[0])
    else:
      sel = asc.selection(string=" or ".join( ["(%s)"%s for s in sel_strs] ))
    h_selected = hierarchy.select(sel)
    h_selected.write_pdb_file(file_name='helix_%d.pdb'%i)

if (__name__ == "__main__"):
  run()