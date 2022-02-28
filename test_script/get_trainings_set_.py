from __future__ import absolute_import, division, print_function
import sys, os
import mmtbx.secondary_structure
from scitbx.array_family import flex
from libtbx.utils import null_out
import iotbx.pdb
from libtbx import group_args
from mmtbx.validation import ramalyze
from mmtbx_validation_ramachandran_ext import rama_eval
from mmtbx.conformation_dependent_library import generate_protein_threes

fmt_helix="HELIX%5d%4d %3s %s%5d  %3s %s%5s%3s                                 %3d"

def get_ss(hierarchy,
           sec_str_from_pdb_file=None,
           method="ksdssp",
           use_recs=False):
  if(use_recs): 
    params = None
    method = None
  else:
    params = mmtbx.secondary_structure.manager.get_default_ss_params()
    params.secondary_structure.protein.search_method=method
    params = params.secondary_structure
  ssm = mmtbx.secondary_structure.manager(
    pdb_hierarchy         = hierarchy,
    sec_str_from_pdb_file = sec_str_from_pdb_file,
    params                = params,
    log                   = null_out())
  return ssm.helix_selection()
  
def as_list_of_helies(array):
  result = []
  i=0
  while i < array.size():
    vi = array[i]
    if(not vi): 
      i+=1
      continue
    j = i
    v2 = True
    sel = flex.size_t()
    while v2:
      sel.append(j)
      j+=1
      if(j>=array.size()): break
      v2 = array[j]
    result.append(sel)
    i = j
  return result
  
def reject_or_trim(helix, hierarchy):
  # Ramachandran plot
  favored = ramalyze.RAMALYZE_FAVORED
  allowed = ramalyze.RAMALYZE_ALLOWED
  outlier = ramalyze.RAMALYZE_OUTLIER
  re = rama_eval()
  hh = hierarchy.select(helix)
  for three in generate_protein_threes(hierarchy=hh, geometry=None):
    rc = three.get_phi_psi_atoms()
    assert rc is not None
    rama_key = three.get_ramalyze_key()
    angles = three.get_phi_psi_angles()
    rama_score = re.get_score(rama_key, angles[0], angles[1])
    r_eval = re.evaluate_score(rama_key, rama_score)
    phi_atoms, psi_atoms = rc
    i_seqs = [atom.i_seq for atom in phi_atoms] + [psi_atoms[-1].i_seq]
    resnames = three.get_resnames()
    r_name = resnames[1]
    assert rama_key in range(6)
    text_rama_key = ramalyze.res_types[rama_key]
    assert text_rama_key in ["general", "glycine", "cis-proline",
      "trans-proline", "pre-proline", "isoleucine or valine"]
    
    print()
    print (rama_key)
    print (angles)
    print (rama_score)
    print (r_eval)
    print (resnames)
    print (text_rama_key)
    if(  r_eval is favored): print("favored")
    elif(r_eval is allowed): print("allowed")
    elif(r_eval is outlier): print("outlier")
    print()
  #
  return helix
  
def extend(hierarchy, selections):
  asc = hierarchy.atom_selection_cache()
  atoms = hierarchy.atoms()
  for i, selection in enumerate(selections):
    start = get_ids(atoms=atoms, i=selection[0])
    stop  = get_ids(atoms=atoms, i=selection[-1])
    sel_str1 = "chain %s and resseq %d"%(start.cid, (start.resseq-1))
    sel_str2 = "chain %s and resseq %d"%(stop.cid, (stop.resseq+1))
    for sel_str in [sel_str1, sel_str2]:
      sel = asc.selection(sel_str).iselection()
      selection.extend(sel)
    selection = selection.select(flex.sort_permutation(selection))
    selections[i] = selection
  return selections

def get_ids(atoms, i):
  ai = atoms[i]
  return group_args(
    resname = ai.parent().resname,
    resseq  = ai.parent().parent().resseq_as_int(),
    cid     = ai.parent().parent().parent().id)

def run(args):
  # Inputs
  file_name = args[0]
  assert os.path.isfile(file_name)
  pdb_inp = iotbx.pdb.input(file_name=file_name)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  pdb_hierarchy.atoms().reset_i_seq()
  atoms = pdb_hierarchy.atoms()
  sec_str_from_pdb_file = pdb_inp.extract_secondary_structure()
  # SS annotations
  from_header = get_ss(
    hierarchy             = pdb_hierarchy,
    sec_str_from_pdb_file = sec_str_from_pdb_file)
  from_ca     = get_ss(hierarchy = pdb_hierarchy, method = "from_ca")
  from_ksdssp = get_ss(hierarchy = pdb_hierarchy, method = "ksdssp")
  from_cablam = get_ss(hierarchy = pdb_hierarchy, method = "cablam")
  # Merged SS annotation
  merged_sel  = from_header | from_ca | from_ksdssp | from_cablam
  from_all    = as_list_of_helies(merged_sel)
  # Artificially extend helices at both ends so one can elavluate Rama plot
  from_all = extend(hierarchy = pdb_hierarchy, selections = from_all)
  # Analyze and massage helices
  helix_recs = []
  cntr = 0
  for helix in from_all:
    #helix = reject_or_trim(helix = helix, hierarchy = pdb_hierarchy)
    if(helix is None): continue
    cntr += 1
    start = get_ids(atoms=atoms, i=helix[0])
    stop  = get_ids(atoms=atoms, i=helix[-1])
    length = stop.resseq-start.resseq+1
    if(length<4): continue
    h_rec = fmt_helix%(
      cntr,cntr, 
      start.resname, start.cid, start.resseq,
      stop.resname,  stop.cid,  stop.resseq, 
      cntr, length)
    helix_recs.append(h_rec)
  helix_recs_str = "\n".join(helix_recs)
  print(helix_recs_str)
  with open("%s_newSS.pdb"%file_name[:-4],"w") as fo:
    fo.write(helix_recs_str)
    fo.write("\n")
    fo.write(pdb_hierarchy.as_pdb_string())    

if __name__ == '__main__':
  run(args=sys.argv[1:])