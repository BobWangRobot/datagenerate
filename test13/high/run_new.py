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

def trim(residue, hierarchy, helix):
  # delete abnormal residues of a helix
  asc = hierarchy.atom_selection_cache()
  atoms = hierarchy.atoms()
  atom_id = residue.atoms()[-1].i_seq
  target = get_ids(atoms=atoms, i=atom_id)
  sel_str = "chain %s and resseq %d" % (target.cid, (target.resseq))
  sel = asc.selection(sel_str).iselection()
  helix = flex.size_t(list(set(helix) - set(sel)))
  helix = helix.select(flex.sort_permutation(helix))
  return helix

def reject_or_trim(helix, hierarchy):
  # Ramachandran plot
  favored = ramalyze.RAMALYZE_FAVORED
  allowed = ramalyze.RAMALYZE_ALLOWED
  outlier = ramalyze.RAMALYZE_OUTLIER
  re = rama_eval()
  hh = hierarchy.select(helix)
  double_delete = 0
  altloc_dict = {}
  for three in generate_protein_threes(hierarchy=hh, geometry=None):
    rc = three.get_phi_psi_atoms()
    # assert rc is not None
    if rc is not None:
      rama_key = three.get_ramalyze_key()
      angles = three.get_phi_psi_angles()
      resseq = rc[0][0].parent().parent().resseq_as_int()
      altloc_dict.setdefault(resseq, 0)
      # print(resseq)
      rama_score = re.get_score(rama_key, angles[0], angles[1])
      r_eval = re.evaluate_score(rama_key, rama_score)
      phi_atoms, psi_atoms = rc
      # i_seqs = [atom.i_seq for atom in phi_atoms] + [psi_atoms[-1].i_seq]
      resnames = three.get_resnames()
      # r_name = resnames[1]
      # print(three)
      # print("phi and psi:", angles[0], angles[1])
      assert rama_key in range(6)
      text_rama_key = ramalyze.res_types[rama_key]
      assert text_rama_key in ["general", "glycine", "cis-proline",
        "trans-proline", "pre-proline", "isoleucine or valine"]
      if not (-84 < angles[0] < -53 and -53 < angles[1] < -25):
        double_delete += 1
        altloc_dict[resseq] += 1
        if double_delete > 1 and altloc_dict[resseq]==2:
          helix = trim(three[1], hierarchy, helix)
          # print(three)
          # print("phi and psi:", angles[0], angles[1])
          # print("delete")
          double_delete = 0
      if(  r_eval is favored):
        pass
      elif(r_eval is allowed):
        pass
        # helix = trim(three[1], hierarchy, helix)
      elif(r_eval is outlier) :
        # print("bad")
        helix = trim(three[1], hierarchy, helix)
  return helix

def check_continuity(helix, hierarchy):
  # checing continuity of a helix, cut the helix if it is not continue.
  atoms = hierarchy.atoms()
  atom_id = helix[0]
  start = get_ids(atoms=atoms, i=atom_id).resseq
  helix_sub = []
  new_helix = []
  assert start == get_ids(atoms=atoms, i=helix[0]).resseq
  for id in helix:
    atom_reseeq = get_ids(atoms=atoms, i=id).resseq
    if atom_reseeq != start:
      start += 1
      if atom_reseeq != start:
        start = atom_reseeq
        new_helix.append(flex.size_t(helix_sub))
        helix_sub = []
    helix_sub.append(id)
  if helix_sub not in new_helix:
    new_helix.append(flex.size_t(helix_sub))
  return new_helix

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
  simialrity = 0
  all_atoms = 0
  for atom1, atom2 in zip(from_ca, from_header):
    all_atoms += 1
    if atom1==atom2:
      simialrity += 1
  sim = simialrity/all_atoms
  print(sim)
  # Merged SS annotation
  merged_sel  = from_header | from_ca | from_ksdssp | from_cablam
  from_all    = as_list_of_helies(merged_sel)
  # Artificially extend helices at both ends so one can elavluate Rama plot
  from_all = extend(hierarchy = pdb_hierarchy, selections = from_all)
  # Analyze and massage helices
  helix_recs = []
  cntr = 0
  for helix in from_all:
    helix = reject_or_trim(helix = helix, hierarchy = pdb_hierarchy)
    if(helix is None): continue
    helix = check_continuity(helix=helix, hierarchy=pdb_hierarchy)
    for new_helix in helix:
      cntr += 1
      start = get_ids(atoms=atoms, i=new_helix[0])
      stop  = get_ids(atoms=atoms, i=new_helix[-1])
      length = stop.resseq-start.resseq+1
      if(length<6):
        cntr = cntr - 1
        continue
      h_rec = fmt_helix%(
        cntr, cntr,
        start.resname, start.cid, start.resseq,
        stop.resname,  stop.cid,  stop.resseq,
        1, length)
      helix_recs.append(h_rec)
  helix_recs_str = "\n".join(helix_recs)
  print(helix_recs_str)
  # print(pdb_hierarchy.as_pdb_string())
  with open("%s_newSS.pdb"%file_name[:-4],"w") as fo:
    fo.write(helix_recs_str)
    fo.write("\n")
    fo.write(pdb_hierarchy.as_pdb_string(crystal_symmetry=pdb_inp.crystal_symmetry()))

if __name__ == '__main__':
  run(args=sys.argv[1:])