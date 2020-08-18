from __future__ import division
from __future__ import print_function
import os
import ccp4io_adaptbx

def run(args):
  assert len(args) == 0
  ccp4io_adaptbx.mmdb.Manager()
  ccp4io_adaptbx.ssm.XAlignText()
  tst_ssm()
  print("OK")

def tst_ssm():
  import libtbx.load_env
  from iotbx import pdb
  from libtbx.test_utils import approx_equal
  ssm_pdb1 = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/ssm1.pdb",
    test=os.path.isfile)
  if (ssm_pdb1 is None):
    print("Skipping exercise_regression(): input pdb (ssm1.pdb) not available")
    return
  ssm_pdb2 = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/ssm2.pdb",
    test=os.path.isfile)
  if (ssm_pdb2 is None):
    print("Skipping exercise_regression(): input pdb (ssm2.pdb) not available")
    return
  pdb_io1 = pdb.input(file_name=ssm_pdb1)
  pdb_hierarchy1 = pdb_io1.construct_hierarchy()
  pdb_hierarchy1.reset_i_seq_if_necessary()
  pdb_io2 = pdb.input(file_name=ssm_pdb2)
  pdb_hierarchy2 = pdb_io2.construct_hierarchy()
  pdb_hierarchy2.reset_i_seq_if_necessary()

  reference_chain = None
  moving_chain = None

  #get reference_chain
  for chain in pdb_hierarchy1.chains():
    if chain.id.strip() == "H":
      reference_chain = chain
      break

  #get moving chain
  for chain in pdb_hierarchy2.chains():
    if chain.id.strip() == "B":
      moving_chain = chain
      break

  assert reference_chain is not None
  assert moving_chain is not None

  ssm = ccp4io_adaptbx.SecondaryStructureMatching(
          reference=reference_chain,
          moving=moving_chain)
  assert ssm.ssm.n_align == 53
  assert approx_equal(ssm.ssm.rmsd, 4.68912136041)

  ssm_alignment = ccp4io_adaptbx.SSMAlignment.residue_groups(match=ssm)
  assert len(ssm_alignment.stats) == 230

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
