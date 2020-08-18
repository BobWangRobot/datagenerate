import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from scitbx.array_family import flex
from mmtbx.dynamics import cartesian_dynamics

def run_dynamics(xrs, gc, stop_at_diff):
  o = cartesian_dynamics.run(
    xray_structure                   = xrs,
    gradients_calculator             = gc,
    temperature                      = 5000,
    n_steps                          = 10000,
    time_step                        = 0.0005,
    initial_velocities_zero_fraction = 0,
    vxyz                             = None,
    n_print                          = 20,
    n_collect                        = 10,
    interleaved_minimization         = False,
    reset_velocities                 = True,
    stop_cm_motion                   = False,
    log                              = None,
    stop_at_diff                     = stop_at_diff,
    random_seed                      = None,
    states_collector                 = None,
    verbose                          = -1)
  return xrs.sites_cart()

def run(file_name):
  pdb_inp = iotbx.pdb.input(file_name = file_name)
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  #
  p.pdb_interpretation.nonbonded_weight=100
  #
  p.pdb_interpretation.c_beta_restraints=False
  p.pdb_interpretation.restraints_library.cdl=False
  p.pdb_interpretation.peptide_link.ramachandran_restraints=False
  
  p.pdb_interpretation.use_neutron_distances=True
  p.pdb_interpretation.secondary_structure.enabled=False
  #
  p.pdb_interpretation.c_beta_restraints
  model = mmtbx.model.manager(
    model_input               = pdb_inp,
    pdb_interpretation_params = p,
    build_grm                 = True,
    log                       = null_out())
  #
  sites_cart_start = model.get_sites_cart()
  gc = cartesian_dynamics.gradients_calculator_geometry_restraints(
    restraints_manager = model.get_restraints_manager())
  for it in range(0,100):
    root = iotbx.pdb.hierarchy.root()
    sites_cart = run_dynamics(
      xrs          = model.get_xray_structure().deep_copy_scatterers(), 
      gc           = gc, 
      stop_at_diff = 5)
    dist = flex.mean(flex.sqrt((sites_cart_start - sites_cart).dot()))
    print "trial: %3d rmsd: %5.3f"%(it, dist)
    ph = model.get_hierarchy().deep_copy()
    ph.atoms().set_xyz(sites_cart)
    #
    chain = list(ph.chains())[0]
    m = iotbx.pdb.hierarchy.model(id=str(it))
    chain_ = chain.detached_copy()
    chain_.id="A"
    m.append_chain(chain_)
    root.append_model(m)
    root.write_pdb_file(
      file_name        = "%s.pdb"%it,
      crystal_symmetry = pdb_inp.crystal_symmetry())

if (__name__ == "__main__"):
  run(file_name="polyGly-helix.pdb")

