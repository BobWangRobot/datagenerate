from __future__ import division
from __future__ import print_function
import boost.python
ext = boost.python.import_ext("ccp4io_adaptbx_ext")
from ccp4io_adaptbx_ext import *

import libtbx.load_env
import operator
import os
op = os.path

for _ in ["libccp4/data/syminfo.lib",
          "lib/data/syminfo.lib"]:
  _ = libtbx.env.under_dist("ccp4io", _)
  if (op.isfile(_)):
    os.environ["SYMINFO"] = _
    break

else:
  if os.getenv( "CLIBD" ) and op.isfile( op.join( os.getenv( "CLIBD" ), "syminfo.lib" ) ):
    os.environ["SYMINFO"] = op.join( os.getenv( "CLIBD" ), "syminfo.lib" )

  else:
    import warnings
    warnings.warn("ccp4io_adaptbx: cannot locate syminfo.lib")

def to_mmdb(root, flags = [], remove_duplicates = True):
  "Converts iotbx.pdb.hierarchy object to an MMDB Manager"

  manager = mmdb.Manager()

  if remove_duplicates:
    selector = select_atoms_with_same_resname

  else:
    selector = lambda rg: rg.atoms()

  if flags:
    manager.SetFlag( reduce( operator.or_, flags ) )

  for rg in root.residue_groups():
    for atom in selector( rg = rg ):
      rc = manager.PutPDBString( atom.format_atom_record() )

      if 0 < rc:
        raise RuntimeError(mmdb.GetErrorDescription( rc ))

  return manager


def select_atoms_with_same_resname(rg):

  atoms_in = {}

  for a in rg.atoms():
    atoms_in.setdefault( a.parent().resname, [] ).append( a )

  return max([(len(_),_) for _ in atoms_in.values()])[1]


ssm.ERROR_DESCRIPTION_FOR = {
  ssm.RETURN_CODE.NoHits: "secondary structure does not match",
  ssm.RETURN_CODE.NoSuperposition: "structures are too remote",
  ssm.RETURN_CODE.NoGraph: "can't make graph for first structure",
  ssm.RETURN_CODE.NoVertices: "empty graph for first structure",
  ssm.RETURN_CODE.NoGraph2: "can't make graph for second structure",
  ssm.RETURN_CODE.NoVertices2: "empty graph for second structure",
  }

ssm.GetErrorDescription = lambda rc: (
  ssm.ERROR_DESCRIPTION_FOR.get( rc, "undocumented return code" )
  )

class SecondaryStructureMatching(object):
  "SSM matching with iotbx.pdb.hierarchy objects"

  def __init__(
    self,
    moving,
    reference,
    precision = ssm.PRECISION.Normal,
    connectivity = ssm.CONNECTIVITY.Flexible
    ):

    self.chains = []
    self.managers = []
    self.handles = []
    self.qvalues = []

    for chain in [ moving, reference ]:
      manager = to_mmdb( root = chain )
      handle = manager.NewSelection()
      manager.Select(
        selHnd = handle,
        selType = mmdb.SELECTION_TYPE.ATOM,
        cid = "*",
        selKey = mmdb.SELECTION_KEY.NEW
        )

      if manager.GetSelLength( selHnd = handle ) <= 0:
        raise RuntimeError(
          "Empty atom selection for structure %s" % ( len( self.managers ) + 1 ))

      self.chains.append( chain )
      self.managers.append( manager )
      self.handles.append( handle )

    assert len( self.chains ) == 2
    assert len( self.managers ) == 2
    assert len(self.handles ) == 2

    self.ssm = ssm.SSMAlign()
    rc = self.ssm.Align(
      manager1 = self.managers[0],
      manager2 = self.managers[1],
      precision = precision,
      connectivity = connectivity,
      selHnd1 = self.handles[0],
      selHnd2 = self.handles[1],
      )
    if rc != ssm.RETURN_CODE.Ok:
      raise RuntimeError(ssm.GetErrorDescription( rc = rc ))

    self.qvalues = self.ssm.GetQvalues()



  def GetQvalues(self):
      return self.qvalues


  def AlignSelectedMatch(self, nselected):
    if nselected >= len(self.qvalues):
        print("Not that many matches available")
        return

    rc = self.ssm.AlignSelectedMatch(
      manager1 = self.managers[0],
      manager2 = self.managers[1],
      precision = ssm.PRECISION.Normal,
      connectivity = ssm.CONNECTIVITY.Flexible,
      selHnd1 = self.handles[0],
      selHnd2 = self.handles[1],
      nselected = nselected
      )
    if rc != ssm.RETURN_CODE.Ok:
      raise RuntimeError(ssm.GetErrorDescription( rc = rc ))



  def get_matrix(self):
    return self.ssm.t_matrix



class SSMAlignment(object):
  "SSM alignment from SSM match"

  def __init__(self, match, indexer):

    align = ssm.XAlignText()
    align.XAlign(
      manager1 = match.managers[0],
      manager2 = match.managers[1],
      ssm_align = match.ssm
      )
    indexer1 = indexer( chain = match.chains[0] )
    indexer2 = indexer( chain = match.chains[1] )
    self.pairs = []
    self.stats = []
    self.stats2 = []

    for ( ( f, s ), a ) in align.get_blocks():
      def get(rgi, indexer):
        if (rgi is None):
            return None

        identifier = ( rgi.chain_id, rgi.resseq, rgi.inscode )
        assert identifier in indexer, "Id %s missing" % str( identifier )
        return indexer[ identifier ]

      self.pairs.append( ( get( f, indexer1 ), get( s, indexer2 ) ) )
      self.stats2.append( (f, s) ) # also contains hydropathy, resname, restype
      self.stats.append( a )
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )



  def GetSSMSequenceIdentity(self):
    gaplessalignlen = 0
    s5 = 0.0; s4 = 0.0; s3 = 0.0; s2 = 0.0; s1 = 0.0; s0 = 0.0
    for e in self.stats:
        if len(e) > 0:
            gaplessalignlen += 1
            """ count number of similar residues in alignment according to score
              S/H   residue belongs to a strand/helix
              +/-/. hydrophylic/hydrophobic/neutral residue
              **    identical residues matched: similarity 5
              ++    similarity 4
              ==    similarity 3
              --    similarity 2
              ::    similarity 1
              ..    dissimilar residues: similarity 0
            """
            if e[2] <= 0: s0 += 1.0
            if e[2] == 1: s1 += 1.0
            if e[2] == 2: s2 += 1.0
            if e[2] == 3: s3 += 1.0
            if e[2] == 4: s4 += 1.0
            if e[2] == 5: s5 += 1.0

    glen = gaplessalignlen/100.0
    #       identical sim4      sim3    sim2     sim1    dissimilar
    return (s5/glen, s4/glen, s3/glen, s2/glen, s1/glen, s0/glen, gaplessalignlen)



  def residue_groups(cls, match):

    indexer = lambda chain: dict(
      [ ( ( chain.id, rg.resseq_as_int(), rg.icode ), rg )
        for rg in chain.residue_groups() ]
      )

    return cls( match = match, indexer = indexer )

  residue_groups = classmethod( residue_groups )

ssm.MALIGN_ERROR_DESCRIPTION_FOR = {
  ssm.MALIGN.BadInput: "bad input",
  ssm.MALIGN.NoStructure: "NULL structure",
  ssm.MALIGN.NoAlignment: "multiple alignment was not achieved",
  }

ssm.GetMultAlignErrorDescription = lambda rc: (
  ssm.MALIGN_ERROR_DESCRIPTION_FOR.get( rc, "unknown return code" ) if rc < ssm.MALIGN.NoGraph
  else "can't make graph for %s" % ( rc - ssm.MALIGN.NoGraph )
  )

def reformat_ssm_t_matrix(mx):

  return (
    ( mx[0], mx[1], mx[2], mx[4], mx[5], mx[6], mx[8], mx[9], mx[10] ),
    ( mx[3], mx[7], mx[11] ),
    )


class SSMMultipleAlignment(object):
  """
  A convenience object to hold the results
  """

  def __init__(self, managers, selstrings):

    if len( managers ) != len( selstrings ):
      raise ValueError("Iterables are not the same length")

    multalign = ssm.MultipleAlignment( managers = managers, selstrings = selstrings )

    if multalign.get_return_code() != ssm.MALIGN.Ok:
      raise RuntimeError(ssm.GetMultAlignErrorDescription(
        rc = multalign.get_return_code(),
        ))

    self.alignment = multalign.get_alignment()
    self.t_matrices = multalign.get_matrices()
    self.rmsd = multalign.get_rmsd()


  @property
  def rt_matrices(self):

    return [ reformat_ssm_t_matrix( mx = mx ) for mx in self.t_matrices ]

