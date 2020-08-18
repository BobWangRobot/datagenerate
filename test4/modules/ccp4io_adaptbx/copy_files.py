from __future__ import division
from __future__ import print_function
import sys, os
op = os.path

def run(args):
  assert len(args) == 1
  source_root = args[0]
  relative_paths = """\
libccp4/data/syminfo.lib
libccp4/data/symop.lib
libccp4/config.h.in
libccp4/ccp4/ccp4_array.c
libccp4/ccp4/ccp4_array.h
libccp4/fortran/ccp4_diskio_f.c
libccp4/ccp4/ccp4_errno.h
libccp4/ccp4/ccp4_file_err.h
libccp4/ccp4/ccp4_fortran.h
libccp4/ccp4/ccp4_general.c
libccp4/ccp4/ccp4_general.h
libccp4/fortran/ccp4_general_f.c
libccp4/ccp4/ccp4_parser.c
libccp4/ccp4/ccp4_parser.h
libccp4/fortran/ccp4_parser_f.c
libccp4/ccp4/ccp4_program.c
libccp4/ccp4/ccp4_program.h
libccp4/ccp4/ccp4_spg.h
libccp4/ccp4/ccp4_sysdep.h
libccp4/ccp4/ccp4_types.h
libccp4/ccp4/ccp4_unitcell.c
libccp4/ccp4/ccp4_unitcell.h
libccp4/fortran/ccp4_unitcell_f.c
libccp4/ccp4/ccp4_utils.h
libccp4/ccp4/ccp4_vars.h
libccp4/ccp4/cmap_accessor.c
libccp4/ccp4/cmap_close.c
libccp4/ccp4/cmap_data.c
libccp4/ccp4/cmap_data.h
libccp4/ccp4/cmap_errno.h
libccp4/ccp4/cmap_header.c
libccp4/ccp4/cmap_header.h
libccp4/ccp4/cmap_labels.c
libccp4/ccp4/cmap_labels.h
libccp4/ccp4/cmap_open.c
libccp4/ccp4/cmap_skew.c
libccp4/ccp4/cmap_skew.h
libccp4/ccp4/cmap_stats.c
libccp4/ccp4/cmap_stats.h
libccp4/ccp4/cmap_symop.c
libccp4/ccp4/cmaplib.h
libccp4/fortran/cmaplib_f.c
libccp4/ccp4/cmaplib_f.h
libccp4/ccp4/cmtzlib.c
libccp4/ccp4/cmtzlib.h
libccp4/fortran/cmtzlib_f.c
libccp4/ccp4/csymlib.c
libccp4/ccp4/csymlib.h
libccp4/fortran/csymlib_f.c
libccp4/ccp4/cvecmat.c
libccp4/ccp4/cvecmat.h
libccp4/fortran/fftlib.f
libccp4/ccp4/library_err.c
libccp4/fortran/library_f.c
libccp4/ccp4/library_file.c
libccp4/ccp4/library_file.h
libccp4/ccp4/library_utils.c
libccp4/ccp4/mtzdata.h
libccp4/ccp4/overview.h
libccp4/ccp4/pack_c.c
libccp4/ccp4/pack_c.h
mmdb/mmdb/hybrid_36.cpp
mmdb/mmdb/hybrid_36.h
mmdb/mmdb/mmdb_atom.cpp
mmdb/mmdb/mmdb_atom.h
mmdb/mmdb/mmdb_bondmngr.cpp
mmdb/mmdb/mmdb_bondmngr.h
mmdb/mmdb/mmdb_chain.cpp
mmdb/mmdb/mmdb_chain.h
mmdb/mmdb/mmdb_cifdefs.cpp
mmdb/mmdb/mmdb_cifdefs.h
mmdb/mmdb/mmdb_coormngr.cpp
mmdb/mmdb/mmdb_coormngr.h
mmdb/mmdb/mmdb_cryst.cpp
mmdb/mmdb/mmdb_cryst.h
mmdb/mmdb/mmdb_defs.h
mmdb/mmdb/mmdb_ficif.cpp
mmdb/mmdb/mmdb_ficif.h
mmdb/mmdb/mmdb_io_file.cpp
mmdb/mmdb/mmdb_io_file.h
mmdb/mmdb/mmdb_io_stream.cpp
mmdb/mmdb/mmdb_io_stream.h
mmdb/mmdb/mmdb_machine_.cpp
mmdb/mmdb/mmdb_machine_.h
mmdb/mmdb/mmdb_manager.cpp
mmdb/mmdb/mmdb_manager.h
mmdb/mmdb/mmdb_mask.cpp
mmdb/mmdb/mmdb_mask.h
mmdb/mmdb/mmdb_math_align.cpp
mmdb/mmdb/mmdb_math_align.h
mmdb/mmdb/mmdb_math_bfgsmin.cpp
mmdb/mmdb/mmdb_math_bfgsmin.h
mmdb/mmdb/mmdb_math_.cpp
mmdb/mmdb/mmdb_math_fft.cpp
mmdb/mmdb/mmdb_math_fft.h
mmdb/mmdb/mmdb_math_graph.cpp
mmdb/mmdb/mmdb_math_graph.h
mmdb/mmdb/mmdb_math_.h
mmdb/mmdb/mmdb_math_linalg.cpp
mmdb/mmdb/mmdb_math_linalg.h
mmdb/mmdb/mmdb_math_rand.cpp
mmdb/mmdb/mmdb_math_rand.h
mmdb/mmdb/mmdb_mattype.cpp
mmdb/mmdb/mmdb_mattype.h
mmdb/mmdb/mmdb_mmcif_.cpp
mmdb/mmdb/mmdb_mmcif_.h
mmdb/mmdb/mmdb_model.cpp
mmdb/mmdb/mmdb_model.h
mmdb/mmdb/mmdb_root.cpp
mmdb/mmdb/mmdb_root.h
mmdb/mmdb/mmdb_rwbrook.cpp
mmdb/mmdb/mmdb_rwbrook.h
mmdb/mmdb/mmdb_selmngr.cpp
mmdb/mmdb/mmdb_selmngr.h
mmdb/mmdb/mmdb_seqsuperpose.cpp
mmdb/mmdb/mmdb_seqsuperpose.h
mmdb/mmdb/mmdb_symop.cpp
mmdb/mmdb/mmdb_symop.h
mmdb/mmdb/mmdb_tables.cpp
mmdb/mmdb/mmdb_tables.h
mmdb/mmdb/mmdb_title.cpp
mmdb/mmdb/mmdb_title.h
mmdb/mmdb/mmdb_uddata.cpp
mmdb/mmdb/mmdb_uddata.h
mmdb/mmdb/mmdb_utils.cpp
mmdb/mmdb/mmdb_utils.h
mmdb/mmdb/mmdb_xml_.cpp
mmdb/mmdb/mmdb_xml_.h
ssm/ssm/ssm_defs.h
ssm/ssm/ssm_csia.cpp
ssm/ssm/ssm_csia.h
ssm/ssm/ssm_graph.cpp
ssm/ssm/ssm_graph.h
ssm/ssm/ssm_vxedge.cpp
ssm/ssm/ssm_vxedge.h
ssm/ssm/ssm_align.cpp
ssm/ssm/ssm_align.h
ssm/ssm/ssm_superpose.cpp
ssm/ssm/ssm_superpose.h
""".splitlines()
  n_makedirs = 0
  n_copied = 0
  n_updated = 0
  n_already_up_to_date = 0
  for path in relative_paths:
    dir = op.split(path)[0]
    if (dir.startswith("mmdb/")) :
      dir = dir[5:]
    elif (dir.startswith("ssm/")) :
      dir = dir[4:]
    if (dir != "" and not op.isdir(dir)):
      os.makedirs(dir)
      n_makedirs += 1
    source = open(op.join(source_root, path), "rb").read() \
      .replace(
        "#undef PACKAGE_ROOT",
        "#define PACKAGE_ROOT NULL")
    tpath = path
    if (tpath.startswith("mmdb/")) :
      tpath = tpath[5:]
    elif (tpath.startswith("ssm/")) :
      tpath = tpath[4:]
    if (tpath.endswith(".in")):
      tpath = tpath[:-3]
    if (not op.isfile(tpath)):
      n_copied += 1
    else:
      target = open(tpath, "rb").read()
      if (target == source):
        n_already_up_to_date += 1
        source = None
      else:
        n_updated += 1
    if (source is not None):
      open(tpath, "wb").write(source)
  print("Directories created:", n_makedirs)
  print("Files copied:", n_copied)
  print("Files updated:", n_updated)
  print("Files already up-to-date:", n_already_up_to_date)
  print("Done.")

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
