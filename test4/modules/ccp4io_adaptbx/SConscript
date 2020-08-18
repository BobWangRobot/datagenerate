from libtbx.str_utils import show_string
import libtbx.load_env
import re
import sys, os
op = os.path

Import("env_base", "env_etc")

env_etc.ccp4io_dist = libtbx.env.dist_path("ccp4io")

# FIXME
# Set lower size limits for MTZ file contents; this is a fix for an obscure
# bug which causes Phenix to crash on MacOS 10.5.  Probably not necessary in
# all circumstances though.
if env_etc.compiler == "win32_cl":
  # XXX for reasons that are unclear to me, defining CALL_LIKE_SUN is necessary
  # to get the proper function names in ccp4_fortran.h using VC++
  env_etc.ccp4io_defines = ["/Di386", "/D_MVS", "/DMXTALS=32", "/DMSETS=32",
    "/DMCOLUMNS=128", "/DCALL_LIKE_SUN"]
else:
  env_etc.ccp4io_defines = ["-DMSETS=32", "-DMXTALS=32", "-DMCOLUMNS=128"]

path_lib_src = op.join(env_etc.ccp4io_dist, "libccp4", "ccp4")
ccp4_src = op.join("libccp4", "ccp4")
env_etc.ccp4io_include = libtbx.env.under_dist(
  module_name="ccp4io", path="libccp4/ccp4")

build_ccp4io_adaptbx = libtbx.env.under_build("ccp4io_adaptbx")
if (not op.isdir(build_ccp4io_adaptbx)):
  os.mkdir(build_ccp4io_adaptbx)
  assert op.isdir(build_ccp4io_adaptbx)

def replace_printf(file_name):
  full_path = op.join(path_lib_src, file_name)
  if (not op.isfile(full_path)):
    full_path = op.join(op.dirname(path_lib_src), "fortran", file_name)
  result = ["#include <ccp4io_adaptbx/printf_wrappers.h>"]
  for line in open(full_path).read().splitlines():
    for key in ["printf", "fprintf"]:
      matches = list(re.finditer(
        pattern="[^A-Za-z0-9_]%s[^A-Za-z0-9_]" % key, string=line))
      if (len(matches) != 0):
        for m in reversed(matches):
          s,e = m.start(), m.end()
          line = line[:s] \
               + line[s:e].replace(key, "ccp4io_%s" % key) \
               + line[e:]
    result.append(line)
  return "\n".join(result)

env = env_base.Clone(
  SHLINKFLAGS=env_etc.shlinkflags)
env.Append(CCFLAGS=env_etc.ccp4io_defines)
env.Append(SHCCFLAGS=env_etc.ccp4io_defines)
env_etc.include_registry.append(
  env=env,
  paths=[
    "#",
    op.dirname(env_etc.ccp4io_include),
    env_etc.ccp4io_include,
    op.join(env_etc.ccp4io_dist)])
env.Append(LIBS=env_etc.libm)
# XXX 2012-06-16: is this actually necessary here, or just in code that links to
# ccp4io.lib?
if (os.name == "nt") :
  env.Prepend(LIBS=["Advapi32"])
if (   op.normcase(op.dirname(env_etc.ccp4io_dist))
    != op.normcase("ccp4io")):
  env.Repository(op.dirname(env_etc.ccp4io_dist))
source = []

c_files = []
c_files.extend(["%s/%s" % (ccp4_src, bn ) for bn in """\
library_err.c
library_file.c
library_utils.c
ccp4_array.c
ccp4_parser.c
ccp4_unitcell.c
cvecmat.c
cmtzlib.c
""".splitlines()])
open(op.join(build_ccp4io_adaptbx, "csymlib.c"), "w").write(
  open(op.join(path_lib_src, "csymlib.c")).read()
    .replace(
      "static int reported_syminfo = 0",
      "static int reported_syminfo = 1"))
source.append(op.join("#ccp4io_adaptbx", "csymlib.c"))

probe_file_name = op.join(path_lib_src, "cmaplib.h")
env_etc.ccp4io_has_cmaplib = op.isfile(probe_file_name)
if (env_etc.ccp4io_has_cmaplib):
  c_files.extend(["%s/%s" % (ccp4_src, bn ) for bn in """\
cmap_accessor.c
cmap_close.c
cmap_data.c
cmap_header.c
cmap_labels.c
cmap_open.c
cmap_skew.c
cmap_stats.c
cmap_symop.c
""".splitlines()])

c_files.extend(["%s/%s.cpp" % ( "mmdb", bn ) for bn in """\
hybrid_36
mmdb_atom
mmdb_bondmngr
mmdb_chain
mmdb_cifdefs
mmdb_coormngr
mmdb_cryst
mmdb_ficif
mmdb_io_file
mmdb_io_stream
mmdb_machine_
mmdb_manager
mmdb_mask
mmdb_math_align
mmdb_math_bfgsmin
mmdb_math_
mmdb_math_fft
mmdb_math_graph
mmdb_math_linalg
mmdb_math_rand
mmdb_mattype
mmdb_mmcif_
mmdb_model
mmdb_root
mmdb_rwbrook
mmdb_selmngr
mmdb_seqsuperpose
mmdb_symop
mmdb_tables
mmdb_title
mmdb_uddata
mmdb_utils
mmdb_xml_
""".splitlines()])
prefix = "#"+op.join(op.basename(env_etc.ccp4io_dist), "")
for file_name in c_files:
  source.append(op.join(prefix, file_name))

ssm_prefix = "#"+op.join(op.basename(env_etc.ccp4io_dist), "ssm")
ssm_sources = """\
ssm_csia.cpp
ssm_graph.cpp
ssm_vxedge.cpp
ssm_align.cpp
ssm_superpose.cpp
ssm_malign.cpp
""".splitlines()
source.extend( [ op.join( ssm_prefix, f ) for f in ssm_sources ] )

need_f_c = (
     libtbx.env.has_module("solve_resolve")
  or libtbx.env.find_in_repositories(relative_path="mosflm_fable"))
if (need_f_c or os.name != "nt"):
  source.append(op.join("#ccp4io_adaptbx", "fortran_call_stubs.c"))
  for file_name in """\
ccp4_diskio_f.c
ccp4_general.c
ccp4_general_f.c
ccp4_parser_f.c
ccp4_program.c
ccp4_unitcell_f.c
cmaplib_f.c
cmtzlib_f.c
csymlib_f.c
library_f.c
""".splitlines():
    open(op.join(build_ccp4io_adaptbx, file_name), "w").write(
      replace_printf(file_name=file_name))
    source.append(op.join("#ccp4io_adaptbx", file_name))
  source.append(op.join("#ccp4io_adaptbx", "printf_wrappers.c"))

# static library for solve_resolve
env.StaticLibrary(target='#lib/ccp4io', source=source)
env_etc.ccp4io_lib = "ccp4io"

if (    libtbx.env.has_module("boost")
    and not env_etc.no_boost_python):
  Import( "env_no_includes_boost_python_ext" )
  sources = [ "#ccp4io_adaptbx/ext.cpp" ]
  env_ext = env_no_includes_boost_python_ext.Clone()
  env_ext.Prepend( LIBS = "ccp4io" )
  env_etc.include_registry.append(
    env = env_ext,
    paths = [
      os.path.join( env_etc.ccp4io_dist),
      env_etc.boost_include,
      env_etc.python_include,
      ]
    )
  env_ext.SharedLibrary( target = "#lib/ccp4io_adaptbx_ext", source = sources )
