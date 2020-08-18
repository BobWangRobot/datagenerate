#ifndef CCP4IO_ADAPTBX_CCP4_OTHER_FEM_HPP
#define CCP4IO_ADAPTBX_CCP4_OTHER_FEM_HPP

#include <fem.hpp> // Fortran EMulation library of fable module

namespace ccp4_other_fem {

using namespace fem::major_types;

inline
bool
vaxvms()
{
  return false;
}

inline
void
ugerr(
  int const& /*status*/,
  str_ref errstr)
{
  throw TBXX_NOT_IMPLEMENTED();
}

inline
void
ccp4h_summary_beg()
{
}

inline
void
ccp4h_summary_end()
{
}

inline
void
ubytes(
  int& inum,
  str_ref string)
{
  inum = 4;
  string = "BYTES";
}

inline
void
ugtenv(
  str_cref namenv,
  str_ref valenv)
{
  fem::getenv(namenv, valenv);
}

inline
void
imsiz(
  str_cref /* odfile */,
  int& nrec,
  int& iylen)
{
  throw TBXX_NOT_IMPLEMENTED();
}

inline
void
readpack_word(
  arr_ref<fem::integer_star_2> datapack,
  str_cref /* filn */)
{
  throw TBXX_NOT_IMPLEMENTED();
}

} // namespace ccp4_other_fem

#endif // GUARD
