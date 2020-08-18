#ifndef CCP4IO_ADAPTBX_CCP4_LIBRARY_FEM_HPP
#define CCP4IO_ADAPTBX_CCP4_LIBRARY_FEM_HPP

#include <fem.hpp> // Fortran EMulation library of fable module

extern "C" {

void
qnan_(
  float* value);

} // extern "C"

namespace ccp4_library_fem {

using namespace fem::major_types;

inline
void
ccpupc(
  str_ref string)
{
  throw TBXX_NOT_IMPLEMENTED();
}

inline
int
lenstr(
  str_cref string)
{
  throw TBXX_NOT_IMPLEMENTED();
}

inline
void
qnan(
  float& value)
{
  qnan_(&value);
}

} // namespace ccp4_library_fem

#endif // GUARD
