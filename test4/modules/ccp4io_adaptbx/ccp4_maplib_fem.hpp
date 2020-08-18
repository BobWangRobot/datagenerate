#ifndef CCP4IO_ADAPTBX_CCP4_MAPLIB_FEM_HPP
#define CCP4IO_ADAPTBX_CCP4_MAPLIB_FEM_HPP

#include <fem.hpp> // Fortran EMulation library of fable module

extern "C" {

void
mwrhdl_(
  int const* iunit,
  char const* mapnam,
  char const* title,
  int const* nsec,
  int const* iuvw,
  int const* mxyz,
  int const* nw1,
  int const* nu1,
  int const* nu2,
  int const* nv1,
  int const* nv2,
  float const* cell,
  int const* lspgrp,
  int* lmode,
  int mapnam_len,
  int title_len);

void
mwrsec_(
  int const* iunit,
  /*arr_cref<float, 2>*/ float const* x,
  int const* mu,
  int const* mv,
  int const* iu1,
  int const* iu2,
  int const* iv1,
  int const* iv2);

void
mwclose_(
  int const* iunit);

void
msywrt_(
  int const* iunit,
  int const* nsym,
  /*arr_cref<float, 3>*/ float const* rot);

} // extern "C"

namespace ccp4_maplib_fem {

using namespace fem::major_types;

inline
void
mwrhdl(
  int const& iunit,
  str_cref mapnam,
  str_cref title,
  int const& nsec,
  arr_cref<int> iuvw,
  arr_cref<int> mxyz,
  int const& nw1,
  int const& nu1,
  int const& nu2,
  int const& nv1,
  int const& nv2,
  arr_cref<float> cell,
  int const& lspgrp,
  int& lmode)
{
  mwrhdl_(
    &iunit,
    mapnam.elems(),
    title.elems(),
    &nsec,
    iuvw.begin(),
    mxyz.begin(),
    &nw1,
    &nu1,
    &nu2,
    &nv1,
    &nv2,
    cell.begin(),
    &lspgrp,
    &lmode,
    mapnam.len(),
    title.len());
}

inline
void
mwrsec(
  int const& iunit,
  arr_cref<float, 2> x,
  int const& mu,
  int const& mv,
  int const& iu1,
  int const& iu2,
  int const& iv1,
  int const& iv2)
{
  mwrsec_(&iunit, x.begin(), &mu, &mv, &iu1, &iu2, &iv1, &iv2);
}

inline
void
mwclose(
  int const& iunit)
{
  mwclose_(&iunit);
}

inline
void
msywrt(
  int const& iunit,
  int const& nsym,
  arr_cref<float, 3> rot)
{
  msywrt_(&iunit, &nsym, rot.begin());
}

} // namespace ccp4_maplib_fem

#endif // GUARD
