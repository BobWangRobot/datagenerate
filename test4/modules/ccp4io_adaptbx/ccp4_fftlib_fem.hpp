#ifndef CCP4IO_ADAPTBX_CCP4_FFTLIB_FEM_HPP
#define CCP4IO_ADAPTBX_CCP4_FFTLIB_FEM_HPP

#include <fem.hpp> // Fortran EMulation library of fable module

namespace ccp4_fftlib_fem {

using namespace fem::major_types;

void
srfp(
  int const& pts,
  int const& pmax,
  int const& twogrp,
  arr_ref<int> factor,
  arr_ref<int> sym,
  int& psym,
  arr_ref<int> unsym,
  bool& error);

void
r2cftk(
  int const& n,
  int const& m,
  arr_ref<float> x0,
  arr_ref<float> y0,
  arr_ref<float> x1,
  arr_ref<float> y1,
  arr_cref<int> dim);

void
r3cftk(
  int const& n,
  int const& m,
  arr_ref<float> x0,
  arr_ref<float> y0,
  arr_ref<float> x1,
  arr_ref<float> y1,
  arr_ref<float> x2,
  arr_ref<float> y2,
  arr_cref<int> dim);

void
r4cftk(
  int const& n,
  int const& m,
  arr_ref<float> x0,
  arr_ref<float> y0,
  arr_ref<float> x1,
  arr_ref<float> y1,
  arr_ref<float> x2,
  arr_ref<float> y2,
  arr_ref<float> x3,
  arr_ref<float> y3,
  arr_cref<int> dim);

void
r5cftk(
  int const& n,
  int const& m,
  arr_ref<float> x0,
  arr_ref<float> y0,
  arr_ref<float> x1,
  arr_ref<float> y1,
  arr_ref<float> x2,
  arr_ref<float> y2,
  arr_ref<float> x3,
  arr_ref<float> y3,
  arr_ref<float> x4,
  arr_ref<float> y4,
  arr_cref<int> dim);

void
r8cftk(
  int const& n,
  int const& m,
  arr_ref<float> x0,
  arr_ref<float> y0,
  arr_ref<float> x1,
  arr_ref<float> y1,
  arr_ref<float> x2,
  arr_ref<float> y2,
  arr_ref<float> x3,
  arr_ref<float> y3,
  arr_ref<float> x4,
  arr_ref<float> y4,
  arr_ref<float> x5,
  arr_ref<float> y5,
  arr_ref<float> x6,
  arr_ref<float> y6,
  arr_ref<float> x7,
  arr_ref<float> y7,
  arr_cref<int> dim);

void
rpcftk(
  int const& n,
  int const& m,
  int const& p,
  int const& r,
  arr_ref<float, 2> x,
  arr_ref<float, 2> y,
  arr_cref<int> dim);

void
mdftkd(
  int const& n,
  arr_cref<int> factor,
  arr_cref<int> dim,
  arr_ref<float> x,
  arr_ref<float> y);

void
diprp(
  int const& pts,
  arr_cref<int> sym,
  int const& psym,
  arr_cref<int> unsym,
  arr_cref<int> dim,
  arr_ref<float> x,
  arr_ref<float> y);

void
cmplft(
  arr_ref<float> x,
  arr_ref<float> y,
  int const& n,
  arr_cref<int> d);

void
hermft(
  arr_ref<float> x,
  arr_ref<float> y,
  int const& n,
  arr_cref<int> dim);

void
realft(
  arr_ref<float> even,
  arr_ref<float> odd,
  int const& n,
  arr_cref<int> dim);

} // namespace ccp4_fftlib_fem

#endif // GUARD
