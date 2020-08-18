#ifndef CCP4IO_ADAPTBX_CCP4_PARSER_FEM_HPP
#define CCP4IO_ADAPTBX_CCP4_PARSER_FEM_HPP

#include <fem.hpp> // Fortran EMulation library of fable module

extern "C" {

typedef unsigned int ccp4_ftn_logical;

void
parser_(
  char* key,
  char* line,
  int* ibeg,
  int* iend,
  int* ityp,
  float* fvalue,
  /*str_arr_ref<>*/ char* cvalue,
  int* idec,
  int* ntok,
  ccp4_ftn_logical* lend,
  ccp4_ftn_logical const* print,
  int key_len,
  int line_len,
  int cvalue_len);

} // extern "C"

namespace ccp4_parser_fem {

using namespace fem::major_types;

inline
void
lerror(
  int const& errflg,
  int const& ifail,
  str_cref errmsg)
{
  throw TBXX_NOT_IMPLEMENTED();
}

inline
void
parser(
  str_ref key,
  str_ref line,
  arr_ref<int> ibeg,
  arr_ref<int> iend,
  arr_ref<int> ityp,
  arr_ref<float> fvalue,
  str_arr_ref<> cvalue,
  arr_ref<int> idec,
  int& ntok,
  bool& lend,
  bool const& print)
{
  ccp4_ftn_logical lend_ = lend;
  ccp4_ftn_logical print_ = print;
  parser_(
    key.elems(),
    line.elems(),
    ibeg.begin(),
    iend.begin(),
    ityp.begin(),
    fvalue.begin(),
    cvalue.begin(),
    idec.begin(),
    &ntok,
    &lend_,
    &print_,
    key.len(),
    line.len(),
    cvalue.len());
  lend = static_cast<bool>(lend_);
}

inline
void
putlin(
  str_cref strout,
  str_cref outwin)
{
  throw TBXX_NOT_IMPLEMENTED();
}

} // namespace ccp4_parser_fem

#endif // GUARD
