#ifndef CCP4IO_ADAPTBX_CCP4_CCPLIB_FEM_HPP
#define CCP4IO_ADAPTBX_CCP4_CCPLIB_FEM_HPP

#include <fem.hpp> // Fortran EMulation library of fable module

extern "C" {

void
ccp4_version_(
  char* version,
  int version_len);

void
ccpfyp_();

void
ccpdat_(
  char* caldat,
  int caldat_len);

}

namespace ccp4_ccplib_fem {

using namespace fem::major_types;

inline
void
ccpupc(
  str_ref string)
{
  int n = string.len();
  char* s = string.elems();
  for(int i=0;i<n;i++) {
    s[i] = fem::utils::to_upper(s[i]);
  }
}

inline
int
lenstr(
  str_cref string)
{
  return fem::len_trim(string);
}

inline
void
ccperr(
  int const& istat,
  str_cref errstr)
{
  throw TBXX_NOT_IMPLEMENTED();
}

inline
void
ccpdat(
  str_ref caldat)
{
  ccpdat_(caldat.elems(), caldat.len());
}

inline
void
ccplwc(
  str_ref string)
{
  throw TBXX_NOT_IMPLEMENTED();
}

inline
bool
litend(
  int const& /* idum */)
{
  bool return_value = fem::bool0;
  throw TBXX_NOT_IMPLEMENTED();
}

inline
bool
litend()
{
  // true if the system byte-order is little endian
  short int word = 0x0001;
  char* byte = (char*) &word;
  return byte[0]==0x01;
}

inline
bool
ccponl(
  int const& /* idum */)
{
  return false;
}

inline
int
ccpe2i(
  str_cref name,
  int const& defval)
{
  std::string k = fem::utils::strip_leading_and_trailing_blank_padding(name);
  char* v = std::getenv(k.c_str());
  if (v == 0) return defval;
  int result;
  int n = std::sscanf(v, "%d", result);
  if (n != 1) {
    std::ostringstream o;
    o << "Environment variable \""
      << std::string(k)
      << "\" is not an integer: \""
      << v << "\"";
    throw std::runtime_error(o.str());
  }
  return result;
}

inline
void
ccp4_version(
  str_ref version)
{
  ccp4_version_(version.elems(),version.len());
}

inline
void
ccpfyp()
{
  ccpfyp_();
}

inline
int
lunsto(
  int const& /*idum*/)
{
  return 6;
}

inline
bool
hkleq(
  arr_cref<int> ih,
  arr_cref<int> kh)
{
  return (ih(1) == kh(1)) && (ih(2) == kh(2)) && (ih(3) == kh(3));
}

} // namespace ccp4_ccplib_fem

#endif // GUARD
