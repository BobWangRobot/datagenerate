#ifndef CCP4IO_ADAPTBX_CCP4_DISKIO_FEM_HPP
#define CCP4IO_ADAPTBX_CCP4_DISKIO_FEM_HPP

#include <fem.hpp> // Fortran EMulation library of fable module

namespace ccp4_diskio_fem {

extern "C" {

void qopen_(
  int const* iunit,
  char const* lognam,
  char const* atbuta,
  int lognam_len,
  int atbuta_len);

void qmode_(
  int const* iunit,
  int const* mode,
  int* size);

void qqinq_(
  int const* iunit,
  char const* lfn,
  char* filnam,
  int const* length,
  int lfn_len,
  int filnam_len);

void qread_(
  int const* iunit,
  unsigned long* buffer,
  int const* nitems,
  int* result);

void qclose_(
  int const* iunit);

void qseek_(
  int const* iunit,
  int const* irec,
  int const* iel,
  int const* lrecl);

void
qwrite_(
  int const* iunit,
  const unsigned long* buffer,
  int const* nitems);

bool
qisnan_(
  float const* value);

void
qback_(
  int const* iunit,
  int const* lrecl);

}

using namespace fem::major_types;

inline
void
qopen(
  int const& iunit,
  str_cref lognam,
  str_cref atbuta)
{
  qopen_(&iunit, lognam.elems(), atbuta.elems(), lognam.len(), atbuta.len());
}

inline
void
qwrite(
  int const& iunit,
  arr_cref<int> buffer,
  int const& nitems)
{
  qwrite_(&iunit, (const unsigned long*)(const void*)buffer.begin(),&nitems);
}

inline
void
qwrite(
  int const& iunit,
  fem::integer_star_2 const& buffer,
  int const& nitems)
{
  qwrite_(&iunit, (const unsigned long*)(const void*)(&buffer),&nitems);
}

inline
void
qwrite(
  int const& iunit,
  arr_ref<fem::integer_star_2> buffer,
  int const& nitems)
{
  qwrite_(&iunit, (const unsigned long*)(const void*)buffer.begin(),&nitems);
}

inline
void
qwrite(
  int const& iunit,
  arr_cref<fem::integer_star_8> buffer,
  int const& nitems)
{
  throw TBXX_NOT_IMPLEMENTED();
}

inline
bool
qisnan(
  float const& value)
{
  return qisnan_(&value);
}

inline
void
qclose(
  int const& iunit)
{
  qclose_(&iunit);
}

inline
void
qseek(
  int const& iunit,
  int const& irec,
  int const& iel,
  int const& lrecl)
{
  qseek_(&iunit, &irec, &iel, &lrecl);
}

inline
void
qmode(
  int const& iunit,
  int const& mode,
  int& size)
{
  qmode_(&iunit, &mode, &size);
}

inline
void
qread(
  int const& /* iunit */,
  arr_ref<float> buffer,
  int const& /* nitems */,
  int& result)
{
  throw TBXX_NOT_IMPLEMENTED();
}

inline
void
qread(
  int const& iunit,
  arr_ref<fem::integer_star_1> buffer,
  int const& nitems,
  int& result)
{
  qread_(&iunit, (unsigned long*)(void*)buffer.begin(),
         &nitems, &result);
}

inline
void
qread(
  int const& iunit,
  arr_ref<fem::integer_star_2> buffer,
  int const& nitems,
  int& result)
{
  qread_(&iunit, (unsigned long*)(void*)buffer.begin(),
         &nitems, &result);
}

inline
void
qread(
  int const& iunit,
  arr_ref<int> buffer,
  int const& nitems,
  int& result)
{
  qread_(&iunit, (unsigned long*)(void*)buffer.begin(),
         &nitems, &result);
}

inline
void
qqinq(
  int const& iunit,
  str_cref lfn,
  str_ref filnam,
  int const& length)
{
  qqinq_(&iunit, lfn.elems(), filnam.elems(), &length, lfn.len(), filnam.len());
}

inline
void
qback(
  int const& iunit,
  int const& lrecl)
{
  qback_(&iunit, &lrecl);
}

inline
void
qreadi(
  int const& iunit,
  arr_ref<int> buffer,
  int const& nitems,
  int& result)
{
  throw TBXX_NOT_IMPLEMENTED();
}

inline
void
qreadi2(
  int const& /* iunit */,
  arr_ref<fem::integer_star_2> buffer,
  int const& /* nitems */,
  int& result)
{
  throw TBXX_NOT_IMPLEMENTED();
}

} // namespace ccp4_diskio_fem

#endif // GUARD
