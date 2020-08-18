// NOT USED -- INCOMPLETE -- DOES NOT COMPILE
// here for future reference

#ifndef CCP4IO_ADAPTBX_DEV_FFTLIB_HPP
#define CCP4IO_ADAPTBX_DEV_FFTLIB_HPP

namespace ccp4io_dev {

void
cmplft(
  FArray1A_float x,
  FArray1A_float y,
  int const n,
  int5a d
);

void
hermft(
  FArray1A_float x,
  FArray1A_float y,
  int const n,
  int5a dim
);

void
realft(
  FArray1A_float even,
  FArray1A_float odd,
  int const n,
  int5a dim
);

void
rsymft(
  FArray1A_float x,
  int const n,
  int5a dim
);

void
sdiad(
  FArray1A_float x,
  FArray1A_float y,
  int const n,
  int5a dim
);

void
inv21(
  FArray1A_float x,
  FArray1A_float y,
  int const n,
  int5a d
);

} // namespace ccp4io_dev

#endif // GUARD
