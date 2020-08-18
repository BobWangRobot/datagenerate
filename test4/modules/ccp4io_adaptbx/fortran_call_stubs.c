#include "ccp4_fortran.h"
#include <stdio.h>

static void
bail(const char* function_name)
{
  fprintf(stderr, "FATAL ERROR: %s() not implemented.", function_name);
  exit(1);
}

FORTRAN_FUN ( int, IARGC, iargc,
              (),
              (),
              ())
{
  /*bail("FORTRAN iargc");*/
  return 0;
}

FORTRAN_SUBR ( GETARG, getarg,
               (int *i, char *arg, int arg_len),
               (int *i, char *arg, int arg_len),
               (int *i, char *arg, int arg_len))
{ /*bail("FORTRAN getarg");*/
  *arg = '\0';
}

FORTRAN_SUBR ( CCP4H_INIT, ccp4h_init,
               (),
               (),
               ())
{ /*bail("ccp4h_init");*/
  /*access environment variables and construct output*/}

FORTRAN_SUBR ( CCP4H_INIT_LIB, ccp4h_init_lib,
               (int* ihtml, int* isumm),
               (int* ihtml, int* isumm),
               (int* ihtml, int* isumm))
{
  // ignoring calls from mtzini() in cmtzlib_f.c
}

FORTRAN_SUBR ( CCP4H_PRE_BEG, ccp4h_pre_beg,
               (),
               (),
               ())
{ bail("ccp4h_pre_beg"); }

FORTRAN_SUBR ( CCP4H_SUMMARY_BEG, ccp4h_summary_beg,
               (),
               (),
               ())
{ bail("ccp4h_summary_beg"); }

FORTRAN_SUBR ( CCP4H_SUMMARY_END, ccp4h_summary_end,
               (),
               (),
               ())
{ bail("ccp4h_summary_end"); }
