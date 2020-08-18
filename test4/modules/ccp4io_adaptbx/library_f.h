#ifndef CCP4IO_ADAPTBX_LIBRARY_F_H
#define CCP4IO_ADAPTBX_LIBRARY_F_H

#include "ccp4_utils.h"
#include "ccp4_errno.h"
#include "ccp4_fortran.h"

#ifdef  __cplusplus
extern "C" {
#endif

/* MOST FUNCTIONS OMITTED (just to save time) */

FORTRAN_SUBR ( QNAN, qnan,
    (union float_uint_uchar *realnum),
    (union float_uint_uchar *realnum),
    (union float_uint_uchar *realnum));

FORTRAN_FUN (int, QISNAN, qisnan,
             (union float_uint_uchar *realnum),
             (union float_uint_uchar *realnum),
             (union float_uint_uchar *realnum));

#ifdef __cplusplus
}
#endif

#endif /* GUARD */
