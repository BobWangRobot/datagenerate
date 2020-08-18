#ifndef CCP4IO_ADAPTBX_CCP4_PARSER_F_H
#define CCP4IO_ADAPTBX_CCP4_PARSER_F_H

#include "ccp4_sysdep.h"
#include "ccp4_parser.h"
#include "ccp4_general.h"

#ifdef  __cplusplus
namespace CCP4 {
extern "C" {
#endif

FORTRAN_SUBR(PARSER,parser,
             (fpstr key, fpstr line, int *ibeg, int *iend, int *ityp, float *fvalue,
              fpstr cvalue, int *idec, int *ntok, ftn_logical *lend,
              const ftn_logical *print, int key_len, int line_len, int cvalue_len),
             (fpstr key, fpstr line, int *ibeg, int *iend, int *ityp, float *fvalue,
              fpstr cvalue, int *idec, int *ntok, ftn_logical *lend,
              const ftn_logical *print),
             (fpstr key, int key_len, fpstr line, int line_len, int *ibeg, int *iend,
              int *ityp, float *fvalue, fpstr cvalue, int cvalue_len, int *idec,
              int *ntok, ftn_logical *lend, const ftn_logical *print));

FORTRAN_SUBR(PARSE,parse,
             (fpstr line, int *ibeg, int *iend, int *ityp, float *fvalue,
              fpstr cvalue, int *idec, int *n, int line_len, int cvalue_len),
             (fpstr line, int *ibeg, int *iend, int *ityp, float *fvalue,
              fpstr cvalue, int *idec, int *n),
             (fpstr line, int line_len, int *ibeg, int *iend,
              int *ityp, float *fvalue, fpstr cvalue, int cvalue_len, int *idec,
              int *n));

FORTRAN_SUBR(PARSDL,parsdl,
             (fpstr newdlm, int *nnewdl, int *nspecd, int newdlm_len),
             (fpstr newdlm, int *nnewdl, int *nspecd),
             (fpstr newdlm, int newdlm_len, int *nnewdl, int *nspecd));

#ifdef __cplusplus
}}
#endif

#endif // GUARD
