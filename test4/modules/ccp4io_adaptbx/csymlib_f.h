#ifndef CCP4IO_ADAPTBX_CSYMLIB_F_H
#define CCP4IO_ADAPTBX_CSYMLIB_F_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ccp4_fortran.h"
#include "ccp4_general.h"
#include "ccp4_parser.h"
#include "csymlib.h"
#include "cmtzlib.h"
#include "cvecmat.h"

#define MSPAC 4
#define MAXSYM 192
#define MAXSYMOPS 20
#define MAXLENSYMOPSTR 80

#ifdef  __cplusplus
namespace CSym {
extern "C" {
#endif

FORTRAN_SUBR ( INVSYM, invsym,
               (const float a[4][4], float ai[4][4]),
               (const float a[4][4], float ai[4][4]),
               (const float a[4][4], float ai[4][4]));

FORTRAN_SUBR ( SYMFR3, symfr3,
               (const fpstr icol, const int *i1, int *nsym, float rot[MAXSYM][4][4],
                     int *eflag, int icol_len),
               (const fpstr icol, const int *i1, int *nsym, float rot[MAXSYM][4][4],
                     int *eflag),
               (const fpstr icol, int icol_len, const int *i1, int *nsym,
                     float rot[MAXSYM][4][4], int *eflag));

FORTRAN_SUBR( SYMFR2, symfr2,
              (fpstr symchs, int *icol, int *nsym, float rot[MAXSYM][4][4], int symchs_len),
              (fpstr symchs, int *icol, int *nsym, float rot[MAXSYM][4][4]),
              (fpstr symchs, int symchs_len, int *icol, int *nsym, float rot[MAXSYM][4][4]));

FORTRAN_SUBR ( SYMTR3, symtr3,
               (const int *nsm, const float rsm[MAXSYM][4][4],
                     fpstr symchs, const int *iprint, int symchs_len),
               (const int *nsm, const float rsm[MAXSYM][4][4],
                     fpstr symchs, const int *iprint),
               (const int *nsm, const float rsm[MAXSYM][4][4],
                     fpstr symchs, int symchs_len, const int *iprint));

FORTRAN_SUBR ( SYMTR4, symtr4,
               (const int *nsm, const float rsm[MAXSYM][4][4],
                     fpstr symchs, int symchs_len),
               (const int *nsm, const float rsm[MAXSYM][4][4],
                     fpstr symchs),
               (const int *nsm, const float rsm[MAXSYM][4][4],
                     fpstr symchs, int symchs_len));

FORTRAN_SUBR ( PGMDF, pgmdf,
               (int *jlass, int*jcentr, int jscrew[3]),
               (int *jlass, int*jcentr, int jscrew[3]),
               (int *jlass, int*jcentr, int jscrew[3]));

FORTRAN_SUBR ( PGDEFN, pgdefn,
               (fpstr nampg, int *nsymp, const int *nsym, float rsmt[192][4][4],
                const ftn_logical *lprint, int nampg_len),
               (fpstr nampg, int *nsymp, const int *nsym, float rsmt[192][4][4],
                const ftn_logical *lprint),
               (fpstr nampg, int nampg_len, int *nsymp, const int *nsym,
                float rsmt[192][4][4], const ftn_logical *lprint));

FORTRAN_SUBR ( PGNLAU, pgnlau,
               (const fpstr nampg, int *nlaue, fpstr launam,
                int nampg_len, int launam_len),
               (const fpstr nampg, int *nlaue, fpstr launam),
               (const fpstr nampg, int nampg_len, int *nlaue,
                fpstr launam, int launam_len));

FORTRAN_SUBR ( CCP4SPG_F_GET_LAUE, ccp4spg_f_get_laue,
               (const int *sindx, int *nlaue, fpstr launam, int launam_len),
               (const int *sindx, int *nlaue, fpstr launam),
               (const int *sindx, int *nlaue, fpstr launam, int launam_len));

FORTRAN_SUBR ( HKLRANGE, hklrange,
               (int *ihrng0, int *ihrng1, int *ikrng0, int *ikrng1, int *ilrng0, int *ilrng1),
               (int *ihrng0, int *ihrng1, int *ikrng0, int *ikrng1, int *ilrng0, int *ilrng1),
               (int *ihrng0, int *ihrng1, int *ikrng0, int *ikrng1, int *ilrng0, int *ilrng1));

FORTRAN_SUBR ( PATSGP, patsgp,
               (const fpstr spgnam, const fpstr pgname, fpstr patnam, int *lpatsg,
                int spgnam_len, int pgname_len, int patnam_len),
               (const fpstr spgnam, const fpstr pgname, fpstr patnam, int *lpatsg),
               (const fpstr spgnam, int spgnam_len, const fpstr pgname,
                int pgname_len, fpstr patnam, int patnam_len, int *lpatsg));

FORTRAN_SUBR ( ASUSET, asuset,
               (fpstr spgnam, int *numsgp, fpstr pgname,
                int *msym, float rrsym[192][4][4], int *msymp,
                int *mlaue, ftn_logical *lprint, int spgnam_len, int pgname_len),
               (fpstr spgnam, int *numsgp, fpstr pgname,
                int *msym, float rrsym[192][4][4], int *msymp,
                int *mlaue, ftn_logical *lprint),
               (fpstr spgnam, int spgnam_len, int *numsgp,
                fpstr pgname,int pgname_len,
                int *msym, float rrsym[192][4][4], int *msymp,
                int *mlaue, ftn_logical *lprint));

FORTRAN_SUBR ( ASUSYM, asusym,
               (float rassym[384][4][4], float rinsym[384][4][4], int *nisym),
               (float rassym[384][4][4], float rinsym[384][4][4], int *nisym),
               (float rassym[384][4][4], float rinsym[384][4][4], int *nisym));

FORTRAN_SUBR ( ASUPUT, asuput,
               (const int ihkl[3], int jhkl[3], int *isym),
               (const int ihkl[3], int jhkl[3], int *isym),
               (const int ihkl[3], int jhkl[3], int *isym));

FORTRAN_SUBR ( ASUGET, asuget,
               (const int ihkl[3], int jhkl[3], const int *isym),
               (const int ihkl[3], int jhkl[3], const int *isym),
               (const int ihkl[3], int jhkl[3], const int *isym));

FORTRAN_SUBR ( ASUPHP, asuphp,
               (const int jhkl[3], const int *lsym, const int *isign,
                const float *phasin, float *phasout),
               (const int jhkl[3], const int *lsym, const int *isign,
                const float *phasin, float *phasout),
               (const int jhkl[3], const int *lsym, const int *isign,
                const float *phasin, float *phasout));

FORTRAN_SUBR ( CCP4SPG_F_LOAD_BY_NAME, ccp4spg_f_load_by_name,
               (const int *sindx, fpstr namspg, int namspg_len),
               (const int *sindx, fpstr namspg),
               (const int *sindx, fpstr namspg, int namspg_len));

FORTRAN_SUBR ( CCP4SPG_F_LOAD_BY_OPS, ccp4spg_f_load_by_ops,
               (const int *sindx, int *msym, float rrsym[192][4][4]),
               (const int *sindx, int *msym, float rrsym[192][4][4]),
               (const int *sindx, int *msym, float rrsym[192][4][4]));

FORTRAN_FUN (int, CCP4SPG_F_EQUAL_OPS_ORDER, ccp4spg_f_equal_ops_order,
               (int *msym1, float rrsym1[192][4][4],int *msym2, float rrsym2[192][4][4]),
               (int *msym1, float rrsym1[192][4][4],int *msym2, float rrsym2[192][4][4]),
               (int *msym1, float rrsym1[192][4][4],int *msym2, float rrsym2[192][4][4]));

FORTRAN_SUBR ( CCP4SPG_F_ASUPUT, ccp4spg_f_asuput,
               (const int *sindx, const int ihkl[3], int jhkl[3], int *isym),
               (const int *sindx, const int ihkl[3], int jhkl[3], int *isym),
               (const int *sindx, const int ihkl[3], int jhkl[3], int *isym));

FORTRAN_FUN (int, INASU, inasu,
               (const int ihkl[3], const int *nlaue),
               (const int ihkl[3], const int *nlaue),
               (const int ihkl[3], const int *nlaue));

FORTRAN_FUN (int, CCP4SPG_F_INASU, ccp4spg_f_inasu,
               (const int *sindx, const int ihkl[3]),
               (const int *sindx, const int ihkl[3]),
               (const int *sindx, const int ihkl[3]));

FORTRAN_SUBR ( PRTRSM, prtrsm,
               (const fpstr pgname, const int *nsymp,
                const float rsymiv[192][4][4], int pgname_len),
               (const fpstr pgname, const int *nsymp,
                const float rsymiv[192][4][4]),
               (const fpstr pgname, int pgname_len, const int *nsymp,
                const float rsymiv[192][4][4]));

FORTRAN_SUBR ( MSYMLB3, msymlb3,
               (const int *ist, int *lspgrp, fpstr namspg_cif,
                fpstr namspg_cifs, fpstr nampg, int *nsymp, int *nsym,
                float rlsymmmatrx[192][4][4], int namspg_cif_len,
                int namspg_cifs_len, int nampg_len),
               (const int *ist, int *lspgrp, fpstr namspg_cif,
                fpstr namspg_cifs, fpstr nampg, int *nsymp, int *nsym,
                float rlsymmmatrx[192][4][4]),
               (const int *ist, int *lspgrp, fpstr namspg_cif, int namspg_cif_len,
                fpstr namspg_cifs, int namspg_cifs_len, fpstr nampg, int nampg_len,
                int *nsymp, int *nsym, float rlsymmmatrx[192][4][4]));

FORTRAN_SUBR ( MSYMLB, msymlb,
               (const int *ist, int *lspgrp, fpstr namspg_cif,
                fpstr nampg, int *nsymp, int *nsym,
                float rlsymmmatrx[192][4][4], int namspg_cif_len,
                int nampg_len),
               (const int *ist, int *lspgrp, fpstr namspg_cif,
                fpstr nampg, int *nsymp, int *nsym,
                float rlsymmmatrx[192][4][4]),
               (const int *ist, int *lspgrp, fpstr namspg_cif, int namspg_cif_len,
                fpstr nampg, int nampg_len,
                int *nsymp, int *nsym, float rlsymmmatrx[192][4][4]));

FORTRAN_SUBR ( MSYMLB2, msymlb2,
               (const int *ist, int *lspgrp, fpstr namspg_cif,
                fpstr nampg, int *nsymp, int *nsym,
                float rlsymmmatrx[192][4][4], int namspg_cif_len,
                int nampg_len),
               (const int *ist, int *lspgrp, fpstr namspg_cif,
                fpstr nampg, int *nsymp, int *nsym,
                float rlsymmmatrx[192][4][4]),
               (const int *ist, int *lspgrp, fpstr namspg_cif, int namspg_cif_len,
                fpstr nampg, int nampg_len,
                int *nsymp, int *nsym, float rlsymmmatrx[192][4][4]));

FORTRAN_SUBR ( MSYGET, msyget,
               (const int *ist, int *lspgrp, int *nsym,
                float rlsymmmatrx[192][4][4]),
               (const int *ist, int *lspgrp, int *nsym,
                float rlsymmmatrx[192][4][4]),
               (const int *ist, int *lspgrp, int *nsym,
                float rlsymmmatrx[192][4][4]));

FORTRAN_SUBR ( EPSLN, epsln,
               (const int *nsm, const int *nsmp, const float rsm[192][4][4],
                const int *iprint),
               (const int *nsm, const int *nsmp, const float rsm[192][4][4],
                const int *iprint),
               (const int *nsm, const int *nsmp, const float rsm[192][4][4],
                const int *iprint));

FORTRAN_SUBR ( EPSLON, epslon,
               (const int ih[3], float *epsi, int *isysab),
               (const int ih[3], float *epsi, int *isysab),
               (const int ih[3], float *epsi, int *isysab));

FORTRAN_SUBR ( CCP4SPG_F_EPSLON, ccp4spg_f_epslon,
               (const int *sindx, const int ih[3], float *epsi, int *isysab),
               (const int *sindx, const int ih[3], float *epsi, int *isysab),
               (const int *sindx, const int ih[3], float *epsi, int *isysab));

FORTRAN_SUBR ( SYSAB, sysab,
               (const int in[3], int *isysab),
               (const int in[3], int *isysab),
               (const int in[3], int *isysab));

FORTRAN_SUBR ( CCP4SPG_F_IS_SYSABS, ccp4spg_f_is_sysabs,
               (const int *sindx, const int in[3], int *isysab),
               (const int *sindx, const int in[3], int *isysab),
               (const int *sindx, const int in[3], int *isysab));

FORTRAN_SUBR ( CENTRIC, centric,
               (const int *nsm, const float rsm[192][4][4],
                const int *iprint),
               (const int *nsm, const float rsm[192][4][4],
                const int *iprint),
               (const int *nsm, const float rsm[192][4][4],
                const int *iprint));

FORTRAN_SUBR ( CENTR, centr,
               (const int ih[3], int *ic),
               (const int ih[3], int *ic),
               (const int ih[3], int *ic));

FORTRAN_SUBR ( CCP4SPG_F_IS_CENTRIC, ccp4spg_f_is_centric,
               (const int *sindx, const int ih[3], int *ic),
               (const int *sindx, const int ih[3], int *ic),
               (const int *sindx, const int ih[3], int *ic));

FORTRAN_SUBR ( CENTPHASE, centphase,
               (const int ih[3], float *cenphs),
               (const int ih[3], float *cenphs),
               (const int ih[3], float *cenphs));

FORTRAN_SUBR ( CCP4SPG_F_CENTPHASE, ccp4spg_f_centphase,
               (const int *sindx, const int ih[3], float *cenphs),
               (const int *sindx, const int ih[3], float *cenphs),
               (const int *sindx, const int ih[3], float *cenphs));

FORTRAN_SUBR ( SETLIM, setlim,
               (const int *lspgrp, float xyzlim[3][2]),
               (const int *lspgrp, float xyzlim[3][2]),
               (const int *lspgrp, float xyzlim[3][2]));

FORTRAN_SUBR ( SETLIM_ZERO, setlim_zero,
               (const int *lspgrp, float xyzlim[3][2]),
               (const int *lspgrp, float xyzlim[3][2]),
               (const int *lspgrp, float xyzlim[3][2]));

FORTRAN_SUBR ( SETGRD, setgrd,
               (const int *nlaue, const float *sample, const int *nxmin,
                const int *nymin, const int *nzmin, int *nx, int *ny, int *nz),
               (const int *nlaue, const float *sample, const int *nxmin,
                const int *nymin, const int *nzmin, int *nx, int *ny, int *nz),
               (const int *nlaue, const float *sample, const int *nxmin,
                const int *nymin, const int *nzmin, int *nx, int *ny, int *nz));

FORTRAN_SUBR ( FNDSMP, fndsmp,
               (const int *minsmp, const int *nmul, const float *sample, int *nsampl),
               (const int *minsmp, const int *nmul, const float *sample, int *nsampl),
               (const int *minsmp, const int *nmul, const float *sample, int *nsampl));

FORTRAN_SUBR ( CALC_ORIG_PS, calc_orig_ps,
               (fpstr namspg_cif, int *nsym, float rsym[192][4][4], int *norig,
                float orig[96][3], ftn_logical *lpaxisx, ftn_logical *lpaxisy,
                ftn_logical *lpaxisz, int namspg_cif_len),
               (fpstr namspg_cif, int *nsym, float rsym[192][4][4], int *norig,
                float orig[96][3], ftn_logical *lpaxisx, ftn_logical *lpaxisy,
                ftn_logical *lpaxisz, int namspg_cif_len),
               (fpstr namspg_cif, int namspg_cif_len, int *nsym, float rsym[192][4][4],
                int *norig, float orig[96][3], ftn_logical *lpaxisx,
                ftn_logical *lpaxisy, ftn_logical *lpaxisz));

FORTRAN_SUBR ( SETRSL, setrsl,
               (const float *a, const float *b, const float *c,
                const float *alpha, const float *beta, const float *gamma),
               (const float *a, const float *b, const float *c,
                const float *alpha, const float *beta, const float *gamma),
               (const float *a, const float *b, const float *c,
                const float *alpha, const float *beta, const float *gamma));

FORTRAN_SUBR (STHLSQ1, sthlsq1,
             (float *reso, const int *ih, const int *ik, const int *il),
             (float *reso, const int *ih, const int *ik, const int *il),
             (float *reso, const int *ih, const int *ik, const int *il));

FORTRAN_SUBR (STS3R41, sts3r41,
             (float *reso, const int *ih, const int *ik, const int *il),
             (float *reso, const int *ih, const int *ik, const int *il),
             (float *reso, const int *ih, const int *ik, const int *il));

FORTRAN_SUBR ( HANDCHANGE, handchange,
               (const int *lspgrp, float *cx, float *cy, float *cz),
               (const int *lspgrp, float *cx, float *cy, float *cz),
               (const int *lspgrp, float *cx, float *cy, float *cz));

#ifdef __cplusplus
} }
#endif

#endif // GUARD
