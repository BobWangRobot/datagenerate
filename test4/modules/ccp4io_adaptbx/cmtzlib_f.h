#ifndef CCP4IO_ADAPTBX_CMTZLIB_F_H
#define CCP4IO_ADAPTBX_CMTZLIB_F_H

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ccp4_fortran.h"
#include "ccp4_utils.h"
#include "cmtzlib.h"
#include "csymlib.h"
#include "ccp4_program.h"
#include "ccp4_general.h"

#define MAXSYM 192

#ifdef  __cplusplus
namespace CMtz {
extern "C" {
#endif

typedef char char_3[3];
typedef char char_31[31];
typedef char char_mtzrl[MTZRECORDLENGTH];
typedef char char_2_31[2][31];

FORTRAN_SUBR ( MTZINI, mtzini,
               ( ),
               ( ),
               ( ));

FORTRAN_SUBR ( LROPEN, lropen,
               (int *mindx, fpstr filename, int *iprint, int *ifail, int filename_len),
               (int *mindx, fpstr filename, int *iprint, int *ifail),
               (int *mindx, fpstr filename, int filename_len, int *iprint, int *ifail));

FORTRAN_SUBR ( LRTITL, lrtitl,
               (int *mindx, fpstr ftitle, int *len, int ftitle_len),
               (int *mindx, fpstr ftitle, int *len),
               (int *mindx, fpstr ftitle, int ftitle_len, int *len));

FORTRAN_SUBR ( LRHIST, lrhist,
               (int *mindx, fpstr hstrng, int *nlines, int hstrng_len),
               (int *mindx, fpstr hstrng, int *nlines),
               (int *mindx, fpstr hstrng, int hstrng_len, int *nlines));

FORTRAN_SUBR ( LRINFO, lrinfo,
               (int *mindx, fpstr versnx, int *ncolx, int *nreflx, float *ranges, int versnx_len),
               (int *mindx, fpstr versnx, int *ncolx, int *nreflx, float *ranges),
               (int *mindx, fpstr versnx, int versnx_len, int *ncolx, int *nreflx, float *ranges));

FORTRAN_SUBR ( LRNCOL, lrncol,
               (int *mindx, int *ncolx),
               (int *mindx, int *ncolx),
               (int *mindx, int *ncolx));

FORTRAN_SUBR ( LRNREF, lrnref,
               (int *mindx, int *nreflx),
               (int *mindx, int *nreflx),
               (int *mindx, int *nreflx));

FORTRAN_SUBR ( LRSORT, lrsort,
               (int *mindx, int sortx[5]),
               (int *mindx, int sortx[5]),
               (int *mindx, int sortx[5]));

FORTRAN_SUBR ( LRBATS, lrbats,
               (int *mindx, int *nbatx, int batchx[]),
               (int *mindx, int *nbatx, int batchx[]),
               (int *mindx, int *nbatx, int batchx[]));

FORTRAN_SUBR ( LRCLAB, lrclab,
               (int *mindx, fpstr clabs, fpstr ctyps, int *ncol, int clabs_len, int ctyps_len),
               (int *mindx, fpstr clabs, fpstr ctyps, int *ncol),
               (int *mindx, fpstr clabs, int clabs_len, fpstr ctyps, int ctyps_len, int *ncol));

FORTRAN_SUBR ( LRCLID, lrclid,
               (int *mindx, int csetid[], int *ncol),
               (int *mindx, int csetid[], int *ncol),
               (int *mindx, int csetid[], int *ncol));

FORTRAN_SUBR ( LRCELL, lrcell,
               (int *mindx, float cell[]),
               (int *mindx, float cell[]),
               (int *mindx, float cell[]));

FORTRAN_SUBR ( LRRSOL, lrrsol,
               (int *mindx, float *minres, float *maxres),
               (int *mindx, float *minres, float *maxres),
               (int *mindx, float *minres, float *maxres));

FORTRAN_SUBR ( LRSYMI, lrsymi,
               (int *mindx, int *nsympx, fpstr ltypex, int *nspgrx, fpstr spgrnx,
                  fpstr pgnamx, int ltypex_len, int spgrnx_len, int pgnamx_len),
               (int *mindx, int *nsympx, fpstr ltypex, int *nspgrx, fpstr spgrnx,
                  fpstr pgnamx),
               (int *mindx, int *nsympx, fpstr ltypex, int ltypex_len, int *nspgrx,
                  fpstr spgrnx, int spgrnx_len, fpstr pgnamx, int pgnamx_len));

FORTRAN_SUBR ( LRSYMM, lrsymm,
               (int *mindx, int *nsymx, float rsymx[MAXSYM][4][4]),
               (int *mindx, int *nsymx, float rsymx[MAXSYM][4][4]),
               (int *mindx, int *nsymx, float rsymx[MAXSYM][4][4]));

FORTRAN_SUBR ( LKYIN, lkyin,
               (const int *mindx, const fpstr lsprgi, const int *nlprgi,
                   const int *ntok, const fpstr labin_line, const int ibeg[],
                   const int iend[], int lsprgi_len, int labin_line_len),
               (const int *mindx, const fpstr lsprgi, const int *nlprgi,
                   const int *ntok, const fpstr labin_line, const int ibeg[],
                   const int iend[]),
               (const int *mindx, const fpstr lsprgi, int lsprgi_len,
                   const int *nlprgi, const int *ntok, const fpstr labin_line,
                   int labin_line_len, const int ibeg[], const int iend[] ));

FORTRAN_SUBR ( LKYOUT, lkyout,
               (const int *mindx, const fpstr lsprgo, const int *nlprgo,
                   const int *ntok, const fpstr labin_line, const int ibeg[],
                   const int iend[], int lsprgo_len, int labin_line_len),
               (const int *mindx, const fpstr lsprgo, const int *nlprgo,
                   const int *ntok, const fpstr labin_line, const int ibeg[],
                   const int iend[]),
               (const int *mindx, const fpstr lsprgo, int lsprgo_len,
                   const int *nlprgo, const int *ntok, const fpstr labin_line,
                   int labin_line_len, const int ibeg[], const int iend[] ));

FORTRAN_SUBR ( LKYSET, lkyset,
               (const fpstr lsprgi, const int *nlprgi,
                   fpstr lsusrj, int kpoint[],
                   const int *itok, const int *ntok, const fpstr labin_line,
                   const int ibeg[], const int iend[],
                   int lsprgi_len, int lsusrj_len, int labin_line_len),
               (const fpstr lsprgi, const int *nlprgi,
                   fpstr lsusrj, int kpoint[],
                   const int *itok, const int *ntok, const fpstr labin_line,
                   const int ibeg[], const int iend[]),
               (const fpstr lsprgi, int lsprgi_len, const int *nlprgi,
                   fpstr lsusrj, int lsusrj_len, int kpoint[],
                   const int *itok, const int *ntok, const fpstr labin_line,
                   int labin_line_len, const int ibeg[], const int iend[] ));

FORTRAN_SUBR ( LRASSN, lrassn,
               (const int *mindx, fpstr lsprgi, int *nlprgi, int lookup[], fpstr ctprgi,
                      int lsprgi_len, int ctprgi_len),
               (const int *mindx, fpstr lsprgi, int *nlprgi, int lookup[], fpstr ctprgi),
               (const int *mindx, fpstr lsprgi, int lsprgi_len, int *nlprgi,
                      int lookup[], fpstr ctprgi, int ctprgi_len));

FORTRAN_SUBR ( LRIDX, lridx,
               (const int *mindx, fpstr project_name,
                  fpstr crystal_name, fpstr dataset_name,
                  int *isets, float *datcell, float *datwave,
                  int *ndatasets, int project_name_len,
                  int crystal_name_len, int dataset_name_len),
               (const int *mindx, fpstr project_name,
                  fpstr crystal_name, fpstr dataset_name,
                  int *isets, float *datcell, float *datwave,
                  int *ndatasets),
               (const int *mindx, fpstr project_name, int project_name_len,
                  fpstr crystal_name, int crystal_name_len,
                  fpstr dataset_name, int dataset_name_len,
                  int *isets, float *datcell, float *datwave,
                  int *ndatasets));

FORTRAN_SUBR ( LRCELX, lrcelx,
               (const int *mindx, const int *iset, float *mtzcell),
               (const int *mindx, const int *iset, float *mtzcell),
               (const int *mindx, const int *iset, float *mtzcell));

FORTRAN_SUBR ( LRIDC, lridc,
               (const int *mindx, fpstr project_name, fpstr dataset_name,
                  int *isets, float *datcell, float *datwave,
                  int *ndatasets, int project_name_len, int dataset_name_len),
               (const int *mindx, fpstr project_name, fpstr dataset_name,
                  int *isets, float *datcell, float *datwave,
                  int *ndatasets),
               (const int *mindx, fpstr project_name, int project_name_len,
                  fpstr dataset_name, int dataset_name_len,
                  int *isets, float *datcell, float *datwave,
                  int *ndatasets));

FORTRAN_SUBR ( LRID, lrid,
               (const int *mindx, fpstr project_name, fpstr dataset_name,
                  int *isets, int *ndatasets,
                  int project_name_len, int dataset_name_len),
               (const int *mindx, fpstr project_name, fpstr dataset_name,
                  int *isets, int *ndatasets),
               (const int *mindx, fpstr project_name, int project_name_len,
                  fpstr dataset_name, int dataset_name_len,
                  int *isets, int *ndatasets));

FORTRAN_SUBR ( LRSEEK, lrseek,
               (const int *mindx, int *nrefl),
               (const int *mindx, int *nrefl),
               (const int *mindx, int *nrefl));

FORTRAN_SUBR ( LRREFL, lrrefl,
               (const int *mindx, float *resol, float adata[], ftn_logical *eof),
               (const int *mindx, float *resol, float adata[], ftn_logical *eof),
               (const int *mindx, float *resol, float adata[], ftn_logical *eof));

FORTRAN_SUBR ( LRREFF, lrreff,
               (const int *mindx, float *resol, float adata[], ftn_logical *eof),
               (const int *mindx, float *resol, float adata[], ftn_logical *eof),
               (const int *mindx, float *resol, float adata[], ftn_logical *eof));

FORTRAN_SUBR ( LRREFM, lrrefm,
               (const int *mindx, ftn_logical logmiss[]),
               (const int *mindx, ftn_logical logmiss[]),
               (const int *mindx, ftn_logical logmiss[]));

FORTRAN_SUBR ( MTZ_CHECK_FOR_MNF, mtz_check_for_mnf,
               (const int *mindx, const int *ndata, float adata[], ftn_logical logmiss[]),
               (const int *mindx, const int *ndata, float adata[], ftn_logical logmiss[]),
               (const int *mindx, const int *ndata, float adata[], ftn_logical logmiss[]));

FORTRAN_SUBR ( LHPRT, lhprt,
               (const int *mindx, const int *iprint),
               (const int *mindx, const int *iprint),
               (const int *mindx, const int *iprint));

FORTRAN_SUBR ( LHPRT_ADV, lhprt_adv,
               (const int *mindx, const int *iprint),
               (const int *mindx, const int *iprint),
               (const int *mindx, const int *iprint));

FORTRAN_SUBR ( LRBAT, lrbat,
               (const int *mindx, int *batno, float rbatch[], fpstr cbatch,
                        const int *iprint, int cbatch_len),
               (const int *mindx, int *batno, float rbatch[], fpstr cbatch,
                        const int *iprint),
               (const int *mindx, int *batno, float rbatch[], fpstr cbatch,
                        int cbatch_len, const int *iprint));

FORTRAN_SUBR ( LBPRT, lbprt,
               (const int *ibatch, const int *iprint, float rbatch[], fpstr cbatch,
                        int cbatch_len),
               (const int *ibatch, const int *iprint, float rbatch[], fpstr cbatch),
               (const int *ibatch, const int *iprint, float rbatch[], fpstr cbatch,
                        int cbatch_len));

FORTRAN_SUBR ( LRBRES, lrbres,
               (const int *mindx, const int *batno),
               (const int *mindx, const int *batno),
               (const int *mindx, const int *batno));

FORTRAN_SUBR ( LRBTIT, lrbtit,
               (const int *mindx, const int *batno, fpstr tbatch,
                        const int *iprint, int tbatch_len),
               (const int *mindx, const int *batno, fpstr tbatch,
                        const int *iprint),
               (const int *mindx, const int *batno, fpstr tbatch,
                        int tbatch_len, const int *iprint));

FORTRAN_SUBR ( LRBSCL, lrbscl,
               (const int *mindx, const int *batno, float batscl[], int *nbatsc),
               (const int *mindx, const int *batno, float batscl[], int *nbatsc),
               (const int *mindx, const int *batno, float batscl[], int *nbatsc));

FORTRAN_SUBR ( LRBSETID, lrbsetid,
               (const int *mindx, const int *batno, int *bsetid),
               (const int *mindx, const int *batno, int *bsetid),
               (const int *mindx, const int *batno, int *bsetid));

FORTRAN_SUBR ( LRREWD, lrrewd,
               (const int *mindx),
               (const int *mindx),
               (const int *mindx));

FORTRAN_SUBR ( LSTRSL, lstrsl,
               (const int *mindx, const float *a, const float *b, const float *c,
                    const float *alpha, const float *beta, const float *gamma ),
               (const int *mindx, const float *a, const float *b, const float *c,
                    const float *alpha, const float *beta, const float *gamma ),
               (const int *mindx, const float *a, const float *b, const float *c,
                    const float *alpha, const float *beta, const float *gamma ));

FORTRAN_SUBR (LSTLSQ1, lstlsq1,
               (float *reso, const int *mindx, const int *ih, const int *ik, const int *il),
               (float *reso, const int *mindx, const int *ih, const int *ik, const int *il),
               (float *reso, const int *mindx, const int *ih, const int *ik, const int *il));

FORTRAN_SUBR ( LRCLOS, lrclos,
               (const int *mindx),
               (const int *mindx),
               (const int *mindx));

FORTRAN_SUBR ( LWOPEN_NOEXIT, lwopen_noexit,
               (const int *mindx, fpstr filename, int *ifail, int filename_len),
               (const int *mindx, fpstr filename, int *ifail),
               (const int *mindx, fpstr filename, int filename_len, int *ifail));

FORTRAN_SUBR ( LWOPEN, lwopen,
               (const int *mindx, fpstr filename, int filename_len),
               (const int *mindx, fpstr filename),
               (const int *mindx, fpstr filename, int filename_len));

FORTRAN_SUBR ( LWTITL, lwtitl,
               (const int *mindx, const fpstr ftitle, const int *flag, int ftitle_len),
               (const int *mindx, const fpstr ftitle, const int *flag),
               (const int *mindx, const fpstr ftitle, int ftitle_len, const int *flag));

FORTRAN_SUBR ( LWSORT, lwsort,
               (const int *mindx, int sortx[5]),
               (const int *mindx, int sortx[5]),
               (const int *mindx, int sortx[5]));

FORTRAN_SUBR ( LWHIST, lwhist,
               (int *mindx, fpstr hstrng, int *nlines, int hstrng_len),
               (int *mindx, fpstr hstrng, int *nlines),
               (int *mindx, fpstr hstrng, int hstrng_len, int *nlines));

FORTRAN_SUBR ( LWHSTL, lwhstl,
               (int *mindx, const fpstr hstrng, int hstrng_len),
               (int *mindx, const fpstr hstrng),
               (int *mindx, const fpstr hstrng, int hstrng_len));

FORTRAN_SUBR ( LWID, lwid,
               (const int *mindx, const fpstr project_name, const fpstr dataset_name,
                  int project_name_len, int dataset_name_len),
               (const int *mindx, const fpstr project_name, const fpstr dataset_name),
               (const int *mindx, const fpstr project_name, int project_name_len,
                  const fpstr dataset_name, int dataset_name_len));

FORTRAN_SUBR ( LWIDC, lwidc,
               (const int *mindx, const fpstr project_name, const fpstr dataset_name,
                  float datcell[6], float *datwave,
                  int project_name_len, int dataset_name_len),
               (const int *mindx, const fpstr project_name, const fpstr dataset_name,
                  float datcell[6], float *datwave),
               (const int *mindx, const fpstr project_name, int project_name_len,
                  const fpstr dataset_name, int dataset_name_len,
                  float datcell[6], float *datwave));

FORTRAN_SUBR ( LWIDX, lwidx,
               (const int *mindx, const fpstr project_name, const fpstr crystal_name,
                  const fpstr dataset_name, float datcell[6], float *datwave,
                  int project_name_len, int crystal_name_len, int dataset_name_len),
               (const int *mindx, const fpstr project_name, const fpstr crystal_name,
                  const fpstr dataset_name, float datcell[6], float *datwave),
               (const int *mindx, const fpstr project_name, int project_name_len,
                  const fpstr crystal_name, int crystal_name_len,
                  const fpstr dataset_name, int dataset_name_len,
                  float datcell[6], float *datwave));

FORTRAN_SUBR ( LWCELL, lwcell,
               (const int *mindx, float cell[6]),
               (const int *mindx, float cell[6]),
               (const int *mindx, float cell[6]));

FORTRAN_SUBR ( LWIDAS, lwidas,
               (const int *mindx, int *nlprgo, fpstr pname, fpstr dname, const int *iappnd,
                      int pname_len, int dname_len),
               (const int *mindx, int *nlprgo, fpstr pname, fpstr dname, const int *iappnd),
               (const int *mindx, int *nlprgo, fpstr pname, int pname_len,
                      fpstr dname, int dname_len, const int *iappnd));

FORTRAN_SUBR ( LWIDASX, lwidasx,
               (const int *mindx, int *nlprgo, fpstr xname, fpstr dname, int *iappnd,
                      int xname_len, int dname_len),
               (const int *mindx, int *nlprgo, fpstr xname, fpstr dname, int *iappnd),
               (const int *mindx, int *nlprgo, fpstr xname, int xname_len,
                      fpstr dname, int dname_len, int *iappnd));

FORTRAN_SUBR ( LWIDALL, lwidall,
               (const int *mindx, fpstr xname, fpstr dname,
                      int xname_len, int dname_len),
               (const int *mindx, fpstr xname, fpstr dname),
               (const int *mindx, fpstr xname, int xname_len,
                      fpstr dname, int dname_len));

FORTRAN_SUBR ( LWSYMM, lwsymm,
               (int *mindx, int *nsymx, int *nsympx, float rsymx[MAXSYM][4][4],
                  fpstr ltypex, int *nspgrx, fpstr spgrnx, fpstr pgnamx,
                  int ltypex_len, int spgrnx_len, int pgnamx_len),
               (int *mindx, int *nsymx, int *nsympx, float rsymx[MAXSYM][4][4],
                  fpstr ltypex, int *nspgrx, fpstr spgrnx, fpstr pgnamx),
               (int *mindx, int *nsymx, int *nsympx, float rsymx[MAXSYM][4][4],
                  fpstr ltypex, int ltypex_len, int *nspgrx, fpstr spgrnx,
                int spgrnx_len, fpstr pgnamx, int pgnamx_len));

FORTRAN_SUBR ( LWASSN, lwassn,
               (const int *mindx, fpstr lsprgo, const int *nlprgo, fpstr ctprgo, int *iappnd,
                      int lsprgo_len, int ctprgo_len),
               (const int *mindx, fpstr lsprgo, const int *nlprgo, fpstr ctprgo, int *iappnd),
               (const int *mindx, fpstr lsprgo, int lsprgo_len, const int *nlprgo,
                      fpstr ctprgo, int ctprgo_len, int *iappnd));

FORTRAN_SUBR ( LWCLAB, lwclab,
               (const int *mindx, fpstr lsprgo, const int *nlprgo, fpstr ctprgo, const int *iappnd,
                      int lsprgo_len, int ctprgo_len),
               (const int *mindx, fpstr lsprgo, const int *nlprgo, fpstr ctprgo, const int *iappnd),
               (const int *mindx, fpstr lsprgo, int lsprgo_len, const int *nlprgo,
                      fpstr ctprgo, int ctprgo_len, const int *iappnd));

FORTRAN_SUBR ( LWBAT, lwbat,
               (const int *mindx, int *batno, float rbatch[], fpstr cbatch,
                        int cbatch_len),
               (const int *mindx, int *batno, float rbatch[], fpstr cbatch),
               (const int *mindx, int *batno, float rbatch[], fpstr cbatch,
                        int cbatch_len));

FORTRAN_SUBR ( LWBTIT, lwbtit,
               (const int *mindx, int *batno, fpstr tbatch, int tbatch_len),
               (const int *mindx, int *batno, fpstr tbatch),
               (const int *mindx, int *batno, fpstr tbatch, int tbatch_len));

FORTRAN_SUBR ( LWBSCL, lwbscl,
               (const int *mindx, int *batno, float batscl[], int *nbatsc),
               (const int *mindx, int *batno, float batscl[], int *nbatsc),
               (const int *mindx, int *batno, float batscl[], int *nbatsc));

FORTRAN_SUBR ( LWBSETID, lwbsetid,
               (const int *mindx, const int *batno, const fpstr project_name,
                  const fpstr dataset_name,
                  int project_name_len, int dataset_name_len),
               (const int *mindx, const int *batno, const fpstr project_name,
                  const fpstr dataset_name),
               (const int *mindx, const int *batno, const fpstr project_name,
                  int project_name_len,
                  const fpstr dataset_name, int dataset_name_len));

FORTRAN_SUBR ( LWBSETIDX, lwbsetidx,
               (const int *mindx, const int *batno, const fpstr crystal_name,
                  const fpstr dataset_name,
                  int crystal_name_len, int dataset_name_len),
               (const int *mindx, const int *batno, const fpstr crystal_name,
                  const fpstr dataset_name),
               (const int *mindx, const int *batno, const fpstr crystal_name,
                  int crystal_name_len,
                  const fpstr dataset_name, int dataset_name_len));

FORTRAN_SUBR ( EQUAL_MAGIC, equal_magic,
               (const int *mindx, float adata[], const int *ncol),
               (const int *mindx, float adata[], const int *ncol),
               (const int *mindx, float adata[], const int *ncol));

FORTRAN_SUBR ( SET_MAGIC, set_magic,
               (const int *mindx, float *val_magic, ftn_logical *setval),
               (const int *mindx, float *val_magic, ftn_logical *setval),
               (const int *mindx, float *val_magic, ftn_logical *setval));

FORTRAN_SUBR ( RESET_MAGIC, reset_magic,
               (const int *mindx, const float adata[], float bdata[],
                 const int *ncol, const float *val_magica, const float *val_magicb),
               (const int *mindx, const float adata[], float bdata[],
                 const int *ncol, const float *val_magica, const float *val_magicb),
               (const int *mindx, const float adata[], float bdata[],
                 const int *ncol, const float *val_magica, const float *val_magicb));

FORTRAN_SUBR ( LWREFL_NOEXIT, lwrefl_noexit,
               (const int *mindx, const float adata[], int *ifail),
               (const int *mindx, const float adata[], int *ifail),
               (const int *mindx, const float adata[], int *ifail));

FORTRAN_SUBR ( LWREFL, lwrefl,
               (const int *mindx, const float adata[]),
               (const int *mindx, const float adata[]),
               (const int *mindx, const float adata[]));

FORTRAN_SUBR ( LWCLOS_NOEXIT, lwclos_noexit,
               (const int *mindx, const int *iprint, int *ifail),
               (const int *mindx, const int *iprint, int *ifail),
               (const int *mindx, const int *iprint, int *ifail));

FORTRAN_SUBR ( LWCLOS, lwclos,
               (const int *mindx, const int *iprint),
               (const int *mindx, const int *iprint),
               (const int *mindx, const int *iprint));

FORTRAN_SUBR ( RBATHD, rbathd,
               (),(),());

FORTRAN_SUBR ( WBATHD, wbathd,
               (),(),());

FORTRAN_SUBR ( LRHDRL, lrhdrl,
               (),(),());

FORTRAN_SUBR ( LABPRT, labprt,
               (),(),());

FORTRAN_SUBR ( LBPRTH, lbprth,
               (),(),());

FORTRAN_SUBR ( SORTUP, sortup,
               (),(),());

FORTRAN_SUBR ( ADDLIN, addlin,
               (),(),());

FORTRAN_FUN (int, NEXTLN, nextln,
               (),(),());

FORTRAN_SUBR ( IS_MAGIC, is_magic,
               (const float *val_magic, const float *valtst, ftn_logical *lvalms),
               (const float *val_magic, const float *valtst, ftn_logical *lvalms),
               (const float *val_magic, const float *valtst, ftn_logical *lvalms));

#ifdef __cplusplus
} }
#endif

#endif // GUARD
