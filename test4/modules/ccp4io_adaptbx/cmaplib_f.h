#ifndef CCP4IO_ADAPTBX_CMAPLIB_F_H
#define CCP4IO_ADAPTBX_CMAPLIB_F_H

#include<string.h>
#include<stdio.h>
#include<fcntl.h>
#include"cmaplib_f.h"
#include"cmap_errno.h"
#include"ccp4_fortran.h"
#include"ccp4_parser.h"
#include"csymlib.h"
#include"ccp4_general.h"

#ifdef  __cplusplus
namespace CMap_io {
extern "C" {
#endif

FORTRAN_SUBR ( MWRHDL, mwrhdl,
    (int *iunit, const fpstr mapnam, const fpstr title, int *nsecs, int iuvw[3],
     int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2,
     float cell[6], int *lspgrp, int *lmode, int mapnam_len, int title_len),
    (int *iunit, const fpstr mapnam, const fpstr title, int *nsecs, int iuvw[3],
     int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2,
     float cell[6], int *lspgrp, int *lmode),
    (int *iunit, const fpstr mapnam, int mapnam_len, const fpstr title,
     int title_len, int *nsecs, int iuvw[3], int mxyz[3], int *nw1, int *nu1,
     int *nu2, int *nv1, int *nv2, float cell[6], int *lspgrp, int *lmode));

FORTRAN_SUBR ( CCP4_MAP_WRITE_OPEN_HEADER_BY_NAME,
               ccp4_map_write_open_header_by_name,
    (int *iunit, const fpstr mapnam, const fpstr title, int *nsecs, int iuvw[3],
     int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2,
     float cell[6], int *lspgrp, int *lmode, int mapnam_len, int title_len),
    (int *iunit, const fpstr mapnam, const fpstr title, int *nsecs, int iuvw[3],
     int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2,
     float cell[6], int *lspgrp, int *lmode),
    (int *iunit, const fpstr mapnam, int mapnam_len, const fpstr title,
     int title_len, int *nsecs, int iuvw[3], int mxyz[3], int *nw1, int *nu1,
     int *nu2, int *nv1, int *nv2, float cell[6], int *lspgrp, int *lmode));;

FORTRAN_SUBR ( MWRHDR, mwrhdr,
    (int *iunit, const fpstr title, int *nsecs, int iuvw[3],
     int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2,
     float cell[6], int *lspgrp, int *lmode, int title_len),
    (int *iunit, const fpstr title, int *nsecs, int iuvw[3],
     int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2,
     float cell[6], int *lspgrp, int *lmode),
    (int *iunit, const fpstr title,
     int title_len, int *nsecs, int iuvw[3], int mxyz[3], int *nw1, int *nu1,
     int *nu2, int *nv1, int *nv2, float cell[6], int *lspgrp, int *lmode));

FORTRAN_SUBR ( CCP4_MAP_WRITE_OPEN_HEADER_BY_ID,
               ccp4_map_write_open_header_by_id,
    (int *iunit, const fpstr title, int *nsecs, int iuvw[3],
     int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2,
     float cell[6], int *lspgrp, int *lmode, int title_len),
    (int *iunit, const fpstr title, int *nsecs, int iuvw[3],
     int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2,
     float cell[6], int *lspgrp, int *lmode),
    (int *iunit, const fpstr title,
     int title_len, int *nsecs, int iuvw[3], int mxyz[3], int *nw1, int *nu1,
     int *nu2, int *nv1, int *nv2, float cell[6], int *lspgrp, int *lmode));

FORTRAN_SUBR( MRDHDS, mrdhds,
              (int *iunit, const fpstr mapnam, fpstr title,
               int *nsec, int iuvw[3], int mxyz[3], int *nw1, int *nu1,
               int *nu2, int *nv1, int *nv2, float cell[6], int *lspgrp,
               int *lmode, float *rhmin, float *rhmax, float *rhmean,
               float *rhrms, int *ifail, int *iprint, int mapnam_len,
               int title_len),
              (int *iunit, const fpstr mapnam, fpstr title, int *nsec,
               int iuvw[3], int mxyz[3], int *nw1, int *nu1, int *nu2,
               int *nv1, int *nv2, float cell[6], int *lspgrp, int *lmode,
               float *rhmin, float *rhmax, float *rhmean, float * rhrms,
               int *ifail, int *iprint),
              (int *iunit, const fpstr mapnam, int mapnam_len,
               fpstr title, int title_len, int *nsec, int iuvw[3],
               int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2,
               float cell[6], int *lspgrp, int *lmode, float *rhmin,
               float *rhmax, float *rhmean, float * rhrms, int *ifail,
               int *iprint));

FORTRAN_SUBR( CCP4_MAP_READ_OPEN_HEADER_CHECK,
              ccp4_map_read_open_header_check,
              (int *iunit, const fpstr mapnam, fpstr title,
               int *nsec, int iuvw[3], int mxyz[3], int *nw1, int *nu1,
               int *nu2, int *nv1, int *nv2, float cell[6], int *lspgrp,
               int *lmode, float *rhmin, float *rhmax, float *rhmean,
               float *rhrms, int *ifail, int *iprint, int mapnam_len,
               int title_len),
              (int *iunit, const fpstr mapnam, fpstr title, int *nsec,
               int iuvw[3], int mxyz[3], int *nw1, int *nu1, int *nu2,
               int *nv1, int *nv2, float cell[6], int *lspgrp, int *lmode,
               float *rhmin, float *rhmax, float *rhmean, float * rhrms,
               int *ifail, int *iprint),
              (int *iunit, const fpstr mapnam, int mapnam_len,
               fpstr title, int title_len, int *nsec, int iuvw[3],
               int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2,
               float cell[6], int *lspgrp, int *lmode, float *rhmin,
               float *rhmax, float *rhmean, float * rhrms, int *ifail,
               int *iprint));

FORTRAN_SUBR( MRDHDR, mrdhdr,
              (int *iunit, const fpstr mapnam, fpstr title,
               int *nsec, int iuvw[3], int mxyz[3], int *nw1, int *nu1,
               int *nu2, int *nv1, int *nv2, float cell[6], int *lspgrp,
               int *lmode, float *rhmin, float *rhmax, float *rhmean,
               float *rhrms, int mapnam_len, int title_len),
              (int *iunit, const fpstr mapnam, fpstr title, int *nsec,
               int iuvw[3], int mxyz[3], int *nw1, int *nu1, int *nu2,
               int *nv1, int *nv2, float cell[6], int *lspgrp, int *lmode,
               float *rhmin, float *rhmax, float *rhmean, float * rhrms),
              (int *iunit, const fpstr mapnam, int mapnam_len,
               fpstr title, int title_len, int *nsec, int iuvw[3],
               int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2,
               float cell[6], int *lspgrp, int *lmode, float *rhmin,
               float *rhmax, float *rhmean, float * rhrms));

FORTRAN_SUBR( CCP4_MAP_READ_OPEN_HEADER,
              ccp4_map_read_open_header,
              (int *iunit, const fpstr mapnam, fpstr title,
               int *nsec, int iuvw[3], int mxyz[3], int *nw1, int *nu1,
               int *nu2, int *nv1, int *nv2, float cell[6], int *lspgrp,
               int *lmode, float *rhmin, float *rhmax, float *rhmean,
               float *rhrms, int mapnam_len, int title_len),
              (int *iunit, const fpstr mapnam, fpstr title, int *nsec,
               int iuvw[3], int mxyz[3], int *nw1, int *nu1, int *nu2,
               int *nv1, int *nv2, float cell[6], int *lspgrp, int *lmode,
               float *rhmin, float *rhmax, float *rhmean, float * rhrms),
              (int *iunit, const fpstr mapnam, int mapnam_len,
               fpstr title, int title_len, int *nsec, int iuvw[3],
               int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2,
               float cell[6], int *lspgrp, int *lmode, float *rhmin,
               float *rhmax, float *rhmean, float * rhrms));

FORTRAN_SUBR( MWCLOSE, mwclose, (int *iunit), (int *iunit), (int *iunit));

FORTRAN_SUBR( CCP4_MAP_WRITE_CLOSE_AUTO,
              ccp4_map_write_close_auto,
              (int *iunit), (int *iunit), (int *iunit));

FORTRAN_SUBR( MCLOSE, mclose,
              (int *iunit, float *min, float *max, float *mean, float *rms),
              (int *iunit, float *min, float *max, float *mean, float *rms),
              (int *iunit, float *min, float *max, float *mean, float *rms));

FORTRAN_SUBR( CCP4_MAP_WRITE_CLOSE_USER_SUM,
              ccp4_map_write_close_user_sum,
              (int *iunit, float *min, float *max, float *mean, float *rms),
              (int *iunit, float *min, float *max, float *mean, float *rms),
              (int *iunit, float *min, float *max, float *mean, float *rms));

FORTRAN_SUBR( MCLOSC, mclosc,
              (int *iunit, float *min, float *max, float *mean, float *rms),
              (int *iunit, float *min, float *max, float *mean, float *rms),
              (int *iunit, float *min, float *max, float *mean, float *rms));

FORTRAN_SUBR( CCP4_MAP_WRITE_CLOSE_USER_MEAN,
              ccp4_map_write_close_user_mean,
              (int *iunit, float *min, float *max, float *mean, float *rms),
              (int *iunit, float *min, float *max, float *mean, float *rms),
              (int *iunit, float *min, float *max, float *mean, float *rms));

FORTRAN_SUBR( MRCLOS, mrclos, (int *iunit), (int *iunit), (int *iunit));

FORTRAN_SUBR( CCP4_MAP_READ_CLOSE,
              ccp4_map_read_close,
              (int *iunit), (int *iunit), (int *iunit));

FORTRAN_SUBR( MSYWRT, msywrt,
              (int *iunit, int *nsym, float *rot),
              (int *iunit, int *nsym, float *rot),
              (int *iunit, int *nsym, float *rot));

FORTRAN_SUBR( CCP4_MAP_WRITE_SYMM_MATRIX, ccp4_map_write_symm_matrix,
              (int *iunit, int *nsym, float *rot),
              (int *iunit, int *nsym, float *rot),
              (int *iunit, int *nsym, float *rot));

FORTRAN_SUBR( MSYMOP, msymop,
              (int *iunit, int *nsym, float *rot),
              (int *iunit, int *nsym, float *rot),
              (int *iunit, int *nsym, float *rot));

FORTRAN_SUBR( CCP4_MAP_READ_SYMM_MATRIX,
              ccp4_map_read_symm_matrix,
              (int *iunit, int *nsym, float *rot),
              (int *iunit, int *nsym, float *rot),
              (int *iunit, int *nsym, float *rot));

FORTRAN_SUBR( MSYPUT, msyput,
              (int *sunit, int *spacegroup, int *iunit),
              (int *sunit, int *spacegroup, int *iunit),
              (int *sunit, int *spacegroup, int *iunit));

FORTRAN_SUBR( CCP4_MAP_WRITE_SPGNAME,
              ccp4_map_write_spgname,
              (int *sunit, int *spacegroup, int *iunit),
              (int *sunit, int *spacegroup, int *iunit),
              (int *sunit, int *spacegroup, int *iunit));

FORTRAN_SUBR( MSPEW, mspew,
              (int *iunit, float *section),
              (int *iunit, float *section),
              (int *iunit, float *section));

FORTRAN_SUBR( CCP4_MAP_WRITE_ALL_SECTION,
              ccp4_map_write_all_section,
              (int *iunit, float *section),
              (int *iunit, float *section),
              (int *iunit, float *section));

FORTRAN_SUBR( MGULP, mgulp,
              (int *iunit, float *section, int *ier),
              (int *iunit, float *section, int *ier),
              (int *iunit, float *section, int *ier));

FORTRAN_SUBR( CCP4_MAP_READ_WHOLE_SECTION_AS_MODE,
              ccp4_map_read_whole_section_as_mode,
              (int *iunit, float *section, int *ier),
              (int *iunit, float *section, int *ier),
              (int *iunit, float *section, int *ier));

FORTRAN_SUBR( MWRSEC, mwrsec,
              (int *iunit, float *section, int *mu, int *mv, int *iu1,
               int *iu2, int *iv1, int *iv2),
              (int *iunit, float *section, int *mu, int *mv, int *iu1,
               int *iu2, int *iv1, int *iv2),
              (int *iunit, float *section, int *mu, int *mv, int *iu1,
               int *iu2, int *iv1, int *iv2));

FORTRAN_SUBR( CCP4_MAP_WRITE_PART_SECTION,
              ccp4_map_write_part_section,
              (int *iunit, float *section, int *mu, int *mv, int *iu1,
               int *iu2, int *iv1, int *iv2),
              (int *iunit, float *section, int *mu, int *mv, int *iu1,
               int *iu2, int *iv1, int *iv2),
              (int *iunit, float *section, int *mu, int *mv, int *iu1,
               int *iu2, int *iv1, int *iv2));

FORTRAN_SUBR( MRDLIN, mrdlin,
              (int *iunit, float *line, int *ier),
              (int *iunit, float *line, int *ier),
              (int *iunit, float *line, int *ier));

FORTRAN_SUBR( CCP4_MAP_READ_LINE_AS_MODE,
              ccp4_map_read_line_as_mode,
              (int *iunit, float *line, int *ier),
              (int *iunit, float *line, int *ier),
              (int *iunit, float *line, int *ier));

FORTRAN_SUBR( MGULPR, mgulpr,
              (int *iunit, float *section, int *ier),
              (int *iunit, float *section, int *ier),
              (int *iunit, float *section, int *ier));

FORTRAN_SUBR( CCP4_MAP_READ_WHOLE_SECT_AS_REAL,
              ccp4_map_read_whole_sect_as_real,
              (int *iunit, float *section, int *ier),
              (int *iunit, float *section, int *ier),
              (int *iunit, float *section, int *ier));

FORTRAN_SUBR( MPOSN, mposn,
              (int *iunit, int *sec),
              (int *iunit, int *sec),
              (int *iunit, int *sec));

FORTRAN_SUBR( CCP4_MAP_READ_POSITION_SELECTION,
              ccp4_map_read_position_selection,
              (int *iunit, int *sec),
              (int *iunit, int *sec),
              (int *iunit, int *sec));

FORTRAN_SUBR( MPOSNW, mposnw,
              (int *iunit, int *sec),
              (int *iunit, int *sec),
              (int *iunit, int *sec));

FORTRAN_SUBR( CCP4_MAP_WRITE_POSITION_SECTION,
              ccp4_map_write_position_section,
              (int *iunit, int *sec),
              (int *iunit, int *sec),
              (int *iunit, int *sec));

FORTRAN_SUBR( MSYCPY, msycpy,
              (int *iunit, int *ounit),
              (int *iunit, int *ounit),
              (int *iunit, int *ounit));

FORTRAN_SUBR( CCP4_MAP_COPY_SYMMETRY,
              ccp4_map_copy_symmetry,
              (int *iunit, int *ounit),
              (int *iunit, int *ounit),
              (int *iunit, int *ounit));

FORTRAN_SUBR( MTTREP, mttrep,
              (const fpstr label, int *posn, int label_len),
              (const fpstr label, int *posn),
              (const fpstr label, int label_len, int *posn));

FORTRAN_SUBR( CCP4_MAP_WRITE_REPLACE_TITLE,
              ccp4_map_write_replace_title,
              (const fpstr label, int *posn, int label_len),
              (const fpstr label, int *posn),
              (const fpstr label, int label_len, int *posn));

FORTRAN_SUBR( MTTCPY, mttcpy,
              (const fpstr label, int label_len),
              (const fpstr label),
              (const fpstr label, int label_len));

FORTRAN_SUBR( CCP4_MAP_COPY_TITLE,
              ccp4_map_copy_title,
              (const fpstr label, int label_len),
              (const fpstr label),
              (const fpstr label, int label_len));

FORTRAN_SUBR( MSKPUT, mskput,
              (float *skew_rot, float *skew_trans),
              (float *skew_rot, float *skew_trans),
              (float *skew_rot, float *skew_trans));
;

FORTRAN_SUBR( CCP4_MAP_WRITE_SKEW_INFO,
              ccp4_map_write_skew_info,
              (float *skew_rot, float *skew_trans),
              (float *skew_rot, float *skew_trans),
              (float *skew_rot, float *skew_trans));

FORTRAN_FUN(int, MSKGET, mskget,
            (float *mask_rot, float *mask_trans),
            (float *mask_rot, float *mask_trans),
            (float *mask_rot, float *mask_trans));

FORTRAN_SUBR( CCP4_MAP_WRITE_SECTION_HEADER,
              ccp4_map_write_section_header,
              (int *iunit, float *section, const fpstr local_hdr, int local_hdr_len),
              (int *iunit, float *section, const fpstr local_hdr),
              (int *iunit, float *section, const fpstr local_hdr, int local_hdr_len));

FORTRAN_SUBR( CCP4_MAP_READ_SECTION_HEADER,
              ccp4_map_read_section_header,
              (int *iunit, float *section, const fpstr local_hdr,
               int *ier, int local_hdr_len),
              (int *iunit, float *section, const fpstr local_hdr,
               int *ier),
              (int *iunit, float *section, const fpstr local_hdr,
               int local_hdr_len, int *ier));

FORTRAN_SUBR( CCP4_MAP_SET_LOCAL_HEADER,
              ccp4_map_set_local_header,
              (int *iunit, int *size),
              (int *iunit, int *size),
              (int *iunit, int *size));

FORTRAN_SUBR( CCP4_MAP_GET_LOCAL_HEADER,
              ccp4_map_get_local_header,
              (int *iunit, int *size),
              (int *iunit, int *size),
              (int *iunit, int *size));

FORTRAN_SUBR( MWFNAM, mwfnam,
              (fpstr logname, int logname_len),
              (fpstr logname),
              (fpstr logname, int logname_len));

FORTRAN_SUBR( CCP4_MAP_GET_LAST_WRITE_FILENAME,
              ccp4_map_get_last_write_filename,
              (fpstr logname, int logname_len),
              (fpstr logname),
              (fpstr logname, int logname_len));

FORTRAN_SUBR( MRFNAM, mrfnam,
              (fpstr logname, int logname_len),
              (fpstr logname),
              (fpstr logname, int logname_len));

FORTRAN_SUBR( CCP4_MAP_GET_LAST_READ_FILENAME,
              ccp4_map_get_last_read_filename,
              (fpstr logname, int logname_len),
              (fpstr logname),
              (fpstr logname, int logname_len));

FORTRAN_SUBR( MSTMST, mstmst, (int *map), (int *map), (int *map));

FORTRAN_SUBR( MODECV, modecv,
              (float *output, float *input, int *number, int *mode),
              (float *output, float *input, int *number, int *mode),
              (float *output, float *input, int *number, int *mode));

FORTRAN_SUBR( CCP4_MAP_MODE_TO_REAL,
              ccp4_map_mode_to_real,
              (float *output, float *input, int *number, int *mode),
              (float *output, float *input, int *number, int *mode),
              (float *output, float *input, int *number, int *mode));

FORTRAN_FUN(int, NBYTXX, nbytxx,  (int *nwords), (int *nwords),(int *nwords));

#ifdef __cplusplus
} }
#endif

#endif // GUARD
