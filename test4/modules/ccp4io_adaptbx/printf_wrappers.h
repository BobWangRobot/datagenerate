#ifndef CCP4IO_ADAPTBX_PRINTF_WRAPPERS_H
#define CCP4IO_ADAPTBX_PRINTF_WRAPPERS_H

#include <stdio.h>

#ifdef  __cplusplus
extern "C" {
#endif

extern int ccp4io_printf_mode;

void
ccp4io_printf(
  const char* format,
  ...)
;

void
ccp4io_fprintf(
  FILE* stream,
  const char* format,
  ...)
;

#ifdef __cplusplus
}
#endif

#endif /* GUARD */
