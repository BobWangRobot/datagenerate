#include <ccp4io_adaptbx/printf_wrappers.h>
#include <stdarg.h>

int ccp4io_printf_mode = 1;

void
ccp4io_printf(
  const char* format,
  ...)
{
  va_list arg;
  va_start(arg, format);
  if (ccp4io_printf_mode != 0) {
    vprintf(format, arg);
  }
  va_end(arg);
}

void
ccp4io_fprintf(
  FILE* stream,
  const char* format,
  ...)
{
  va_list arg;
  va_start(arg, format);
  if (ccp4io_printf_mode != 0) {
    vprintf(format, arg);
  }
  va_end(arg);
}
