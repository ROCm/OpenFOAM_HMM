#ifndef MAKEDEPEND
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#endif

#include "cgnslib.h"
#include "cgns_header.h"
#include "libadf/ADF.h"

char cgns_error_mess[200] = "no CGNS error reported";

void cgi_error(char *format, ...) {
    va_list arg;
    va_start(arg, format);
    vsprintf(cgns_error_mess,format, arg);
    va_end(arg);
}

void cgi_warning(char *format, ...) {
    va_list arg;
    fprintf(stdout,"*** Warning:");
    va_start(arg, format);
    vfprintf(stdout,format,arg);
    va_end(arg);
    fprintf(stdout," ***\n");
}

char const *cg_get_error() {
    return cgns_error_mess;
}

void cg_error_exit() {
    fprintf(stderr,"%s\n",cgns_error_mess);
    exit(1);
}

void cg_error_print() {
    fprintf(stderr,"%s\n",cgns_error_mess);
}

void adf_error(char *routine_name, int ier) {
    char adf_err[80];

    ADF_Error_Message(ier, adf_err);
    cgi_error("Error in routine '%s':\n '%s'",routine_name,adf_err);
}

