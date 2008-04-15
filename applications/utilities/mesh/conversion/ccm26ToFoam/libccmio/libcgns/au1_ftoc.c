#ifndef MAKEDEPEND
#include <stdlib.h>
#include <string.h>
#endif

#include "fortran_macros.h"
#include "libadf/ADF_fbind.h"
#include "ADF_U1_internals.h"

void FMNAME(au1cn, AU1CN) (double *id, STR_PSTR(name), STR_PSTR(label),
	int *depth, int *numnodes, int *ier STR_PLEN(name) STR_PLEN(label)) {
	int len_name = (int)STR_LEN(name);
	int len_label = (int)STR_LEN(label);
	FMCALL(au1cn2, AU1CN2)(id, STR_PTR(name), &len_name, STR_PTR(label), &len_label,
        	depth, numnodes, ier);
}

void FMNAME(au1rnid, AU1RNID) (double *id, STR_PSTR(name), STR_PSTR(label),
	int *istart, int *ilen, int *depth, int *ilen_ret, double *node_id,
        int *ier STR_PLEN(name) STR_PLEN(label)) {
	int len_name = (int)STR_LEN(name);
	int len_label = (int)STR_LEN(label);
	FMCALL(au1rnid2, AU1RNID2)(id, STR_PTR(name), &len_name, STR_PTR(label), &len_label,
        	istart, ilen, depth, ilen_ret, node_id, ier);
}

void FMNAME(au1can, AU1CAN) (double *id, STR_PSTR(name), STR_PSTR(label),
	int *nancstr, int *ier STR_PLEN(name) STR_PLEN(label)) {
	int len_name = (int)STR_LEN(name);
	int len_label = (int)STR_LEN(label);
	FMCALL(au1can2, AU1CAN2)(id, STR_PTR(name), &len_name, STR_PTR(label), &len_label,
        	nancstr, ier);
}

void FMNAME(au1raid, AU1RAID) (double *id, STR_PSTR(name), STR_PSTR(label),
	int *istart, int *ilen, int *ilen_ret, double *node_id,
        int *ier STR_PLEN(name) STR_PLEN(label)) {
	int len_name = (int)STR_LEN(name);
	int len_label = (int)STR_LEN(label);
	FMCALL(au1raid2, AU1RAID2)(id, STR_PTR(name), &len_name, STR_PTR(label), &len_label,
        	istart, ilen, ilen_ret, node_id, ier);
}

