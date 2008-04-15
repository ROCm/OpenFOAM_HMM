#ifndef MAKEDEPEND
#include <stdlib.h>
#include <string.h>
#endif

#include "fortran_macros.h"
#include "libadf/ADF_internals.h"

/* ADF Fortran to C interface */

void FMNAME(adfcnam, ADFCNAM) (double *ID, int *istart, int *inum,
	int *inamlen, int *inumret, STR_PSTR(names),
 	int *ier STR_PLEN(names)) {
	int len_names = (int)STR_LEN(names);
	FMCALL(adfcna2,ADFCNA2)(ID, istart, inum, inamlen, &len_names,
		inumret, STR_PTR(names), ier);
}

void FMNAME(adfcids, ADFCIDS) (double *ID, int *istart, int *inum,
	int *inumret, double *cIDs, int *ier) {
	FMCALL(adfcid2,ADFCID2)(ID, istart, inum, inumret, cIDs, ier);
}

void FMNAME(adfcre, ADFCRE) (double *PID, STR_PSTR(name), double *ID,
	int *ier STR_PLEN(name)) {
	int len_name = (int)STR_LEN(name);
	FMCALL(adfcre2, ADFCRE2)(PID, STR_PTR(name), &len_name, ID, ier);
}

void FMNAME(adfdclo, ADFDCLO) (double *RootID, int *ier) {
	FMCALL(adfdcl2, ADFDCL2)(RootID, ier);
}

void FMNAME(adfddel, ADFDDEL) (STR_PSTR(filename),
	int *ier STR_PLEN(filename)) {
	int len_filename = (int)STR_LEN(filename);
	FMCALL(adfdde2, ADFDDE2)(STR_PTR(filename), &len_filename, ier);
}

void FMNAME(adfdel, ADFDEL) (double *PID, double *ID, int *ier) {
	FMCALL(adfdel2, ADFDEL2)(PID, ID, ier);
}

void FMNAME(adfdgc, ADFDGC) (double *ID, int *ier) {
	FMCALL(adfdgc2, ADFDGC2)(ID, ier);
}

void FMNAME(adfdgf, ADFDGF) (double *RootID, STR_PSTR(format),
	int *ier STR_PLEN(format)) {
	int len_format = (int)STR_LEN(format);
	FMCALL(adfdgf2, ADFDGF2)(RootID, STR_PTR(format), &len_format, ier);
}

void FMNAME(adfdopn, ADFDOPN) (STR_PSTR(filename), STR_PSTR(status),
	STR_PSTR(format), double *RootID, int *ier STR_PLEN(filename)
 	STR_PLEN(status) STR_PLEN(format)) {
	int len_filename = (int)STR_LEN(filename);
	int len_status = (int)STR_LEN(status);
	int len_format = (int)STR_LEN(format);
	FMCALL(adfdop2, ADFDOP2)(STR_PTR(filename), &len_filename,
        	STR_PTR(status), &len_status, STR_PTR(format),
                &len_format, RootID, ier);
}

void FMNAME(adfdsf, ADFDSF) (double *RootID, STR_PSTR(format),
	int *ier STR_PLEN(format)) {
	int len_format = (int)STR_LEN(format);
	FMCALL(adfdsf2, ADFDSF2)(RootID, STR_PTR(format), &len_format, ier);
}

void FMNAME(adfdver, ADFDVER) (double *RootID, STR_PSTR(version),
	STR_PSTR(cdate), STR_PSTR(mdate), int *ier STR_PLEN(version)
 	STR_PLEN(cdate) STR_PLEN(mdate)) {
	int len_version = (int)STR_LEN(version);
	int len_cdate = (int)STR_LEN(cdate);
	int len_mdate = (int)STR_LEN(mdate);
	FMCALL(adfdve2, ADFDVE2)(RootID, STR_PTR(version), STR_PTR(cdate),
        	STR_PTR(mdate), &len_version, &len_cdate, &len_mdate, ier);
}

void FMNAME(adferr, ADFERR) (int *ier, STR_PSTR(errstr) STR_PLEN(errstr)) {
	int len_errstr = (int)STR_LEN(errstr);
	FMCALL(adferr2, ADFERR2)(ier, STR_PTR(errstr), &len_errstr);
}

void FMNAME(adfftd, ADFFTD) (double *ID, int *ier) {
	FMCALL(adfftd2, ADFFTD2)(ID, ier);
}

void FMNAME(adfgdt, ADFGDT) (double *ID, STR_PSTR(dtype),
	int *ier STR_PLEN(dtype)) {
	int len_dtype = (int)STR_LEN(dtype);
	FMCALL(adfgdt2, ADFGDT2)(ID, STR_PTR(dtype), &len_dtype, ier);
}

void FMNAME(adfgdv, ADFGDV) (double *ID, int *dvals, int *ier) {
	FMCALL(adfgdv2, ADFGDV2)(ID, dvals, ier);
}

void FMNAME(adfges, ADFGES) (int *estate, int *ier) {
	FMCALL(adfges2, ADFGES2)(estate, ier);
}

void FMNAME(adfglb, ADFGLB) (double *ID, STR_PSTR(label),
	int *ier STR_PLEN(label)) {
	int len_label = (int)STR_LEN(label);
	FMCALL(adfglb2, ADFGLB2)(ID, STR_PTR(label), &len_label, ier);
}

void FMNAME(adfglkp, ADFGLKP) (double *ID, STR_PSTR(file), STR_PSTR(name),
	int *ier STR_PLEN(file) STR_PLEN(name)) {
	int len_file = (int)STR_LEN(file);
	int len_name = (int)STR_LEN(name);
	FMCALL(adfglk2, ADFGLK2)(ID, STR_PTR(file), &len_file, STR_PTR(name),
        	&len_name, ier);
}

void FMNAME(adfgnam, ADFGNAM) (double *ID, STR_PSTR(name),
	int *ier STR_PLEN(name)) {
	int len_name = (int)STR_LEN(name);
	FMCALL(adfgna2, ADFGNA2)(ID, STR_PTR(name), &len_name, ier);
}

void FMNAME(adfgnd, ADFGND) (double *ID, int *ndims, int *ier) {
	FMCALL(adfgnd2, ADFGND2)(ID, ndims, ier);
}

void FMNAME(adfgnid, ADFGNID) (double *PID, STR_PSTR(name), double *ID,
	int *ier STR_PLEN(name)) {
	int len_name = (int)STR_LEN(name);
	FMCALL(adfgni2, ADFGNI2)(PID, STR_PTR(name), &len_name, ID, ier);
}

void FMNAME(adfgrid, ADFGRID) (double *ID, double *RootID, int *ier) {
	FMCALL(adfgri2, ADFGRI2)(ID, RootID, ier);
}

void FMNAME(adfislk, ADFISLK) (double *ID, int *lplen, int *ier) {
	FMCALL(adfisl2, ADFISL2)(ID, lplen, ier);
}

void FMNAME(adflink, ADFLINK) (double *PID, STR_PSTR(name), STR_PSTR(file),
	STR_PSTR(nfile), double *ID, int *ier STR_PLEN(name)
	STR_PLEN(file) STR_PLEN(nfile)) {
	int len_name = (int)STR_LEN(name);
	int len_file = (int)STR_LEN(file);
	int len_nfile = (int)STR_LEN(nfile);
	FMCALL(adflin2, ADFLIN2)(PID, STR_PTR(name), STR_PTR(file),
        	STR_PTR(nfile), &len_name, &len_file, &len_nfile, ID, ier);
}

void FMNAME(adflver, ADFLVER) (STR_PSTR(version), int *ier STR_PLEN(version)) {
	int len_version = (int)STR_LEN(version);
	FMCALL(adflve2, ADFLVE2)(STR_PTR(version), &len_version, ier);
}

void FMNAME(adfmove, ADFMOVE) (double *PID, double *ID, double *NPID, int *ier) {
	FMCALL(adfmov2, ADFMOV2)(PID, ID, NPID, ier);
}

void FMNAME(adfncld, ADFNCLD) (double *ID, int *numcld, int *ier) {
	FMCALL(adfncl2, ADFNCL2)(ID, numcld, ier);
}

void FMNAME(adfpdim, ADFPDIM) (double *ID, STR_PSTR(dtype), int *dims,
	int *dvals, int *ier STR_PLEN(dtype)) {
	int len_dtype = (int)STR_LEN(dtype);
	FMCALL(adfpdi2, ADFPDI2)(ID, STR_PTR(dtype), &len_dtype, dims,
        	dvals, ier);
}

void FMNAME(adfpnam, ADFPNAM) (double *PID, double *ID, STR_PSTR(name),
	int *ier STR_PLEN(name)) {
	int len_name = (int)STR_LEN(name);
	FMCALL(adfpna2, ADFPNA2)(PID, ID, STR_PTR(name), &len_name, ier);
}

void FMNAME(adfses, ADFSES) (int *estate, int *ier) {
	FMCALL(adfses2, ADFSES2)(estate, ier);
}

void FMNAME(adfslb, ADFSLB) (double *ID, STR_PSTR(label),
	int *ier STR_PLEN(label)) {
	int len_label = (int)STR_LEN(label);
	FMCALL(adfslb2, ADFSLB2)(ID, STR_PTR(label), &len_label, ier);
}

/* added mainly to handle character output under windows */

void FMNAME(adfreadc, ADFREADC) (double *ID, int *s_start, int *s_end,
	int *s_stride, int *m_num_dims, int *m_dims, int *m_start, int *m_end,
	int *m_stride, STR_PSTR(data), int *ier STR_PLEN(data)) {
	FMCALL(adfread,ADFREAD)(ID, s_start, s_end, s_stride, m_num_dims,
                m_dims, m_start, m_end, m_stride, STR_PTR(data), ier);
}

void FMNAME(adfrallc, ADFRALLC) (double *ID, STR_PSTR(data),
	int *ier STR_PLEN(data)) {
        FMCALL(adfrall, ADFRALL)(ID, STR_PTR(data), ier);
}

void FMNAME(adfrblkc, ADFRBLKC) (double *ID, int *b_start, int *b_end,
	STR_PSTR(data), int *ier STR_PLEN(data)) {
        FMCALL(adfrblk, ADFRBLK)(ID, b_start, b_end, STR_PTR(data), ier);
}

void FMNAME(adfwritc, ADFWRITC) (double *ID, int *s_start, int *s_end,
	int *s_stride, int *m_num_dims, int *m_dims, int *m_start, int *m_end,
	int *m_stride, STR_PSTR(data), int *ier STR_PLEN(data)) {
	FMCALL(adfwrit,ADFWRIT)(ID, s_start, s_end, s_stride, m_num_dims,
                m_dims, m_start, m_end, m_stride, STR_PTR(data), ier);
}

void FMNAME(adfwallc, ADFWALLC) (double *ID, STR_PSTR(data),
	int *ier STR_PLEN(data)) {
        FMCALL(adfwall, ADFWALL)(ID, STR_PTR(data), ier);
}

void FMNAME(adfwblkc, ADFWBLKC) (double *ID, int *b_start, int *b_end,
	STR_PSTR(data), int *ier STR_PLEN(data)) {
        FMCALL(adfwblk, ADFWBLK)(ID, b_start, b_end, STR_PTR(data), ier);
}

