#include <R.h>
#include <Rinternals.h>

SEXP  mmap_dataSymbol;

extern SEXP mmapReplaceReal(SEXP index, SEXP value, SEXP mmap_obj) {
	int i, LEN = length(index);
	double *data;
	double *index_p = REAL(index);
	double *real_value = REAL(value);

	mmap_dataSymbol = install("data");
	data = R_ExternalPtrAddr(findVar(mmap_dataSymbol,mmap_obj));

	for(i=0; i < LEN; i++) {

		Rprintf("i %d p %f p2 %d v %f o %f\n", i, index_p[i], ((long)index_p[i]-1),
				real_value[i], data[((long)index_p[i]-1)]);

		memcpy(
				&(data[((long)index_p[i]-1)]),
				&(real_value[i]),
				sizeof(double));
	}
	return R_NilValue;
}


static const R_CMethodDef CEntries[] = {
		{NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
		{"mmapReplaceReal",       (DL_FUNC) &mmapReplaceReal,       3},
		{NULL, NULL, 0}
};

void R_init_localEM(DllInfo *dll)
{
	R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
}
