void do_gemm(double *x, int nrx, int ncx,
             double *y, int nry, int ncy,
             double *z);
SEXP getListElement(SEXP list, const char *str);
SEXP getListElementIfThere(SEXP list, const char *str);
