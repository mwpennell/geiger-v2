#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void pic_variance(void *, void *, void *, void *, void *, void *);
extern void r_mkn_core(void *, void *, void *, void *, void *, void *, void *, void *);
extern void initmod_mkn_pij(void *);
extern void derivs_mkn_pij(int *, double *, double *, double *, double *, int *);

/* .Call calls */
extern SEXP binary_edges(SEXP);
extern SEXP bm_direct(SEXP, SEXP);
extern SEXP cache_descendants(SEXP);
extern SEXP compile_ancestors(SEXP);
extern SEXP compile_descendants(SEXP);
extern SEXP get_descendants(SEXP);
extern SEXP open_subtree(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"pic_variance", (DL_FUNC) &pic_variance, 6},
    {"r_mkn_core",   (DL_FUNC) &r_mkn_core,   8},
    {"initmod_mkn_pij", (DL_FUNC) &initmod_mkn_pij, 1},
    {"derivs_mkn_pij", (DL_FUNC) &derivs_mkn_pij, 6},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"binary_edges",        (DL_FUNC) &binary_edges,        1},
    {"bm_direct",           (DL_FUNC) &bm_direct,           2},
    {"cache_descendants",   (DL_FUNC) &cache_descendants,   1},
    {"compile_ancestors",   (DL_FUNC) &compile_ancestors,   1},
    {"compile_descendants", (DL_FUNC) &compile_descendants, 1},
    {"get_descendants",     (DL_FUNC) &get_descendants,     1},
    {"open_subtree",        (DL_FUNC) &open_subtree,        2},
    {NULL, NULL, 0}
};

void R_init_geiger(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
