#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void basiserg(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void exber_maxwdth(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void extremfnpreg(void *, void *, void *, void *, void *, void *, void *);
extern void indextremw(void *, void *, void *, void *);
extern void multires(void *, void *, void *);
extern void p7fit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void pkcnt(void *, void *, void *, void *, void *, void *, void *);
extern void wsspoisschngd(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Fortran calls */
extern void F77_NAME(fcnst)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(npreg)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"basiserg",      (DL_FUNC) &basiserg,      10},
  {"exber_maxwdth", (DL_FUNC) &exber_maxwdth, 12},
  {"extremfnpreg",  (DL_FUNC) &extremfnpreg,   7},
  {"indextremw",    (DL_FUNC) &indextremw,     4},
  {"multires",      (DL_FUNC) &multires,       3},
  {"p7fit",         (DL_FUNC) &p7fit,         11},
  {"pkcnt",         (DL_FUNC) &pkcnt,          7},
  {"wsspoisschngd", (DL_FUNC) &wsspoisschngd, 13},
  {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
  {"fcnst", (DL_FUNC) &F77_NAME(fcnst), 10},
  {"npreg", (DL_FUNC) &F77_NAME(npreg), 24},
  {NULL, NULL, 0}
};

void R_init_diffractometry(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
