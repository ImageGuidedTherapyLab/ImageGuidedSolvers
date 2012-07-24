/*****************************************************************************/

/*
 * C/C++ compatability:
 */

/*
 * Use this in prototypes like this:  extern FORTRAN_FUNCTION void foo(...)
 *
 * At present, this is set up to tell a C++ compiler that  foo()  uses
 * a C-compatible calling convention.
 */
#ifdef __cplusplus
  #define FORTRAN_FUNCTION	"C"
#else
  #define FORTRAN_FUNCTION	/* empty */
#endif

/*****************************************************************************/

/* array subscripting offset, i.e. C "arr[k]" is Fortran "arr(k+?)" */
#define FORTRAN_INDEX_ORIGIN	1

/*****************************************************************************/


/*
 * Names of Fortran routines are often altered by the compiler/loader.  The
 * following macro should be used to call a Fortran routine from C code, i.e.
 *	call sgefa(...)			-- Fortran code
 *	FORTRAN_NAME(sgefa)(...);	-- C code to do the same thing
 *
 * Unfortunately, the "alterations" are generally at the token level, and this
 * can''t be done portably in pre-ANSI C.  In ANSI C, the preprocessor "token
 * splicing" facility is designed to handle just this sort of thing, but in
 * pre-ANSI C we have to use rather ugly system-dependent hacks of the sort
 * exemplified below.
 */

  /* C code should reference Fortran names in lower case */
  #define FORTRAN_NAME(x)       x ## _
