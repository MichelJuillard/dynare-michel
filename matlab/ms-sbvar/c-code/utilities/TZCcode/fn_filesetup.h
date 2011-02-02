#ifndef __FN_FILESETUP_H__
#define __FN_FILESETUP_H__
   #include <string.h>
/*     //#include <malloc.h>                  // For malloc, calloc, etc.   ansi-c*/

   #include "tzmatlab.h"
   #include "mathlib.h"  /*  Used for tz_TransposeRegular().   ansi-c*/

   int fn_ParseCommandLine(int n_arg, char **args, char ch);
   char *fn_ParseCommandLine_String(int n_arg, char **args, char ch, char *default_return);
   int fn_ParseCommandLine_Integer(int n_arg, char **args, char ch, int default_return);
   int fn_SetFilePosition(FILE *f, const char *id);

   int fn_ReadVector_int(FILE *fptr, int *x_v, const int d_x_v);
   int fn_ReadVector_lf(FILE *fptr, double *x_v, const int d_x_v);
   int fn_ReadMatrix_int(FILE *fptr, int *x_m, const int r_x_m, const int c_x_m);
   int fn_ReadMatrix_lf(FILE *fptr, double *x_m, const int r_x_m, const int c_x_m);

   int ReadNullTerminatedString(FILE *fptr, TScvector *x_cv);
   int ReadVector_int(FILE *fptr, TSivector *x_iv);
   int ReadVector_lf(FILE *fptr, TSdvector *x_dv);
   int ReadMatrix_int(FILE *fptr, TSimatrix *X_im);
   int ReadMatrix_lf(FILE *fptr, TSdmatrix *x_dm);
   int ReadCell_lf(FILE *fptr, TSdcell *x_dc);
   int ReadCellvec_lf(FILE *fptr, TSdcellvec *x_dcv);

   void fn_WriteMatrix_f(FILE *fprt_debug, const double *x_m, const int r_x_m, const int c_x_m);
   void fn_WriteMatrix_lf(FILE *fprt_debug, const double *x_m, const int r_x_m, const int c_x_m);
   void fn_WriteMatrix_int(FILE *fprt_debug, const int *x_m, const int r_x_m, const int c_x_m);
   void fn_WriteVector_f(FILE *fprt_debug, const double *x_v, const int d_x_v);

   void WriteMatrix_f(FILE *fprt_debug, const TSdmatrix *x_dm);
   void WriteMatrix_lf(FILE *fprt_debug, const TSdmatrix *x_dm);
   void WriteMatrix(FILE *fprt_debug, const TSdmatrix *x_dm, const char *format);
   void WriteMatrixTranspose(FILE *fptr_debug, const TSdmatrix *x_dm, const char *format);
   void WriteCell_lf(FILE *fprt_debug, const TSdcell *x_dc);
   void WriteCell_f(FILE *fprt_debug, const TSdcell *x_dc);
   void WriteCell(FILE *fprt_debug, const TSdcell *x_dc, const char *format);
   void WriteCellTranspose(FILE *fptr_debug, const TSdcell *x_dc, const char *format);
   void WriteCellvec_lf(FILE *fprt_debug, const TSdcellvec *x_dcv);
   void WriteCellvec_f(FILE *fprt_debug, const TSdcellvec *x_dcv);
   void WriteCellvec(FILE *fptr_debug, const TSdcellvec *x_dcv, const char *format);
   void WriteFourth_f(FILE *fptr_debug, const TSdfourth *x_d4);
   void WriteFourth(FILE *fptr_debug, const TSdfourth *x_d4, const char *format);
   void WriteVector_f(FILE *fprt_debug, const TSdvector *x_dv);
   void WriteVector_lf(FILE *fprt_debug, const TSdvector *x_dv);
   void WriteVector(FILE *fprt_debug, const TSdvector *x_dv, const char *format);
   void WriteVector_column(FILE *fptr_debug, const TSdvector *x_dv, const char *format);
   void WriteCellvec_int(FILE *fptr_debug, const TSicellvec *x_icv);
   void WriteMatrix_int(FILE *fprt_debug, const TSimatrix *x_im);
   void WriteVector_int(FILE *fprt_debug, const TSivector *x_iv);


   void PrintVector_int(const TSivector *x_iv);
   void PrintVector(const TSdvector *x_dv, const char *format);
   void PrintVector_f(const TSdvector *x_dv);
   void PrintVector_dz(const TSdzvector *x_dzv);
   void PrintMatrix_int(const TSimatrix *X_im);
   void PrintMatrix_f(const TSdmatrix *x_dm);
   void PrintMatrix(const TSdmatrix *x_dm, const char *format);
   void PrintMatrix_dz(const TSdzmatrix *x_dzm);
   void PrintCellvec_f(const TSdcellvec *x_dcv);
   void PrintCell_f(const TSdcell *x_dc);
   void PrintCell(const TSdcell *x_dc, const char *format);
   void PrintFourthvec_f(TSdfourthvec *x_d4v);


   void ReprintInputData(FILE *fptr_in, FILE *fptr_out);
#endif
