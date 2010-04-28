/***********
 * Reads the input file name and output file names specified by the user from the command line with automatic default to
 *   both input an output files.
***********/

#include "fn_filesetup.h"

#include "modify_for_mex.h"

//-----------------
// For command line.
// Finds /ch in the command line.  If found, returns the args location
//   indexed by int and zero otherwise.
//-----------------
int fn_ParseCommandLine(int n_arg, char **args, char ch) {
   int i;
   for (i=1; i<n_arg; i++)
      if ((args[i][0] == '/') && (args[i][1] == ch)) return i;
   return 0;
}


//-----------------
// For command line.
// Finds /ch in the command line.  If found returns a pointer
//   to the string trailing /ch.  If /ch is not found or there is
//   no trailing string or the trailing string is another argument,
//   then default_return is returned.  No memory is allocated and
//   the calling routine should not free the returned pointer.
//-----------------
char *fn_ParseCommandLine_String(int n_arg, char **args, char ch, char *default_return) {
   int i=fn_ParseCommandLine(n_arg,args,ch);
   if (i > 0)
      if (strlen(args[i]) > 2) return args[i]+2;
           // In case the user forgot typing a space between /ch and string following it, still returns a pointer to the string folloing /ch.
      else if ((i+1 < n_arg) && (args[i+1][0] != '/')) return args[i+1];
           // Returns a pointer to the string that does NOT begin with / and there is a whitespace between /ch and the string.
   return default_return;
}


//-----------------
// For command line.
// Finds /ch in the command line.  If found returns the integer
//   value of the string trailing /ch (e.g, the integer value is
//   sample size or normalization index.  If /ch is not found or there
//   is no trailing string or the trailing string is another argument,
//   then the default_return value is returned.
//-----------------
int fn_ParseCommandLine_Integer(int n_arg, char **args, char ch, int default_return) {
   char *str=fn_ParseCommandLine_String(n_arg,args,ch,(char*)NULL);
   return str ? atoi(str) : default_return;
}


//-----------------
// Finds proper location in the input data file.
// Returns 1 if the NUL-terminated string id is found
//   in the file and 0 otherwise.  The file pointer is set
//   to the line immediately after the line containing id.
//   If the string id has a length (including the new line
//   character \n) more than 1023, it will be cut off at 1023.
//-----------------
int fn_SetFilePosition(FILE *f, const char *id) {
   // As an output, the file pointer f will be reset to the beginning of the line next to the line headed by the string id.
   char buffer[1024];
   size_t n=strlen(id);
   int ch;

   if ( !f )  fn_DisplayError(".../fn_filesetup.c/fn_SetFilePosition():  the file, *f, must be created (opened)");
   if (n>1023)  n=1023;
   rewind(f);                             // Reset a file poiniter to the beginning of the file.  There may be more efficient ways but this is good enough as long as the file is not too long.
   while (fgets(buffer,1024,f)) {        // Reads a line at a time in the file f (including \n and a NUL byte) until it matches id.  fgets returns the pointer to the buffer and is often only used to check for EOF.
      if (buffer[strlen(buffer)-1] != '\n')  // -1 because the first element of the buffer is indexed by buffer[0].
         // If the end of the buffer (excluding the NUL byte) encounters no new line, f points to the next character after
         //   the end of the buffer on the SAME line (i.e., f does not point to the begining of the new line at this point).
         //   The following do loop will take f to point to the beginning of the new line.
         do ch=fgetc(f);  // Gets one character at a time until it reachs the end of the current '\n' or the end of the file EOF.
         while ( (ch != '\n') && (ch != EOF) );
      if (!memcmp(buffer,id,n)) return 1;  // The match is found.
   }
   return 0;                            // No match is found.
}


//-----------------
// Reads a string from the input data file with the NULL-terminated
//   character but without the new line character.
// Returns 1 if the vector of characters is all read without
//   errors and 0 otherwise.  The file pointer is then moved
//   to point to the next non-whitespace character after these
//   characters.
//-----------------
int ReadNullTerminatedString(FILE *fptr, TScvector *x_cv)
{
   //x_cv will have a string without the new line character and with the NULL character.
   //It is the user's responsiblity to ensure the string x_cv has an enough length to use fgets().
   //  If not, it stops after x_cv->n-1 characters have been stored in x_cv->v and a NULL byte is appended to make it a string.
   //  If yet, reading stops after a newline character is read and stored in x_cv->v and a NULL byte is then appended.
   int _n;
   char *cv;
   if (!fptr || !x_cv)  fn_DisplayError(".../fn_filesetup.c/ReadNullTerminatedString(): File or input string must be created (memory-allocated)");
   _n = x_cv->n;
   cv = x_cv->v;
   if ( !fgets(cv, _n, fptr) ) return 0;
   cv[strlen(cv)-1] = '\0';  //Removes the new line character and replaces it with the NULL character.
      //The string length (size_t type) strlen(cv) does NOT count the NULL byte at the end, but it counts the new line character.
   return 1;
}


//-----------------
// Reads a vector of integers from the input data file.
// Returns 1 if the vector of integers is all read without
//   errors and 0 otherwise.  The file pointer is then moved
//   to point to the next non-whitespace character after these
//   integers.
//-----------------
int fn_ReadVector_int(FILE *fptr, int *x_v, const int d_x_v) {
   int ki;
   for (ki=0; ki<d_x_v; ki++)
      if ( fscanf(fptr, " %d ", &x_v[ki]) !=1 ) return 0;
   return 1;
}
int ReadVector_int(FILE *fptr, TSivector *x_iv) {
   int ki, _n,
       *v;
   if (!fptr || !x_iv) fn_DisplayError(".../fn_filesetup.c/ReadVector_int(): File or input matrix must be created (memory-allocated)");
   _n = x_iv->n;
   v = x_iv->v;
   for (ki=0; ki<_n; ki++)
      if ( fscanf(fptr, " %d ", &v[ki]) != 1 ) return 0;
   return 1;
}


//-----------------
// Reads a vector of doubles from the input data file.
// Returns 1 if the vector of doubles is all read without
//   errors and 0 otherwise.  The file pointer is then moved
//   to point to the next non-whitespace character after these
//   doubles.
//-----------------
int fn_ReadVector_lf(FILE *fptr, double *x_v, const int d_x_v) {
   int ki;
   for (ki=0; ki<d_x_v; ki++)
      if ( fscanf(fptr, " %lf ", &x_v[ki]) !=1 ) return 0;
   return 1;
}
int ReadVector_lf(FILE *fptr, TSdvector *x_dv) {
   int ki, _n;
   double *v;
   if (!fptr || !x_dv) fn_DisplayError(".../fn_filesetup.c/ReadVector_lf(): File or input matrix must be created (memory-allocated)");
   _n = x_dv->n;
   v = x_dv->v;
   for (ki=0; ki<_n; ki++)
      if ( fscanf(fptr, " %lf ", &v[ki]) != 1 ) return 0;

   x_dv->flag = V_DEF;
   return 1;
}


//-----------------
// Reads a column-major matrix of integers from the input data file.
// Returns 1 if the matrix of integers is all read without
//   errors and 0 otherwise.  The file pointer is then moved
//   to point to the next non-whitespace character after these
//   integers.
//-----------------
int fn_ReadMatrix_int(FILE *fptr, int *x_m, const int r_x_m, const int c_x_m) {
   int ki, kj;

   for (ki=0; ki<r_x_m; ki++)
      for (kj=0; kj<c_x_m; kj++)
         if ( fscanf(fptr, " %d ", &x_m[kj*r_x_m+ki]) !=1 ) return 0;
   return 1;
}
int ReadMatrix_int(FILE *fptr, TSimatrix *X_im)
{
   int ki, kj;
   int nrows, ncols;
   if (!fptr || !X_im) fn_DisplayError(".../fn_filesetup.c/ReadMatrix_int(): File or input matrix must be created (memory-allocated)");

   nrows = X_im->nrows;
   ncols = X_im->ncols;
   for (ki=0; ki<nrows; ki++)
      for (kj=0; kj<ncols; kj++)
         if ( fscanf(fptr, " %d ", (X_im->M+mos(ki,kj,nrows))) !=1 ) return 0;
   return 1;
}


//-----------------
// Reads a column-major matrix of doubles from the input data file.
// Returns 1 if the matrix of doubles is all read without
//   errors and 0 otherwise.  The file pointer is then moved
//   to point to the next non-whitespace character after these
//   doubles.
//-----------------
int fn_ReadMatrix_lf(FILE *fptr, double *x_m, const int r_x_m, const int c_x_m) {
   int ki, kj;
   for (ki=0; ki<r_x_m; ki++)
      for (kj=0; kj<c_x_m; kj++)
         if ( fscanf(fptr, " %lf ", &x_m[kj*r_x_m+ki]) !=1 ) return 0;
   return 1;
}
int ReadMatrix_lf(FILE *fptr, TSdmatrix *x_dm) {
   //Outputs:
   //  x_dm (whose memory is already allocated): To be filled with the numbers from the file fptr.
   int ki, kj, nrows, ncols;
   double *M;
   if (!fptr || !x_dm) fn_DisplayError(".../fn_filesetup.c/ReadMatrix_lf(): File or input matrix must be created (memory-allocated)");
   nrows = x_dm->nrows;
   ncols = x_dm->ncols;
   M = x_dm->M;
   for (ki=0; ki<nrows; ki++)
      for (kj=0; kj<ncols; kj++)
         if ( fscanf(fptr, " %lf ", &M[mos(ki,kj,nrows)]) !=1 ) return 0;

   x_dm->flag = M_GE;
   return 1;
}



//-----------------
// Reads a column-major cell of double vectors from the input data file.
// Returns 1 if all data are read without errors and 0 otherwise.
// The file pointer is then moved to point to the next non-whitespace character
//   after these doubles.
//-----------------
int ReadCellvec_lf(FILE *fptr, TSdcellvec *x_dcv) {
   //Outputs:
   //  x_dcv (whose memory is already allocated): To be filled with the numbers from the file fptr.
   int ci, kj, _n, ncells;
   double *v;
   if (!fptr || !x_dcv) fn_DisplayError(".../fn_filesetup.c/ReadCellvec_lf(): File or input cell must be created (memory-allocated)");
   ncells = x_dcv->ncells;
   for (ci=0; ci<ncells; ci++) {
      _n = x_dcv->C[ci]->n;
      v = x_dcv->C[ci]->v;
      for (kj=0; kj<_n; kj++)
         if ( fscanf(fptr, " %lf ", &v[kj]) != 1 ) return 0;
   }
   return 1;
}




//-----------------
// Reads a column-major cell of double matrices from the input data file.
// Returns 1 if all data are  read without errors and 0 otherwise.
// The file pointer is then moved to point to the next non-whitespace character
//   after these doubles.
//-----------------
int ReadCell_lf(FILE *fptr, TSdcell *x_dc) {
   //Outputs:
   //  x_dc (whose memory is already allocated): To be filled with the numbers from the file fptr.
   int ci, ki, kj, nrows, ncols, ncells;
   double *M;
   if (!fptr || !x_dc) fn_DisplayError(".../fn_filesetup.c/ReadCell_lf(): File or input cell must be created (memory-allocated)");
   ncells = x_dc->ncells;
   for (ci=0; ci<ncells; ci++) {
      nrows = x_dc->C[ci]->nrows;
      ncols = x_dc->C[ci]->ncols;
      M = x_dc->C[ci]->M;
      for (ki=0; ki<nrows; ki++)
         for (kj=0; kj<ncols; kj++)
            if ( fscanf(fptr, " %lf ", &M[mos(ki,kj,nrows)]) != 1 ) return 0;
   }
   return 1;
}



//-----------------
// Write a column-major matrix of floats to the output file.
// The file pointer is then moved to point to the next
//   non-whitespace character after these doubles.
//-----------------
void fn_WriteMatrix_f(FILE *fptr_debug, const double *x_m, const int r_x_m, const int c_x_m) {
   int _i, _j;

   for (_i=0; _i<r_x_m; _i++) {
      for (_j=0; _j<c_x_m; _j++) {
         fprintf(fptr_debug, " %f ", x_m[_j*r_x_m + _i]);
         if (_j==c_x_m-1) fprintf(fptr_debug, "\n");
      }
      if (_i==r_x_m-1) fprintf(fptr_debug, "\n\n");
   }
}
void WriteMatrix_f(FILE *fptr_debug, const TSdmatrix *x_dm) {
   int _i, _j;
   if (!fptr_debug || !x_dm) fn_DisplayError(".../fn_filesetup.c/WriteMatrix_f(): File or input matrix cannot be NULL (must be created)");
   for (_i=0; _i<x_dm->nrows; _i++) {
      for (_j=0; _j<x_dm->ncols; _j++) {
         fprintf(fptr_debug, " %10.5f ", x_dm->M[_j*x_dm->nrows + _i]);
         if (_j==x_dm->ncols-1) fprintf(fptr_debug, "\n");
      }
      if (_i==x_dm->nrows-1) fprintf(fptr_debug, "\n\n");
   }
}


//-----------------
// Write a column-major matrix of doubles to the output file.
// The file pointer is then moved to point to the next
//   non-whitespace character after these doubles.
//-----------------
void fn_WriteMatrix_lf(FILE *fptr_debug, const double *x_m, const int r_x_m, const int c_x_m) {
   int _i, _j;
   for (_i=0; _i<r_x_m; _i++) {
      for (_j=0; _j<c_x_m; _j++) {
         fprintf(fptr_debug, " %.16e ", x_m[_j*r_x_m + _i]);
         if (_j==c_x_m-1) fprintf(fptr_debug, "\n");
      }
      if (_i==r_x_m-1) fprintf(fptr_debug, "\n\n");
   }
}
void WriteMatrix_lf(FILE *fptr_debug, const TSdmatrix *x_dm) {
   int _i, _j;
   if (!fptr_debug || !x_dm) fn_DisplayError(".../fn_filesetup.c/WriteMatrix_lf(): File or input matrix cannot be NULL (must be created)");
   for (_i=0; _i<x_dm->nrows; _i++) {
      for (_j=0; _j<x_dm->ncols; _j++) {
         fprintf(fptr_debug, " %.16e ", x_dm->M[_j*x_dm->nrows + _i]);
         if (_j==x_dm->ncols-1) fprintf(fptr_debug, "\n");
      }
      if (_i==x_dm->nrows-1) fprintf(fptr_debug, "\n\n");
   }
}
void WriteMatrix(FILE *fptr_debug, const TSdmatrix *x_dm, const char *format) {
   int _i, _j, nrows, ncols;
   double *M;
   if (!fptr_debug || !x_dm) fn_DisplayError(".../fn_filesetup.c/WriteMatrix(): File or input matrix cannot be NULL (must be created)");
   nrows = x_dm->nrows;
   ncols = x_dm->ncols;
   M = x_dm->M;
   if (!format)   format=" %10.5f ";   //Default format.
   for (_i=0; _i<nrows; _i++)
      for (_j=0; _j<ncols; _j++) {
         fprintf(fptr_debug, format, M[_j*x_dm->nrows + _i]);
         if (_j==ncols-1) fprintf(fptr_debug, "\n");
      }
   //fprintf(fptr_debug, "\n");
}
//+
void WriteMatrixTranspose(FILE *fptr_debug, const TSdmatrix *x_dm, const char *format)
{
   int _i, _j, nrows, ncols;
   double *M;
   //===
   TSdmatrix *Xtran_dm = NULL;

   if (!fptr_debug || !x_dm) fn_DisplayError(".../fn_filesetup.c/WriteMatrixTranspose(): File or input matrix cannot be NULL (must be created)");

   Xtran_dm = tz_TransposeRegular((TSdmatrix *)NULL, x_dm);

   nrows = Xtran_dm->nrows;
   ncols = Xtran_dm->ncols;
   M = Xtran_dm->M;
   if (!format)   format=" %10.5f ";   //Default format.
   for (_i=0; _i<nrows; _i++)
      for (_j=0; _j<ncols; _j++) {
         fprintf(fptr_debug, format, M[_j*Xtran_dm->nrows + _i]);
         if (_j==ncols-1) fprintf(fptr_debug, "\n");
      }
   //fprintf(fptr_debug, "\n");

   //===
   DestroyMatrix_lf(Xtran_dm);
}


//-----------------
// Write cells of column-major double matrices to the output file.
// The file pointer is then moved to point to the next
//   non-whitespace character after these doubles.
//-----------------
void WriteCell_lf(FILE *fptr_debug, const TSdcell *x_dc) {
   int _i, _n;
   if (!fptr_debug || !x_dc) fn_DisplayError(".../fn_filesetup.c/WriteCell_lf(): File or input cell cannot be NULL (must be created)");
   _n = x_dc->ncells;
   for (_i=0; _i<_n; _i++) {
      fprintf(fptr_debug, "Cell %d\n", _i);
      WriteMatrix_lf(fptr_debug, x_dc->C[_i]);
   }
}
void WriteCell_f(FILE *fptr_debug, const TSdcell *x_dc) {
   int _i, _n;
   if (!fptr_debug || !x_dc) fn_DisplayError(".../fn_filesetup.c/WriteCell_f(): File or input cell cannot be NULL (must be created)");
   _n = x_dc->ncells;
   for (_i=0; _i<_n; _i++) {
      fprintf(fptr_debug, "Cell %d\n", _i);
      WriteMatrix_f(fptr_debug, x_dc->C[_i]);
   }
}
void WriteCell(FILE *fptr_debug, const TSdcell *x_dc, const char *format) {
   int _i, _n;
   if (!fptr_debug || !x_dc) fn_DisplayError(".../fn_filesetup.c/WriteCell(): File or input cell cannot be NULL (must be created)");
   _n = x_dc->ncells;
   for (_i=0; _i<_n; _i++)
   {
      WriteMatrix(fptr_debug, x_dc->C[_i], format);
      fprintf(fptr_debug, "\n");
   }
}
//+
void WriteCellTranspose(FILE *fptr_debug, const TSdcell *x_dc, const char *format)
{
   int _i, _n;
   if (!fptr_debug || !x_dc) fn_DisplayError(".../fn_filesetup.c/WriteCell(): File or input cell cannot be NULL (must be created)");
   _n = x_dc->ncells;
   for (_i=0; _i<_n; _i++)
   {
      WriteMatrixTranspose(fptr_debug, x_dc->C[_i], format);
      fprintf(fptr_debug, "\n");
   }
}


//-----------------
// Write cells of vectors to the output file.
// The file pointer is then moved to point to the next
//   non-whitespace character after these doubles.
//-----------------
void WriteCellvec_lf(FILE *fptr_debug, const TSdcellvec *x_dcv) {
   int _i;
   if (!fptr_debug || !x_dcv) fn_DisplayError(".../fn_filesetup.c/WriteCellvec_lf(): File or input cell cannot be NULL (must be created)");
   for (_i=0; _i<x_dcv->ncells; _i++) {
      fprintf(fptr_debug, "Cell %d\n", _i);
      WriteVector_lf(fptr_debug, x_dcv->C[_i]);
   }
}
void WriteCellvec_f(FILE *fptr_debug, const TSdcellvec *x_dcv) {
   int _i;
   if (!fptr_debug || !x_dcv) fn_DisplayError(".../fn_filesetup.c/WriteCellvec_lf(): File or input cell cannot be NULL (must be created)");
   for (_i=0; _i<x_dcv->ncells; _i++) {
      fprintf(fptr_debug, "Cell %d\n", _i);
      WriteVector_f(fptr_debug, x_dcv->C[_i]);
   }
}
void WriteCellvec(FILE *fptr_debug, const TSdcellvec *x_dcv, const char *format) {
   int _i, _n;
   if (!fptr_debug || !x_dcv) fn_DisplayError(".../fn_filesetup.c/WriteCellvec(): File or input cell cannot be NULL (must be created)");
   _n = x_dcv->ncells;
   for (_i=0; _i<_n; _i++)  WriteVector(fptr_debug, x_dcv->C[_i], format);
}
void WriteCellvec_int(FILE *fptr_debug, const TSicellvec *x_icv)
{
   int _i, _n;
   if (!fptr_debug || !x_icv) fn_DisplayError(".../fn_filesetup.c/WriteCellvec_int(): File or input cell cannot be NULL (must be created)");
   _n = x_icv->ncells;
   for (_i=0; _i<_n; _i++)  WriteVector_int(fptr_debug, x_icv->C[_i]);
}



//-----------------
// Write fourths of column-major double matrices to an output file.
// The file pointer is then moved to point to the next
//   non-whitespace character after these doubles.
//-----------------
void WriteFourth_f(FILE *fptr_debug, const TSdfourth *x_d4) {
   int _j, _i, _m, _n;
   if (!fptr_debug || !x_d4) fn_DisplayError(".../fn_filesetup.c/WriteFourth_f(): File or input fourth cannot be NULL (must be created)");
   _m = x_d4->ndims;
   for (_j=0; _j<_m; _j++) {
      _n = x_d4->F[_j]->ncells;
      fprintf(fptr_debug, "Fourth %d\n", _j);
      for (_i=0; _i<_n; _i++) {
         fprintf(fptr_debug, "Cell %d\n", _i);
         WriteMatrix_f(fptr_debug, x_d4->F[_j]->C[_i]);
      }
   }
}
void WriteFourth(FILE *fptr_debug, const TSdfourth *x_d4, const char *format) {
   int _j, _i, _m, _n;
   if (!fptr_debug || !x_d4) fn_DisplayError(".../fn_filesetup.c/WriteFourth_f(): File or input fourth cannot be NULL (must be created)");
   _m = x_d4->ndims;
   for (_j=0; _j<_m; _j++) {
      _n = x_d4->F[_j]->ncells;
      for (_i=0; _i<_n; _i++) {
         WriteMatrix(fptr_debug, x_d4->F[_j]->C[_i], format);
      }
   }
}


//-----------------
// Write a column-major matrix of ints to the output file.
// The file pointer is then moved to point to the next
//   non-whitespace character after these doubles.
//-----------------
void fn_WriteMatrix_int(FILE *fptr_debug, const int *x_m, const int r_x_m, const int c_x_m) {
   int _i, _j;
   for (_i=0; _i<r_x_m; _i++) {
      for (_j=0; _j<c_x_m; _j++) {
         fprintf(fptr_debug, " %d ", x_m[_j*r_x_m + _i]);
         if (_j==c_x_m-1) fprintf(fptr_debug, "\n");
      }
      if (_i==r_x_m-1) fprintf(fptr_debug, "\n\n");
   }
}
void WriteMatrix_int(FILE *fptr_debug, const TSimatrix *x_im) {
   int _i, _j;
   if (!fptr_debug || !x_im) fn_DisplayError(".../fn_filesetup.c/WriteMatrix_int(): File or input matrix cannot be NULL (must be created)");
   for (_i=0; _i<x_im->nrows; _i++) {
      for (_j=0; _j<x_im->ncols; _j++) {
         fprintf(fptr_debug, " %d ", x_im->M[_j*x_im->nrows + _i]);
         if (_j==x_im->ncols-1) fprintf(fptr_debug, "\n");
      }
      if (_i==x_im->nrows-1) fprintf(fptr_debug, "\n\n");
   }
}


//-----------------
// Write a vector of doubles to the output file.
// The file pointer is then moved to point to the next
//   non-whitespace character after these doubles.
//-----------------
void fn_WriteVector_lf(FILE *fptr_debug, const double *x_v, const int d_x_v) {
   int _i;
   for (_i=0; _i<d_x_v; _i++) {
      fprintf(fptr_debug, " %20.16f ", x_v[_i]);
      if (_i==d_x_v-1) fprintf(fptr_debug, "\n\n");
   }
}
void WriteVector_lf(FILE *fptr_debug, const TSdvector *x_dv) {
   int _i;
   for (_i=0; _i<x_dv->n; _i++) {
      fprintf(fptr_debug, " %20.16f ", x_dv->v[_i]);
      if (_i==x_dv->n-1) fprintf(fptr_debug, "\n\n");
   }
}
void WriteVector(FILE *fptr_debug, const TSdvector *x_dv, const char *format) {
   int _i, _n;
   double *v;
   if ( !fptr_debug || !x_dv ) fn_DisplayError(".../fn_filesetup.c/WriteVector(): File or input vector cannot be NULL (must be created)");
   _n = x_dv->n;
   v = x_dv->v;
   if (!format)  format=" %10.5f ";   //Default format.
   for (_i=0; _i<_n; _i++)  fprintf(fptr_debug, format, v[_i]);
   fprintf(fptr_debug, "\n");
}
void WriteVector_column(FILE *fptr_debug, const TSdvector *x_dv, const char *format)
{
   int _i, _n;
   double *v;
   if ( !fptr_debug || !x_dv ) fn_DisplayError(".../fn_filesetup.c/WriteVector_column(): File or input vector cannot be NULL (must be created)");
   _n = x_dv->n;
   v = x_dv->v;
   if (!format)  format=" %10.5f ";   //Default format.
   for (_i=0; _i<_n; _i++)
   {
      fprintf(fptr_debug, format, v[_i]);
      fprintf(fptr_debug, "\n");
   }
}


//-----------------
// Write a vector of floats to the output file.
// The file pointer is then moved to point to the next
//   non-whitespace character after these doubles.
//-----------------
void fn_WriteVector_f(FILE *fptr_debug, const double *x_v, const int d_x_v) {
   int _i;
   for (_i=0; _i<d_x_v; _i++)  fprintf(fptr_debug, " %f ", x_v[_i]);
   fprintf(fptr_debug, "\n");
}
void WriteVector_f(FILE *fptr_debug, const TSdvector *x_dv) {
   int _i;
   if (!fptr_debug || !x_dv) fn_DisplayError(".../fn_filesetup.c/WriteVector_f(): File or input vector cannot be NULL (must be created)");
   for (_i=0; _i<x_dv->n; _i++)  fprintf(fptr_debug, " %10.5f ", x_dv->v[_i]);
   fprintf(fptr_debug, "\n");
}



//-----------------
// Write a vector of integers to the output file.
// The file pointer is then moved to point to the next
//   non-whitespace character after these doubles.
//-----------------
void WriteVector_int(FILE *fptr_debug, const TSivector *x_iv)
{
   int _i;
   if (!fptr_debug || !x_iv) fn_DisplayError(".../fn_filesetup.c/WriteVector_int(): File or input vector cannot be NULL (must be created)");
   for (_i=0; _i<x_iv->n; _i++) {
      fprintf(fptr_debug, " %d ", x_iv->v[_i]);
      if (_i==x_iv->n-1) fprintf(fptr_debug, "\n\n");
   }
}



void PrintVector_int(const TSivector *x_iv)
{
   int _i, _n;

   if (!x_iv) fn_DisplayError(".../fn_filesetup.c/PrintVector_int(): Input vector must be created (memory-allocated)");
   _n = x_iv->n;
   // printf("\nVector:\n");
   for (_i=0; _i<_n; _i++) {
      printf("v[%d]=%d\n", _i, x_iv->v[_i]);
   }
}

//-----------------
// Print a vector of doubles to the screen.
//-----------------
void PrintVector(const TSdvector *x_dv, const char *format)
{
   int _i, _n;

   if (!x_dv) fn_DisplayError(".../fn_filesetup.c/PrintVector(): Input vector must be created (memory-allocated)");
   _n = x_dv->n;
   // printf("\n\nVector:\n");
   for (_i=0; _i<_n; _i++) {
      printf(format, x_dv->v[_i]);
   }
}
//+
void PrintVector_f(const TSdvector *x_dv)
{
   int _i, _n;

   if (!x_dv) fn_DisplayError(".../fn_filesetup.c/PrintVector_f(): Input vector must be created (memory-allocated)");
   _n = x_dv->n;
   // printf("\n\nVector:\n");
   for (_i=0; _i<_n; _i++) {
      printf("v[%d]=%6.4f\n", _i, x_dv->v[_i]);
   }
}

void PrintVector_dz(const TSdzvector *x_dzv)
{
   int _i;

   if (!x_dzv) fn_DisplayError(".../fn_filesetup.c/PrintVector_dz(): Input complex vector must be created (memory-allocated)");

   printf("\n\nComplex vector:\n");
   for (_i=0; _i<x_dzv->real->n; _i++) {
      printf("vreal[%d]=%6.4f;  vimag[%d]=%6.4f\n", _i, x_dzv->real->v[_i], _i, x_dzv->imag->v[_i]);
   }
}


void PrintMatrix_int(const TSimatrix *X_im)
{
   int _i, _j, nrows, ncols;
   int *M=X_im->M;

   if (!X_im) fn_DisplayError(".../fn_filesetup.c/PrintMatrix_int(): Input matrix must be created (memory-allocated)");
   else {
       nrows = X_im->nrows;
       ncols = X_im->ncols;
       M = X_im->M;
   }

   printf("\n\nMatrix:\n");
   for (_i=0; _i<nrows; _i++) {
      for (_j=0; _j<ncols; _j++) {
         printf(" %d ", M[_j*nrows + _i]);
         if (_j==ncols-1) printf("\n");
      }
      if (_i==nrows-1) printf("\n");
   }
}

void PrintMatrix_f(const TSdmatrix *x_dm)
{
   int _i, _j, nrows, ncols;
   double *M=x_dm->M;

   if (!x_dm) fn_DisplayError(".../fn_filesetup.c/PrintMatrix_f(): Input matrix must be created (memory-allocated)");
   else {
       nrows = x_dm->nrows;
       ncols = x_dm->ncols;
       M = x_dm->M;
   }

   printf("\n\nMatrix:\n");
   for (_i=0; _i<nrows; _i++) {
      for (_j=0; _j<ncols; _j++) {
         printf(" %6.4f ", M[_j*nrows + _i]);
         if (_j==ncols-1) printf("\n");
      }
      if (_i==nrows-1) printf("\n");
   }
}

void PrintMatrix(const TSdmatrix *x_dm, const char *format)
{
   int _i, _j, nrows, ncols;
   double *M=x_dm->M;

   if (!x_dm) fn_DisplayError(".../fn_filesetup.c/PrintMatrix_f(): Input matrix must be created (memory-allocated)");
   else {
       nrows = x_dm->nrows;
       ncols = x_dm->ncols;
       M = x_dm->M;
   }

   printf("\n\nMatrix:\n");
   if (!format)  format=" %10.5f ";   //Default format.
   for (_i=0; _i<nrows; _i++) {
      for (_j=0; _j<ncols; _j++) {
         printf(format, M[_j*nrows + _i]);
         if (_j==ncols-1) printf("\n");
      }
      if (_i==nrows-1) printf("\n");
   }
}

void PrintMatrix_dz(const TSdzmatrix *x_dzm) {
   int _i, _j, nrows, ncols;
   double *Mr=NULL,
          *Mi=NULL;

   if (!x_dzm) fn_DisplayError(".../fn_filesetup.c/PrintMatrix_dz(): Input complex matrix must be created (memory-allocated)");
   else {
       nrows = x_dzm->real->nrows;
       ncols = x_dzm->real->ncols;
       Mr = x_dzm->real->M,
       Mi = x_dzm->imag->M;
   }

   printf("\n\nReal part of the matrix:\n");
   for (_i=0; _i<nrows; _i++) {
      for (_j=0; _j<ncols; _j++) {
         printf(" %6.4f ", Mr[_j*nrows + _i]);
         if (_j==ncols-1) printf("\n");
      }
      if (_i==nrows-1) printf("\n");
   }

   printf("\n\nImaginary part of the matrix:\n");
   for (_i=0; _i<nrows; _i++) {
      for (_j=0; _j<ncols; _j++) {
         printf(" %6.4f ", Mi[_j*nrows + _i]);
         if (_j==ncols-1) printf("\n");
      }
      if (_i==nrows-1) printf("\n");
   }
}

void PrintCellvec_f(const TSdcellvec *x_dcv) {
   int _i, ci, _n;
   double *v;

   if (!x_dcv) fn_DisplayError(".../fn_filesetup.c/PrintCellvec_f(): Input cell must be created (memory-allocated)");
   for (ci=0; ci<x_dcv->ncells; ci++ ) {
      _n = x_dcv->C[ci]->n;
      v = x_dcv->C[ci]->v;
      printf("\nCellvec %d:\n", ci);
      for (_i=0; _i<_n; _i++) {
         printf("v[%d]=%6.4f\n", _i, v[_i]);
      }
   }
}
void PrintCell_f(const TSdcell *x_dc) {
   int _i, _j, ci, nrows, ncols;
   double *M;

   if (!x_dc) fn_DisplayError(".../fn_filesetup.c/PrintCell_f(): Input cell must be created (memory-allocated)");
   for (ci=0; ci<x_dc->ncells; ci++ ) {
      nrows = x_dc->C[ci]->nrows;
      ncols = x_dc->C[ci]->ncols;
      M = x_dc->C[ci]->M;

      printf("\nCell %d:\n", ci);
      for (_i=0; _i<nrows; _i++) {
         for (_j=0; _j<ncols; _j++) {
            printf(" %6.4f ", M[_j*nrows + _i]);
            if (_j==ncols-1) printf("\n");
         }
         if (_i==nrows-1) printf("\n");
      }
   }
}


void PrintCell(const TSdcell *x_dc, const char *format)
{
   int _i, _j, ci, nrows, ncols;
   double *M;

   if (!x_dc) fn_DisplayError(".../fn_filesetup.c/PrintCell_f(): Input cell must be created (memory-allocated)");
   for (ci=0; ci<x_dc->ncells; ci++ ) {
      nrows = x_dc->C[ci]->nrows;
      ncols = x_dc->C[ci]->ncols;
      M = x_dc->C[ci]->M;

      printf("\nCell %d:\n", ci);
      if (!format)  format=" %10.5f ";   //Default format.
      for (_i=0; _i<nrows; _i++) {
         for (_j=0; _j<ncols; _j++) {
            printf(format, M[_j*nrows + _i]);
            if (_j==ncols-1) printf("\n");
         }
         if (_i==nrows-1) printf("\n");
      }
   }
}


void PrintFourthvec_f(TSdfourthvec *x_d4v) {
   int _j, _i, _k, _m, _n, _o;
   if (!x_d4v) fn_DisplayError(".../fn_filesetup.c/PrintFourthvec_f(): Input fourthvec cannot be NULL (must be created)");
   _m = x_d4v->ndims;
   for (_j=0; _j<_m; _j++) {
      _n = x_d4v->F[_j]->ncells;
      for (_i=0; _i<_n; _i++) {
         printf("\nFourthvec %d and Cell %d:\n", _j, _i);
         _o = x_d4v->F[_j]->C[_i]->n;
         for (_k=0; _k<_o; _k++) {
            printf("v[%d]=%6.4f\n", _k, x_d4v->F[_j]->C[_i]->v[_k]);
         }
      }
   }
}




//-------------------
// Prints entire input data (fptr_in) to the output file (fptr_out)
//   for the user to know what has produced the output.
//   The maximum number of characters in each line of the input file
//   is 4095 (excluding the NUL byte), but the rest of the line will
//   continue to be printed in new lines in the output file.
//-------------------
#define BUFFERLEN 4096
void ReprintInputData(FILE *fptr_in, FILE *fptr_out)
{
   char *inpbuffer;

   inpbuffer = tzMalloc(BUFFERLEN, char);   //@ Allocate memory to the string (including the NUL byte).
   rewind(fptr_in);
   while (fgets(inpbuffer,BUFFERLEN,fptr_in))
      fprintf(fptr_out, "%s", inpbuffer);
   fprintf(fptr_out, "\n\n\n\n\n//------------------------------- Output Data Begin Here -------------------------------\n");
   free(inpbuffer);
}
#undef BUFFERLEN

