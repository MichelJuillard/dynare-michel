/*=========================================================
 * Optimization package for different third-party routines, including csminwel.
 *
=========================================================*/
#include "optpackage.h"

#include "modify_for_mex.h"

#define STRLEN 256
static char filename_sp_vec_minproj[STRLEN];

static struct TSetc_csminwel_tag *CreateTSetc_csminwel(FILE *fptr_input1, const int n,  const int q, const int k);  //Used by CreateTSminpack() only.
static struct TSetc_csminwel_tag *DestroyTSetc_csminwel(struct TSetc_csminwel_tag *etc_csminwel_ps);  //Used by DestroyTSminpack() only.
//------- For csminwel only. -------
static TSminpack *SetMincsminwelGlobal(TSminpack *minpack_csminwel_ps);
static double minobj_csminwelwrap(double *x, int n, double **dummy1, int *dummy2);
static int mingrad_csminwelwrap(double *x, int n, double *g, double **dummy1, int *dummy2);
//------- For IMSL linearly constrainted optimization only. -------
static double GLB_FVALMIN = NEARINFINITY;  //Must be initialized to ba a very big number.
static int GLB_DISPLAY = 1;     //Print out intermediate results on screen.
static TSdvector *XIMSL_DV = NULL;  //To save the minimized value in case the IMSL quits with a higher value.
static struct TStateModel_tag *SetModelGlobalForIMSLconlin(struct TStateModel_tag *smodel_ps);
static void ObjFuncForModel_imslconlin(int d_x0, double *x0_p, double *fret_p);
static void imslconlin_SetPrintFile(char *filename);
static double opt_logOverallPosteriorKernal(struct TStateModel_tag *smodel_ps, TSdvector *xchange_dv);
static void gradcd_imslconlin(int n, double *x, double *g);
static double ObjFuncForModel_congrad(double *x0_p, int d_x0);



////TSminpack *CreateTSminpack(TFminpackage *minfinder_func, TFminobj *minobj_func, TFmingrad *mingrad_func, TFSetPrintFile *printinterresults_func, const int n, const int package)  //, const int indxAnag)
////TSminpack *CreateTSminpack(TFminfinder *minfinder_func, TFminobj *minobj_func, TFmingrad *mingrad_func, char *filename_printout, const int n, const int package)  //, const int indxAnag)
TSminpack *CreateTSminpack(TFminobj *minobj_func, void **etc_project_pps, TFmindestroy_etcproject *etcproject_func, TFmingrad *mingrad_func, char *filename_printout, const int n, const int package)
{
   TSminpack *minpack_ps = tzMalloc(1, TSminpack);

   //$$$$WARNING: Note the vector xtemp_dv->v or gtemp_dv-v itself is not allocated memory, but only the POINTER.
   //$$$$           Within the minimization routine like csminwel(), the temporary array x enters as the argument in
   //$$$$           the objective function to compare with other values.  If we use minpack_ps->x_dv->v = x
   //$$$$           in a wrapper function like minobj_csminwelwrap() where x is a temporay array in csminwel(),
   //$$$$           this tempoary array (e.g., x[0] in csminwel()) within the csminwel minimization routine
   //$$$$           will be freed after the csminwel minimization is done.  Consequently, minpack_ps->x_dv-v, which
   //$$$$           which was re-pointed to this tempoary array, will freed as well.  Thus, no minimization results
   //$$$$           would be stored and trying to access to minpack_ps->x_dv would cause memory leak.
   //$$$$           We don't need, however, to create another temporary pointer within the objective function itself,
   //$$$$           but we must use minpack_ps->xtemp_dv for a *wrapper* function instead and at the end of
   //$$$$           minimization, minpack_ps->x_dv will have the value of minpack_ps->xtemp_dv, which is automatically
   //$$$$           taken care of by csminwel with the lines such as
   //$$$$                 memcpy(xh,x[3],n*sizeof(double));
   //$$$$           where xh and minpack_ps->x_dv->v point to the same memory space.


   minpack_ps->xtemp_dv = tzMalloc(1, TSdvector);
   minpack_ps->gtemp_dv = tzMalloc(1, TSdvector);
   minpack_ps->xtemp_dv->flag = minpack_ps->gtemp_dv->flag = V_DEF;  //Set the flag first but will be assigned legal values in minobj_csminwelwrap().
   minpack_ps->xtemp_dv->n = minpack_ps->gtemp_dv->n = n;


   minpack_ps->x_dv = CreateVector_lf(n);
   minpack_ps->g_dv = CreateVector_lf(n);
   minpack_ps->x0_dv = CreateVector_lf(n);
   //
   minpack_ps->etc_project_ps  = (void *)*etc_project_pps;
   minpack_ps->DestroyTSetc_project = etcproject_func;
   if (etcproject_func) *etc_project_pps = NULL;   //If destroy function makes this structure responsible to free memory of the passing pointer, reset this passing pointer to NULL to avoid double destroying actions and cause memory problem.
   //
   minpack_ps->etc_package_ps = NULL;
   minpack_ps->minobj = minobj_func;
   minpack_ps->mingrad = mingrad_func;
   minpack_ps->filename_printout = filename_printout;
//   minpack_ps->SetPrintFile = printinterresults_func;

   if ( (minpack_ps->package=package) & MIN_CSMINWEL ) {
      minpack_ps->etc_package_ps = (void *)CreateTSetc_csminwel(((struct TSetc_minproj_tag *)minpack_ps->etc_project_ps)->args_blockcsminwel_ps->fptr_input1, n, 0, 0);
//      if (minpack_ps->mingrad)  fn_DisplayError(".../optpackage.c/CreateTSminpack():  Have not got time to deal with analytical gradient situation");
//      if (minpack_ps->indxAnag=indxAnag)  fn_DisplayError(".../optpackage.c/CreateTSminpack():  Have not got time to deal with analytical gradient situation");
   }
   else  fn_DisplayError(".../optpackage.c/CreateTSminpack():  Have not got time to specify other minimization packages than csminwel");

   return (minpack_ps);
}
//---
TSminpack *DestroyTSminpack(TSminpack *minpack_ps)
{
   if (minpack_ps) {
      //$$$$WARNING: Note the following vectors themselves are NOT allocated memory, but only the POINTERs.  Used within the minimization problem.
      //$$$$         See minobj_csminwelwrap() as an example.
      free(minpack_ps->xtemp_dv);
      free(minpack_ps->gtemp_dv);


      DestroyVector_lf(minpack_ps->x_dv);
      DestroyVector_lf(minpack_ps->g_dv);
      DestroyVector_lf(minpack_ps->x0_dv);
      if (minpack_ps->DestroyTSetc_project)  minpack_ps->DestroyTSetc_project(minpack_ps->etc_project_ps);   //If destroy function is active, destroy it here; ohterwise, it will be destroyed somewhere else.
      if ( minpack_ps->package & MIN_CSMINWEL )  DestroyTSetc_csminwel((TSetc_csminwel *)minpack_ps->etc_package_ps);

      //===
      free(minpack_ps);
      return ((TSminpack *)NULL);
   }
   else  return (minpack_ps);
}



//-----------------------------------------------------------------------
// Unconstrained BFGS csminwel package.
//-----------------------------------------------------------------------
static TSetc_csminwel *CreateTSetc_csminwel(FILE *fptr_input1, const int n, const int q, const int k)
{
   //If fptr_input1==NULL or no no values supplied when fptr_input1 != NULL, default values are taken.

   int _i;
   //===
   TSetc_csminwel *etc_csminwel_ps = tzMalloc(1, TSetc_csminwel);

   etc_csminwel_ps->_k = k;
   if (!k) {
      etc_csminwel_ps->args = (double **)NULL;
      etc_csminwel_ps->dims = (int *)NULL;
   }
   else {
      etc_csminwel_ps->dims = tzMalloc(k, int);
      etc_csminwel_ps->args = tzMalloc(k, double *);
      for (_i=k-1; _i>=0; _i--) *(etc_csminwel_ps->args + _i) = tzMalloc(q, double);
   }

   //=== Default values of input arguments.
   etc_csminwel_ps->Hx_dm = CreateMatrix_lf(n, n);  //n-by-n inverse Hessian.
   //+
   etc_csminwel_ps->badg = 1;   //1: numerical gradient will be used.
   etc_csminwel_ps->indxnumgrad_csminwel = INDXNUMGRAD_CSMINWEL;  //Method of the numerical gradient.

   //=== Reads doubles.
   if ( !fptr_input1 || !fn_SetFilePosition(fptr_input1, "//== crit ==//") || fscanf(fptr_input1, " %lf ", &etc_csminwel_ps->crit) != 1 )
      etc_csminwel_ps->crit = CRIT_CSMINWEL;   //Defaut for overall convergence criterion for the function value.
   if ( !fptr_input1 || !fn_SetFilePosition(fptr_input1, "//== ini_h_csminwel ==//") || fscanf(fptr_input1, " %lf ", &etc_csminwel_ps->ini_h_csminwel) != 1 )
      etc_csminwel_ps->ini_h_csminwel = INI_H_CSMINWEL;   //Defaut
   if ( !fptr_input1 || !fn_SetFilePosition(fptr_input1, "//== gradstps_csminwel ==//") || fscanf(fptr_input1, " %lf ", &etc_csminwel_ps->gradstps_csminwel) != 1 )
      etc_csminwel_ps->gradstps_csminwel = GRADSTPS_CSMINWEL;  //Default for step size of the numerical gradient.

   //=== Reads integers.
   if ( !fptr_input1 || !fn_SetFilePosition(fptr_input1, "//== itmax ==//") || fscanf(fptr_input1, " %d ", &etc_csminwel_ps->itmax) != 1 )
      etc_csminwel_ps->itmax = ITMAX_CSMINWEL;  //Default for maximum number of iterations.


   return (etc_csminwel_ps);
}
//#undef CRIT_CSMINWEL
//#undef ITMAX_CSMINWEL
//---
static TSetc_csminwel *DestroyTSetc_csminwel(TSetc_csminwel *etc_csminwel_ps)
{
   int _i;

   if (etc_csminwel_ps) {
      for (_i=etc_csminwel_ps->_k-1; _i>=0; _i--)  tzDestroy(etc_csminwel_ps->args[_i]);
      tzDestroy(etc_csminwel_ps->args);
      tzDestroy(etc_csminwel_ps->dims);
      //---
      DestroyMatrix_lf(etc_csminwel_ps->Hx_dm);

      //===
      free(etc_csminwel_ps);
      return ((TSetc_csminwel *)NULL);
   }
   else  return (etc_csminwel_ps);
}



/*********************************************
 * WARNING:  All the following data structures are declared global because
 *           (1) the minimization package takes only global variables;
 *           (2) these global structures make the existing functions reusable;
 *           (3) modifying the exisiting functions to keep global variables at minimum is NOT really worth the time.
*********************************************/
//---------------------------------
// Begin: This wrapper function makes it conformable to the call of the csminwel package.
//---------------------------------
static struct TSminpack_tag *MINPACK_CSMINWEL_PS = NULL;    //Minimization to find the MLE or posterior peak.
static TSminpack *SetMincsminwelGlobal(TSminpack *minpack_csminwel_ps)
{
   //Returns the old pointer in order to preserve the previous value.
   TSminpack *tmp_ps = MINPACK_CSMINWEL_PS;
   MINPACK_CSMINWEL_PS = minpack_csminwel_ps;
   return (tmp_ps);
}
static double minobj_csminwelwrap(double *x, int n, double **dummy1, int *dummy2)
{
   if (!MINPACK_CSMINWEL_PS || !MINPACK_CSMINWEL_PS->minobj)  fn_DisplayError(".../optpackage.c/minobj_csminwelwrap(): (1) MINPACK_CSMINWEL_PS must be created and (2) there exists an objective function assigned to MINPACK_CSMINWEL_PS->minobj");
   //  if (MINPACK_CSMINWEL_PS->x_dv->n != n)  fn_DisplayError(".../optpackage.c/minobj_csminwelwrap(): Length of passing vector must match minpack_ps->x_dv");
   MINPACK_CSMINWEL_PS->xtemp_dv->v = x;
   return (MINPACK_CSMINWEL_PS->minobj(MINPACK_CSMINWEL_PS));    //This function is specified in the main program.
}
//---
static int mingrad_csminwelwrap(double *x, int n, double *g, double **dummy1, int *dummy2)
{
   if (!MINPACK_CSMINWEL_PS || !MINPACK_CSMINWEL_PS->mingrad)  fn_DisplayError(".../optpackage.c/mingrad_csminwelwrap(): (1) MINPACK_CSMINWEL_PS must be created and (2) there exists an objective function assigned to MINPACK_CSMINWEL_PS->minobj");
   //  if (MINPACK_CSMINWEL_PS->x_dv->n != n)  fn_DisplayError(".../optpackage.c/mingrad_csminwelwrap(): Length of passing vector must match minpack_ps->x_dv");
   MINPACK_CSMINWEL_PS->xtemp_dv->v = x;
   MINPACK_CSMINWEL_PS->gtemp_dv->v = g;
   //>>>>>>>> Inside the following function, make sure to set MINPACK_CSMINWEL_PS->etc_csminwel_ps->badg = 0;   //1: numerical gradient will be used.
   MINPACK_CSMINWEL_PS->mingrad(MINPACK_CSMINWEL_PS);
   //<<<<<<<<

   return (0);
}
//---------------------------------
// End: This wrapper function makes it conformable to the call of the csminwel package.
//---------------------------------


//---------------------------------------------------------------//
//--- New ways to set up the minimization problems. 03/10/06. ---//
//---------------------------------------------------------------//
//------- Step 1. -------
//===
//=== Using blockwise csminwel minimization package.
struct TSargs_blockcsminwel_tag *CreateTSargs_blockcsminwel(FILE *fptr_input1)
{
   //If fptr_input1==NULL or no no values supplied when fptr_input1 != NULL, default values are taken.

   int nvec;
   struct TSargs_blockcsminwel_tag *args_blockcsminwel_ps = tzMalloc(1, struct TSargs_blockcsminwel_tag);


   //=== Reads doubles.
   if ( !fptr_input1 || !fn_SetFilePosition(fptr_input1, "//== criterion_start ==//") || fscanf(fptr_input1, " %lf ", &args_blockcsminwel_ps->criterion_start) != 1 )
      args_blockcsminwel_ps->criterion_start = 1.0e-3; //Default.
   if ( !fptr_input1 || !fn_SetFilePosition(fptr_input1, "//== criterion_end ==//") || fscanf(fptr_input1, " %lf ", &args_blockcsminwel_ps->criterion_end) != 1 )
      args_blockcsminwel_ps->criterion_end = 1.0e-6; //Default.
   if ( !fptr_input1 || !fn_SetFilePosition(fptr_input1, "//== criterion_increment ==//") || fscanf(fptr_input1, " %lf ", &args_blockcsminwel_ps->criterion_increment) != 1 )
      args_blockcsminwel_ps->criterion_increment = 0.1; //Default.
   if ( !fptr_input1 || !fn_SetFilePosition(fptr_input1, "//== max_iterations_increment ==//") || fscanf(fptr_input1, " %lf ", &args_blockcsminwel_ps->max_iterations_increment) != 1 )
      args_blockcsminwel_ps->max_iterations_increment = 1.5; //Default.
   if ( !fptr_input1 || !fn_SetFilePosition(fptr_input1, "//== ini_h_scale ==//") || fscanf(fptr_input1, " %lf ", &args_blockcsminwel_ps->ini_h_scale) != 1 )
      args_blockcsminwel_ps->ini_h_scale = 5.0e-4; //Default.
   if ( !fptr_input1 || !fn_SetFilePosition(fptr_input1, "//== gradstps_csminwel_const ==//") || fscanf(fptr_input1, " %lf ", &args_blockcsminwel_ps->gradstps_csminwel_const) != 1 )
      args_blockcsminwel_ps->gradstps_csminwel_const = 1.0e-4; //Default.

   //=== Reads integers.
   if ( !fptr_input1 || !fn_SetFilePosition(fptr_input1, "//== max_iterations_start ==//") || fscanf(fptr_input1, " %d ", &args_blockcsminwel_ps->max_iterations_start) != 1 )
      args_blockcsminwel_ps->max_iterations_start = 50; //Default.
   if ( !fptr_input1 || !fn_SetFilePosition(fptr_input1, "//== max_block_iterations ==//") || fscanf(fptr_input1, " %d ", &args_blockcsminwel_ps->max_block_iterations) != 1 )
      args_blockcsminwel_ps->max_block_iterations = 70; //Default.

   //=== Reads vectors.
   if (fptr_input1 && fn_SetFilePosition(fptr_input1, "//== gradstps_csminwel_dv ==//"))
   {
      if ( fscanf(fptr_input1, " %d ", &nvec) != 1)
         fn_DisplayError(".../fwz_comfuns.c/CreateTSinput(): check the first integer in the first row below the line //== gradstps_csminwel_dv ==// in the input data file");
      args_blockcsminwel_ps->gradstps_csminwel_dv = CreateVector_lf(nvec);
      if ( !ReadVector_lf(fptr_input1, args_blockcsminwel_ps->gradstps_csminwel_dv) )
         fn_DisplayError(".../fwz_comfuns.c/CreateTSinput(): check the data matrix or vector after the first row below the line //== gradstps_csminwel_dv ==// in the input data file");
      args_blockcsminwel_ps->gradstps_csminwel_dv->flag = V_DEF;
   }
   else  //Default (hard-coded).  fn_DisplayError(".../fwz_comfuns.c/CreateTSinput(): the line with //== gradstps_csminwel_dv ==// in the input data file does not exist");
   {
      args_blockcsminwel_ps->gradstps_csminwel_dv = CreateVector_lf(3);
      args_blockcsminwel_ps->gradstps_csminwel_dv->v[0] = 1.0e-02;
      args_blockcsminwel_ps->gradstps_csminwel_dv->v[1] = 1.0e-03;
      args_blockcsminwel_ps->gradstps_csminwel_dv->v[2] = 1.0e-03;
      args_blockcsminwel_ps->gradstps_csminwel_dv->flag = V_DEF;
   }


   args_blockcsminwel_ps->fptr_input1 = fptr_input1;

   return (args_blockcsminwel_ps);
}
//---
struct TSargs_blockcsminwel_tag *DestroyTSargs_blockcsminwel(struct TSargs_blockcsminwel_tag *args_blockcsminwel)
{
   if (args_blockcsminwel)
   {
      //===
      free(args_blockcsminwel);
      return ((struct TSargs_blockcsminwel_tag *)NULL);
   }
   else
      return (args_blockcsminwel);
}
//===
//=== Sets up a project-specific structure.
struct TSetc_minproj_tag *CreateTSetc_minproj(struct TStateModel_tag **smodel_pps, TFDestroyTStateModel *DestroyTStateModel_func,
               struct TSargs_blockcsminwel_tag **args_blockcsminwel_pps, struct TSargs_blockcsminwel_tag *(*DestroyTSargs_blockcsminwel)(struct TSargs_blockcsminwel_tag *))
{
   struct TSetc_minproj_tag *etc_minproj_ps = tzMalloc(1, struct TSetc_minproj_tag);

   //=== Initialization.
   etc_minproj_ps->smodel_ps = *smodel_pps;
   etc_minproj_ps->DestroyTStateModel = DestroyTStateModel_func;
   if (DestroyTStateModel_func)  *smodel_pps = (struct TStateModel_tag *)NULL;
        //If destroy function makes this structure responsible to free memory of the passing pointer, reset this passing pointer to NULL to avoid double destroying actions and cause memory problem.
        //  In this case, the original pointer *smodel_pps or smodel_ps is no longer valid, while etc_minproj_ps->smodel_ps.
        //  Note that we pass **smodel_pps only when we want to use DestroyTStateModel_func and let this structure take over smodel_ps.
        //  In many other cases, we do not need to pass **smodel_pps, but only *smodel_ps will do.
   //+
   etc_minproj_ps->args_blockcsminwel_ps = *args_blockcsminwel_pps;
   etc_minproj_ps->DestroyTSargs_blockcsminwel = DestroyTSargs_blockcsminwel;
   if (DestroyTSargs_blockcsminwel)  *args_blockcsminwel_pps = (struct TSargs_blockcsminwel_tag *)NULL;
        //If destroy function makes this structure responsible to free memory of the passing pointer, reset this passing pointer to NULL to avoid double destroying actions and cause memory problem.

   return (etc_minproj_ps);
}
//---
struct TSetc_minproj_tag *DestroyTSetc_minproj(struct TSetc_minproj_tag *etc_minproj_ps)
{
   if (etc_minproj_ps)
   {
      if (etc_minproj_ps->DestroyTStateModel)  etc_minproj_ps->DestroyTStateModel(etc_minproj_ps->smodel_ps);
             //If destroy function is active, destroy it here; ohterwise, it will be destroyed somewhere else.
      if (etc_minproj_ps->DestroyTSargs_blockcsminwel)  etc_minproj_ps->DestroyTSargs_blockcsminwel(etc_minproj_ps->args_blockcsminwel_ps);
             //If destroy function is active, destroy it here; ohterwise, it will be destroyed somewhere else.

      //===
      free(etc_minproj_ps);
      return ((struct TSetc_minproj_tag *)NULL);
   }
   else  return (etc_minproj_ps);
}
//------- Step 2. -------
//$$$$$$ 28/Oct/2007: I commented them out because it'd better left to be the user's function because of
//$$$$$$                (1) constant-parameter case without using DW's functions;
//$$$$$$                (2) allowing us to generate parameters randomly, which depends on the specific model.
//$$$$$$ See lwz_est.c in D:\ZhaData\WorkDisk\LiuWZ\Project2_empirical\EstimationOct07
//$$$$$$  or ExamplesForC.prn in D:\ZhaData\CommonFiles\C_Examples_DebugTips.
/**
void InitializeForMinproblem(struct TSminpack_tag *minpack_ps, char *filename_sp, TSdvector *gphi_dv, int indxStartValuesForMin)
{
   //Outputs:
   //  minpack_ps->x_dv and minpack_ps->xtemp_dv:
   //    The 1st gphi_dv->n elements of x_dv are model parameters (excluding those in the transition matrices).
   //    The 2nd-part or rest of the elements of x_dv are the free parameters in the transition matrices.
   //Inputs:
   //  gphi_dv: model free parameters (excluding those in the transition matrices);
   //  indxStartValuesForMin (corresponding to the command option /c in runprog.bat):
   //    0: continuing from the last estimated results contained in filename_sp.
   //    1: starts from the fixed values for gphi_dv, manually keyed in datainpu_setup.prn.
   //    2: randomly or arbitarily selects the initial starting values for the MLE or posterior estimate.
   FILE *fptr_sp = NULL;
   int _n, _i;
   int nqs;
   TSdvector xphi_sdv, xqs_sdv;
   TSdvector *x_dv = minpack_ps->x_dv;
   TSdvector *x0_dv = minpack_ps->x0_dv;
   //---
   struct TStateModel_tag *smodel_ps = (struct TStateModel_tag *)((struct TSetc_minproj_tag *)minpack_ps->etc_project_ps)->smodel_ps;
   int nfreempars = smodel_ps->routines->pNumberFreeParametersTheta(smodel_ps);

   if ( nfreempars != gphi_dv->n )
      fn_DisplayError("optpackage.c/InitializeForMinproblem(): Input vector gphi_dv must be free model parameters only");
   if ( nqs=NumberFreeParametersQ(smodel_ps) != x_dv->n - nfreempars )
      fn_DisplayError("optpackage.c/InitializeForMinproblem(): Minimization vector must have length equal to # of free model parameters plus # of free transition matrix parameters");

   xphi_sdv.flag = V_DEF;
   xphi_sdv.n = nfreempars;
   xphi_sdv.v = x_dv->v;

   xqs_sdv.flag = V_DEF;
   xqs_sdv.n = nqs;
   xqs_sdv.v = x_dv->v + xphi_sdv.n;

   if (indxStartValuesForMin == 1)
   {
      CopyVector0(&xphi_sdv, gphi_dv);
      ConvertQToFreeParameters(smodel_ps, xqs_sdv.v);     //Waggnoer's own function for the transition matrix.
      x_dv->flag = V_DEF;
   }
   else if (!indxStartValuesForMin)
   {
      fptr_sp = tzFopen(filename_sp,"r");
      rewind(fptr_sp);   //Must put the pointer at the beginning of the file.

      for (_n=x_dv->n, _i=0; _i<_n; _i++)
         if (fscanf(fptr_sp, " %lf ", x_dv->v+_i) != 1)
         {
            printf("Error: optpackage.c/InitializeForMinproblem() -- cannot read the number from the file %s.  Check the data file", filename_sp);
            exit(EXIT_FAILURE);
         }
      x_dv->flag = V_DEF;

      tzFclose(fptr_sp);
   }
   else  fn_DisplayError("optpackage.c/InitializeForMinproblem(): the case indxStartValuesForMin = 2 has not been programmed yet");


   //--- Initial or starting values of the parameters.
   CopyVector0(x0_dv, x_dv);
   SetupObjectiveFunction(smodel_ps, xphi_sdv.v, xqs_sdv.v, xphi_sdv.v);  //Must before using PosteriorObjectiveFunction();
   minpack_ps->fret0 = minpack_ps->fret = PosteriorObjectiveFunction(xphi_sdv.v, xphi_sdv.n);
//   minpack_ps->fret0 = minpack_ps->fret = -logOverallPosteriorKernal(smodel_ps, x0_dv);
   if (minpack_ps->fret0 >= NEARINFINITY)
   {
      printf("\nFatal Error:\n");
      printf("  optpackage.c/InitializeForMinproblem(): Bad initialization. All parameters must be in the reasonable range.\n");
//      printf("  optpackage.c/InitializeForMinproblem(): Bad initialization. All parameters must be in the reasonable range.\n"
//             "    Most likely, the parameters get stuck in the following line in swz2_confuns.c:\n"
//             " if ((tmpd1mPhi=1.0-fn_normalcdf(xid * (log(mt-boundthetamt_1) - logdbar))) <= 0.0)  logvalue = -NEARINFINITY;\n");
//      exit(EXIT_FAILURE);
   }
}
/**/
//------- Step 3. -------
void minfinder_blockcsminwel(struct TSminpack_tag *minpack_ps, int indx_findMLE)
{
   //Better version (November 2007)
   //Inputs:
   //  indx_findMLE: 1: find MLE without a prior, 0: find posterior (with a prior).

   //--- Block-csminwel arguments.
   struct TSargs_blockcsminwel_tag *args_blockcsminwel_ps = ((struct TSetc_minproj_tag *)minpack_ps->etc_project_ps)->args_blockcsminwel_ps;
   //--- DW's Markov-switching structure.
   struct TStateModel_tag *smodel_ps = ((struct TSetc_minproj_tag *)minpack_ps->etc_project_ps)->smodel_ps;
   //--- TSminpack arguments.
   TSdvector *x_dv = minpack_ps->x_dv;
   TSdvector *g_dv = minpack_ps->g_dv;
   char *filename_printout = minpack_ps->filename_printout;  //Printing out the intermediate results of x_dv and g_dv.
   double fret = minpack_ps->fret; //Returned value of the objective function.
   struct TSetc_csminwel_tag *etc_csminwel_ps = (TSetc_csminwel *)minpack_ps->etc_package_ps;
   //--- Blockwise arguments.
   int n1, n2;
   int _n = x_dv->n;
   double *x1_pd, *x2_pd, *g1_pd, *g2_pd;
   double fret_last, logvalue;
   TSdvector *gradstps_csminwel_dv = args_blockcsminwel_ps->gradstps_csminwel_dv;
   //--- Blockwise csminwel intput arguments.
   double criterion_start = args_blockcsminwel_ps->criterion_start;
   double criterion_end =  args_blockcsminwel_ps->criterion_end;
   double criterion_increment =  args_blockcsminwel_ps->criterion_increment;
   int max_iterations_start = args_blockcsminwel_ps->max_iterations_start;
   double max_iterations_increment = args_blockcsminwel_ps->max_iterations_increment;
   int max_block_iterations = args_blockcsminwel_ps->max_block_iterations;
   double ini_h_csminwel = args_blockcsminwel_ps->ini_h_scale;
   //+ Other csminwel arguments
   int iteration, total_iteration;
   int niters, fcount, retcodeh, max_niters;
   double crit;
   //=== Blockwise and overall memory creations.
   TSdmatrix *H1_dm = NULL;
   TSdmatrix *H2_dm = NULL;
   TSdmatrix *H_dm = NULL;
   //
   FILE *fptr_interesults = (FILE *)NULL;   //Printing intermediate results to a file.



   if (!x_dv || !x_dv->flag)  fn_DisplayError("swz2_comfuns.c/ minfinder_blockcsminwel(): free parameters x_dv must be initialized");

   n1 = NumberFreeParametersTheta(smodel_ps);       //Number of free model parameters.
   n2 = NumberFreeParametersQ(smodel_ps);   //Number of free transition matrix elements.
   if (_n != (n1 + n2))  fn_DisplayError("optpackage.c/minfinder_blockcsminwel(): total number of free parameters"
              "  must be equal to number of free model parameters + number of free q's");
   H1_dm = CreateMatrix_lf(n1, n1);
   H2_dm = CreateMatrix_lf(n2, n2);
   H_dm = CreateMatrix_lf(_n, _n);
   //
   x1_pd = x_dv->v;
   x2_pd = x_dv->v+n1;
   g1_pd = g_dv->v;
   g2_pd = g_dv->v+n1;

   //---- Refreshing the parameters outside this function.  TZ October 2007.
   SetupObjectiveFunction(smodel_ps, x1_pd, x2_pd, x_dv->v);
   logvalue = -( minpack_ps->fret0 = minpack_ps->fret = PosteriorObjectiveFunction(x_dv->v, x_dv->n) );  //Refreshing. logPosterirPdf.  DW function.
   fprintf(FPTR_OPT, "\n=========== Beginning Blockwise and Overall csminwel Minimizations =======================\nLog Peak Value: %.16e\n", logvalue);
   fflush(FPTR_OPT);


   //======= Minimizing using csminwel =======
   //--- Set up a printout file to record x_dv and g_dv.
   csminwel_SetPrintFile(filename_printout);  //Set the print-out file outputsp_mle_tag.prn.
   for (total_iteration=1, crit=criterion_start, max_niters=max_iterations_start;
        crit >= criterion_end;
        crit*=criterion_increment, max_niters=(int)(max_niters*max_iterations_increment))
   {
      for (iteration=1; iteration <= max_block_iterations; total_iteration++, iteration++)
      {
         fret_last = fret;
         //=== Minimizing the objective function w.r.t. the 1st block of parameters (model parameters).
         printf("\nMinimizing user's specific model parameters at iteration %d\n",iteration);
         InitializeDiagonalMatrix_lf(H1_dm, ini_h_csminwel);
         H1_dm->flag = M_GE | M_SU | M_SL;       //Hessian is symmetric.
         //+
         SetupObjectiveFunction(smodel_ps, x1_pd, x2_pd, x_dv->v);
         GRADSTPS_CSMINWEL = gradstps_csminwel_dv->v[0];
         if (indx_findMLE)
            csminwel(MLEObjectiveFunction_csminwel, x1_pd, n1, H1_dm->M, g1_pd, NULL,
                     &fret, crit, &niters, max_niters, &fcount, &retcodeh,
                     (double **)NULL, (int *)NULL);
         else
            csminwel(PosteriorObjectiveFunction_csminwel, x1_pd, n1, H1_dm->M, g1_pd, NULL,
                     &fret, crit, &niters, max_niters, &fcount, &retcodeh,
                     (double **)NULL, (int *)NULL);

         ConvertFreeParametersToQ(smodel_ps,x2_pd);
         ConvertFreeParametersToTheta(smodel_ps,x1_pd);

         //+
         logvalue = -fret;
         fprintf(FPTR_OPT, "\n=========== Block iteration %d for block 1 at total iteration %d =======================\nLog Peak Value: %.16e\n", iteration, total_iteration, logvalue);
         fflush(FPTR_OPT);



         //=== Minimizing the objective function w.r.t. the 2nd block of parameters (transition matrix).
         printf("\nMinimizing transitiona matrix Q at iteration %d\n",iteration);
         InitializeDiagonalMatrix_lf(H2_dm, ini_h_csminwel);
         H2_dm->flag = M_GE | M_SU | M_SL;       //Hessian is symmetric.
         //+
         SetupObjectiveFunction(smodel_ps, x2_pd, x2_pd, x_dv->v);
         GRADSTPS_CSMINWEL = gradstps_csminwel_dv->v[1];
         if (indx_findMLE)
            csminwel(MLEObjectiveFunction_csminwel, x2_pd, n2, H2_dm->M, g2_pd, NULL,
                     &fret, crit, &niters, max_niters, &fcount, &retcodeh,
                     (double **)NULL, (int *)NULL);
         else
            csminwel(PosteriorObjectiveFunction_csminwel, x2_pd, n2, H2_dm->M, g2_pd, NULL,
                     &fret, crit, &niters, max_niters, &fcount, &retcodeh,
                     (double **)NULL, (int *)NULL);

         ConvertFreeParametersToQ(smodel_ps,x2_pd);
         ConvertFreeParametersToTheta(smodel_ps,x1_pd);

         //+
         logvalue = -fret;
         fprintf(FPTR_OPT, "\n=========== Block iteration %d for block 2 at total iteration %d =======================\nLog Peak Value: %.16e\n", iteration, total_iteration, logvalue);
         fprintf(FPTR_OPT, "--------Numerical gradient---------\n");
         WriteVector(FPTR_OPT, g_dv, " %0.16e ");
         fprintf(FPTR_OPT, "--------Restarting point---------\n");
         WriteVector(FPTR_OPT, x_dv, " %0.16e ");
         fflush(FPTR_OPT);


         if (fabs(fret - fret_last) <= crit)  break;
      }

      //=== Minimizing the overall likelihood or posterior kernel.
      logvalue = -fret;
      fprintf(FPTR_OPT,"\n\n=========== Total iteration %d ===========\n",++total_iteration);
      fprintf(FPTR_OPT,"Criterion/Max_Numer_Iterations:  %le  %d\n",crit,max_niters);
      fprintf(FPTR_OPT,"Log peak value before overall minimization:  %.16e\n", logvalue);
      fflush(FPTR_OPT);
      //---
      InitializeDiagonalMatrix_lf(H_dm, ini_h_csminwel);
      H_dm->flag = M_GE | M_SU | M_SL;       //Hessian is symmetric.
      //+
      SetupObjectiveFunction(smodel_ps, x_dv->v, x2_pd, x_dv->v);
      GRADSTPS_CSMINWEL = gradstps_csminwel_dv->v[2];
      if (indx_findMLE)
         csminwel(MLEObjectiveFunction_csminwel, x_dv->v, _n, H_dm->M, g_dv->v, NULL,
                  &fret, crit, &niters, max_niters, &fcount, &retcodeh,
                  (double **)NULL, (int *)NULL);
      else
         csminwel(PosteriorObjectiveFunction_csminwel, x_dv->v, _n, H_dm->M, g_dv->v, NULL,
                  &fret, crit, &niters, max_niters, &fcount, &retcodeh,
                  (double **)NULL, (int *)NULL);


      //---
      logvalue = -fret;
      fprintf(FPTR_OPT,"Log peak value after overall minimization:  %.16e\n", logvalue);
      fprintf(FPTR_OPT, "--------Numerical gradient---------\n");
      WriteVector(FPTR_OPT, g_dv, " %0.16e ");
      fprintf(FPTR_OPT, "--------Restarting point---------\n");
      WriteVector(FPTR_OPT, x_dv, " %0.16e ");
      fflush(FPTR_OPT);

      //--- Write to the intermediate results file.
      if ( !(fptr_interesults = fopen(filename_printout,"w")) ) {
         printf("\n\nUnable to open the starting point data file %s in minfinder_blockcsminwel() in optpackage.c!\n", filename_printout);
         getchar();
         exit(EXIT_FAILURE);
      }
      fprintf(fptr_interesults, "========= All blocks are reported here. ========== \n");
      fprintf(fptr_interesults, "--------Numerical gradient---------\n");
      WriteVector(fptr_interesults, g_dv, " %0.16e ");
      fprintf(fptr_interesults, "--------Restarting point---------\n");
      WriteVector(fptr_interesults, x_dv, " %0.16e ");
      fflush(fptr_interesults);
      tzFclose(fptr_interesults);


      ConvertFreeParametersToQ(smodel_ps,x2_pd);
      ConvertFreeParametersToTheta(smodel_ps,x1_pd);
   }

   etc_csminwel_ps->niter = niters;   //Number of iterations taken by csminwel.
   etc_csminwel_ps->fcount = fcount;   //Number of function evaluations used by csminwel.
   etc_csminwel_ps->retcode = retcodeh;  //Return code for the terminating condition.
                // 0, normal step (converged). 1, zero gradient (converged).
                // 4,2, back and forth adjustment of stepsize didn't finish.
                // 3, smallest stepsize still improves too slow. 5, largest step still improves too fast.
                // 6, no improvement found.

   DestroyMatrix_lf(H1_dm);
   DestroyMatrix_lf(H2_dm);
   DestroyMatrix_lf(H_dm);
}


//-----------------------------------------------------
// Minimization csminwel for the constant parameter model only.  5/24/04.
//-----------------------------------------------------
//------- Step 2. -------
//--- 28/Oct/07: This function has NOT been used even for the constant-parameter model.
//--- For examples, see lwz_est.c in D:\ZhaData\WorkDisk\LiuWZ\Project2_empirical\EstimationOct07
//---                or ExamplesForC.prn under D:\ZhaData\CommonFiles\C_Examples_DebugTips.
/**
void InitializeForMinproblem_const(struct TSminpack_tag *minpack_ps, char *filename_sp, TSdvector *gphi_dv, int indxStartValuesForMin)
{
   //Outputs:
   //  minpack_ps->x_dv and minpack_ps->xtemp_dv:
   //    The 1st gphi_dv->n elements of x_dv are model parameters (excluding those in the transition matrices).
   //    The 2nd-part or rest of the elements of x_dv are the free parameters in the transition matrices.
   //Inputs:
   //  gphi_dv: model free parameters (excluding those in the transition matrices);
   //  indxStartValuesForMin (corresponding to the command option /c in runprog.bat):
   //    0: continuing from the last estimated results contained in filename_sp.
   //    1: starts from the fixed values for gphi_dv, manually keyed in datainpu_setup.prn.
   //    2: randomly or arbitarily selects the initial starting values for the MLE or posterior estimate.
   FILE *fptr_sp = NULL;
   int _n, _i;
   TSdvector xphi_sdv;
   TSdvector *x_dv = minpack_ps->x_dv;
   TSdvector *x0_dv = minpack_ps->x0_dv;
   //---
   struct TStateModel_tag *smodel_ps = (struct TStateModel_tag *)((struct TSetc_minproj_tag *)minpack_ps->etc_project_ps)->smodel_ps;
   int nfreempars = smodel_ps->routines->pNumberFreeParametersTheta(smodel_ps);

   if ( nfreempars != gphi_dv->n )
      fn_DisplayError("optpackage.c/InitializeForMinproblem_const(): Input vector gphi_dv must be free model parameters only");

   xphi_sdv.flag = V_DEF;
   xphi_sdv.n = nfreempars;
   xphi_sdv.v = x_dv->v;

   if (indxStartValuesForMin == 1)
   {
      CopyVector0(&xphi_sdv, gphi_dv);
      x_dv->flag = V_DEF;
   }
   else if (!indxStartValuesForMin)
   {
      fptr_sp = tzFopen(filename_sp,"r");
      rewind(fptr_sp);   //Must put the pointer at the beginning of the file.

      for (_n=x_dv->n, _i=0; _i<_n; _i++)
         if (fscanf(fptr_sp, " %lf ", x_dv->v+_i) != 1)
         {
            printf("Error: optpackage.c/InitializeForMinproblem_const() -- cannot read the number from the file %s.  Check the data file", filename_sp);
            exit(EXIT_FAILURE);
         }
      x_dv->flag = V_DEF;

      tzFclose(fptr_sp);
   }
   else  fn_DisplayError("optpackage.c/InitializeForMinproblem_const(): the case indxStartValuesForMin = 2 has not been programmed yet");


   //--- Initial or starting values of the parameters.
   CopyVector0(x0_dv, x_dv);
   //The following line does not work because minpack_ps->xtemp_ps will be used in minneglogpost_const(), which has not be initialized.
   //Use instead minpack_ps->fret0 = minpack_ps->fret = logOverallPosteriorKernal_const(smodel_ps, minpack_ps->x0_dv);
   //minpack_ps->fret0 = minpack_ps->fret = minpack_ps->minobj(minpack_ps);  This will not work because
}
/**/

//------- Step 3. -------
void minfinder(TSminpack *minpack_ps)
{
   TSdvector *x_dv = minpack_ps->x_dv;
   //--- For MIN_CSMINWEL only.
   TSdmatrix *Hx_dm;
   TSetc_csminwel *etc_csminwel_ps;

   if (minpack_ps->package & MIN_CSMINWEL) {
      if (!x_dv->flag)  fn_DisplayError("optpackage.c/ minfinder(): Parameter x_dv must be initialized");
      else {
         //=== BFGS (csminwel) method.
         etc_csminwel_ps = (TSetc_csminwel *)minpack_ps->etc_package_ps;
         Hx_dm = etc_csminwel_ps->Hx_dm;
         //Alternative: Hx_dm = ((TSetc_csminwel *)minpack_ps->etc_package_ps)->Hx_dm;
         if (!Hx_dm->flag) {
            InitializeDiagonalMatrix_lf(Hx_dm, INI_H_CSMINWEL);
            Hx_dm->flag = M_GE | M_SU | M_SL;       //Hessian is symmetric.
         }
         //if (minpack_ps->filename_printout)  csminwel_SetPrintFile(minpack_ps->filename_printout);
         csminwel_SetPrintFile(minpack_ps->filename_printout);
         SetMincsminwelGlobal(minpack_ps);
         GRADSTPS_CSMINWEL = etc_csminwel_ps->gradstps_csminwel;
         csminwel(minobj_csminwelwrap, x_dv->v, x_dv->n, Hx_dm->M, minpack_ps->g_dv->v, minpack_ps->mingrad ? mingrad_csminwelwrap : NULL,
                  &minpack_ps->fret, etc_csminwel_ps->crit, &etc_csminwel_ps->niter, etc_csminwel_ps->itmax,
                  &etc_csminwel_ps->fcount, &etc_csminwel_ps->retcode, (double **)NULL, (int *)NULL);
      }
   }
   else  fn_DisplayError("optpackage.c/minfinder():  (1) minpack_ps must be created and (2) I have not got time to specify other minimization packages such as ISML");
}



//-----------------------------------------------------------------------
// Linearly-constrained IMSL minimization package.
//-----------------------------------------------------------------------
struct TSpackage_imslconlin_tag *CreateTSpackagae_imslconlin(const int npars_tot, const int neqs, const int ncons)
{
   //npars_tot: total number of variables (e.g., all the model variables plus free parameters in transition matrix).
   //ncons: total number of constraints (excluding simple bounds) which include the linear equality constraints.
   //neqs:  number of linear equality constrains.  Thus, ncons >= neqs.
   //lh_coefs_dv: ncons*npars_tot-by-1 left-hand-side constraint cofficients with the first neqs rows dealing with equality constraints.
   //rh_constraints_dv: ncons-by-1 right-hand-side the values for all the constraints.
   //lowbounds_dv: npars_tot-by-1 simple lower bounds.
   //upperbounds_dv: npars_tot-by-1 simple upper bounds.
   //===
   struct TSpackage_imslconlin_tag *package_imslconlin_ps = tzMalloc(1, struct TSpackage_imslconlin_tag);

   if (neqs > ncons || npars_tot<=0)
      fn_DisplayError("CreateTSpackage_imslconlin(): make sure (1) # of equality constraints no greater than total # of constraints"
                      "\t and (2) number of free parameters must be greater than 0");
   package_imslconlin_ps->npars_tot = npars_tot;
   package_imslconlin_ps->neqs = neqs;
   package_imslconlin_ps->ncons = ncons;


   if (ncons<=0)
   {
      package_imslconlin_ps->lh_coefs_dv = NULL;
      package_imslconlin_ps->rh_constraints_dv = NULL;
   }
   else
   {
      package_imslconlin_ps->lh_coefs_dv = CreateConstantVector_lf(ncons*npars_tot, 0.0);
      package_imslconlin_ps->rh_constraints_dv = CreateVector_lf(ncons);
   }
   package_imslconlin_ps->lowbounds_dv = CreateConstantVector_lf(npars_tot, -BIGREALNUMBER);
   package_imslconlin_ps->upperbounds_dv = CreateConstantVector_lf(npars_tot, BIGREALNUMBER);
   //-
   package_imslconlin_ps->xsaved_dv = CreateVector_lf(package_imslconlin_ps->npars_tot);
   XIMSL_DV = CreateVector_lf(package_imslconlin_ps->npars_tot);  //Used in ObjFuncForModel_imslconlin() to save the minimized value in case the IMSL quits with a higher value.

   package_imslconlin_ps->crit = CRIT_IMSLCONLIN;
   package_imslconlin_ps->itmax = ITMAX_IMSLCONLIN;


   return (package_imslconlin_ps);
}
//---
struct TSpackage_imslconlin_tag *DestroyTSpackagae_imslconlin(struct TSpackage_imslconlin_tag *package_imslconlin_ps)
{
   if (package_imslconlin_ps)
   {
      DestroyVector_lf(package_imslconlin_ps->lh_coefs_dv);
      DestroyVector_lf(package_imslconlin_ps->rh_constraints_dv);
      DestroyVector_lf(package_imslconlin_ps->lowbounds_dv);
      DestroyVector_lf(package_imslconlin_ps->upperbounds_dv);
      //
      DestroyVector_lf(package_imslconlin_ps->xsaved_dv);
      DestroyVector_lf(XIMSL_DV);

      //===
      free(package_imslconlin_ps);
      return ((struct TSpackage_imslconlin_tag *)NULL);
   }
   else  return (package_imslconlin_ps);
}
//-----------------------------------------------------------------------
// Using Linearly-constrained IMSL minimization package.
//-----------------------------------------------------------------------
void minfinder_noblockimslconlin(struct TSpackage_imslconlin_tag *package_imslconlin_ps, struct TSminpack_tag *minpack_ps, char *filename_printout, int ntheta)
{
   //ntheta: number of free model parameters (NOT including free transition matrix Q parameters).
   //filename_printout: the file that stores the intermediate results.

   //--- Model or project specific structure.
   struct TStateModel_tag *smodel_ps = ((struct TSetc_minproj_tag *)minpack_ps->etc_project_ps)->smodel_ps;
   //---
   TSdvector *x_dv = minpack_ps->x_dv;
   TSdvector *g_dv = minpack_ps->g_dv;
   double *x1_pd, *x2_pd;
   //===
   TSdvector *xguess_dv = CreateVector_lf(x_dv->n);

   x1_pd = x_dv->v;
   x2_pd = x_dv->v + ntheta;  //In the constant parameter model, this will point to invalid,
                              //  but will be taken care of automatically by DW's function ConvertFreeParametersToQ().

   CopyVector0(xguess_dv, x_dv);


   //======= IMSL linearly-constrained optimization, which makes sure that the boundary condition is met.
   imslconlin_SetPrintFile(filename_printout);  //Set the print-out file outputsp_min_tag.prn.
   printf("\n\n======= Starting the IMSL constrained optimization======= \n\n");
   fflush(stdout);
   //====== The following linearly-constrained minimization works well for this kind of model but has a bugger of returning a higher value of the objective function.
   CopyVector0(XIMSL_DV, x_dv); //This is absolutely necessary because once imsl_d_min_con_gen_lin() is called, x_dv will be
                                //  changed before ObjFuncForModel_imslconlin() is evaluated.  It is possible that x_dv is changed
                                //  so much that bad objective is returned and thus XIMSL_DV would be bad from the start, thus
                                //  giving 1.eE+3000 from beginning to end.
   GLB_FVALMIN = -LogPosterior_StatesIntegratedOut(smodel_ps);
   CopyVector0(package_imslconlin_ps->xsaved_dv, XIMSL_DV);
   //+
   SetModelGlobalForIMSLconlin(smodel_ps);
   if (imsl_d_min_con_gen_lin(ObjFuncForModel_imslconlin, x_dv->n, package_imslconlin_ps->ncons, package_imslconlin_ps->neqs,
                              package_imslconlin_ps->lh_coefs_dv->v,
                              package_imslconlin_ps->rh_constraints_dv->v,
                              package_imslconlin_ps->lowbounds_dv->v, package_imslconlin_ps->upperbounds_dv->v,
                              IMSL_XGUESS, xguess_dv->v, IMSL_GRADIENT, gradcd_imslconlin,
                              IMSL_MAX_FCN, package_imslconlin_ps->itmax, IMSL_OBJ, &minpack_ps->fret,
                              IMSL_TOLERANCE, package_imslconlin_ps->crit, IMSL_RETURN_USER, x_dv->v, 0))
   {
      printf("\nFinished: IMSL linearly-constrained optimization is successfully finished with the value of obj. fun.: %.16e\n", minpack_ps->fret);
   }
   else  printf("\nWarning: IMSL linearly-constrained optimization fails, so the results from csminwel and congramin are used.\n");
   printf("\n===Ending the IMSL constrained optimization===\n");

   //=== Printing out messages indicating that IMSL has bugs.
   if (minpack_ps->fret > GLB_FVALMIN)
   {
      //IMSL linearly-constrained optimization returns a higher obj. func.  (a bug).
      printf("\n----------IMSL linearly-constrained minimization finished but with a higher objective function value!----------\n");
      printf("The improperly-returned value is %.10f while the lowest value of the objective function is %.16e.\n\n", minpack_ps->fret, GLB_FVALMIN);
      fflush(stdout);
      fprintf(FPTR_DEBUG, "\n----------IMSL linearly-constrained minimization finished but with a higher objective function value!----------\n");
      fprintf(FPTR_DEBUG, "The improperly-returned value is %.16e while the lowest value of the objective function is %.16e.\n\n", minpack_ps->fret, GLB_FVALMIN);
      fflush(FPTR_DEBUG);
   }

   ConvertFreeParametersToQ(smodel_ps,x2_pd);
         //DW's function, which takes care of the degenerate case where x2_pd points to an
         //  invalid place as in the constant parameter case.
   ConvertFreeParametersToTheta(smodel_ps,x1_pd); //DW's function, which calls TZ's function.  So essentially it's TZ's function.

   //Saved the last best results in case the IMSL quits with a bug.
   CopyVector0(package_imslconlin_ps->xsaved_dv, XIMSL_DV);


   //===
   DestroyVector_lf(xguess_dv);
}
//===
static struct TStateModel_tag *SMODEL_PS = NULL;    //Minimization to find the MLE or posterior peak.
static struct TStateModel_tag *SetModelGlobalForIMSLconlin(struct TStateModel_tag *smodel_ps)
{
   //Returns the old pointer in order to preserve the previous value.
   struct TStateModel_tag *tmp_ps = SMODEL_PS;
   SMODEL_PS = smodel_ps;
   return (tmp_ps);
}
static void ObjFuncForModel_imslconlin(int d_x0, double *x0_p, double *fret_p)
{
   TSdvector x0_sdv;
   //
   FILE *fptr_startingpoint_vec = NULL;
   static int ncnt_fevals = -1;

   // printf("\n----- Entering the objective function. ------");
   // fflush(stdout);
   x0_sdv.v = x0_p;
   x0_sdv.n = d_x0;
   x0_sdv.flag = V_DEF;

   *fret_p = -opt_logOverallPosteriorKernal(SMODEL_PS, &x0_sdv);
   if ( GLB_DISPLAY) {
      printf("\nValue of objective function at the %dth evaluation: %.16e\n", ++ncnt_fevals, *fret_p);
      fflush(stdout);
   }
   if (*fret_p < GLB_FVALMIN) {
      //=== Resets GLB_FVALMIN at *fret_p and then prints the intermediate point to a file.
      fptr_startingpoint_vec = tzFopen(filename_sp_vec_minproj,"w");
      fprintf(fptr_startingpoint_vec, "================= Output from IMSC linear constrained optimization ====================\n");
      fprintf(fptr_startingpoint_vec, "IMSL: Value of objective miminization function at the %dth iteration: %.15f\n", ncnt_fevals, GLB_FVALMIN=*fret_p);
      fprintf(fptr_startingpoint_vec, "--------Restarting point---------\n");
      WriteVector(fptr_startingpoint_vec, &x0_sdv, " %0.16e ");
      CopyVector0(XIMSL_DV, &x0_sdv);  //Saved in case the IMSL quits with a bug.
      //=== Must print this results because imsl_d_min_con_gen_lin() has a bug and quits with a higher value. The printed-out results in the debug file may be used for imsl_d_min_con_nonlin() to continue.
      fprintf(FPTR_DEBUG, "\nIMSL: Value of objective miminization function at the %dth iteration: %.15f\n", ncnt_fevals, GLB_FVALMIN=*fret_p);
      fprintf(FPTR_DEBUG, "--------Restarting point---------\n");
      WriteVector(FPTR_DEBUG, &x0_sdv, " %0.16e ");
      fflush(FPTR_DEBUG);
   }

   //  printf("\n----- Leaving the objective function. ------\n");
   //  fflush(stdout);

   tzFclose(fptr_startingpoint_vec);
}
//------------------------
// Overall posterior kernal for calling Waggoner's regime-switching procedure.
//------------------------
static double opt_logOverallPosteriorKernal(struct TStateModel_tag *smodel_ps, TSdvector *xchange_dv)
{
   double *x1_pd, *x2_pd;

   x1_pd = xchange_dv->v;
   x2_pd = xchange_dv->v + NumberFreeParametersTheta(smodel_ps);
        //Note that NumberFreeParametersTheta() is DW's function, which points to TZ's function.
        //In the constant parameter model, this will point to invalid,
        //  but will be taken care of automatically by DW's function ConvertFreeParametersToQ().

   //======= This is a must step to refresh the value at the new point. =======
   ConvertFreeParametersToTheta(smodel_ps, x1_pd);   //Waggoner's function, which calls TZ's Convertphi2*().
   ConvertFreeParametersToQ(smodel_ps, x2_pd);   //Waggoner's function, which automatically takes care of the constant-parameter situition
   ThetaChanged(smodel_ps); //DW's function, which will also call my function to set a flag for refreshing everything under these new parameters.
   if (1)  //Posterior function.
      return ( LogPosterior_StatesIntegratedOut(smodel_ps) ); //DW's function.
   else //Likelihood (with no prior)
      return ( LogLikelihood_StatesIntegratedOut(smodel_ps) ); //DW's function.
}
//---
static void imslconlin_SetPrintFile(char *filename) {
   if (!filename)   sprintf(filename_sp_vec_minproj, "outdata5imslconlin.prn");  //Default filename.
   else if (STRLEN-1 < strlen(filename))  fn_DisplayError(".../optpackage.c:  the allocated length STRLEN for filename_sp_vec_minproj is too short.  Must increase the string length");
   else  strcpy(filename_sp_vec_minproj, filename);
}
//---
static void gradcd_imslconlin(int n, double *x, double *g)
{
   //Outputs:
   //  g: the gradient n-by-1 g (no need to be initialized).
   //Inputs:
   //  x: the vector point at which the gradient is evaluated.  No change in the end although will be added or
   //       substracted by dh during the function (but in the end the original value will be put back).
   //  n: the dimension of g or x.
   //int _i;
   FILE *fptr_startingpoint_vec = NULL;

   // printf("\n=== Entering the gradient function. ===\n");
   // fflush(stdout);
   GLB_DISPLAY = 0;   //This guarantees that the objective function printouts will not show when the gradient is computed.
   gradcd_gen(g, x, n, ObjFuncForModel_congrad, (double *)NULL, ObjFuncForModel_congrad(x, n));

   //=== Prints the intermediate gradient to a file.
   // fptr_startingpoint_vec = tzFopen(filename_sp_vec_minproj,"r");
   // fprintf(fptr_startingpoint_vec, "--------Numerical gradient---------\n");
   // for (_i=0; _i<n; _i++)  fprintf(fptr_startingpoint_vec, " %0.16e ", g[_i]);
   // tzFclose(fptr_startingpoint_vec);

   GLB_DISPLAY = 1;
   // printf("\n=== Leaving the gradient function. ===\n");
   // fflush(stdout);
}
//--- For conjugate gradient minimization as well.
static double ObjFuncForModel_congrad(double *x0_p, int d_x0)
{
   TSdvector x0_sdv;
   x0_sdv.v = x0_p;
   x0_sdv.n = d_x0;
   x0_sdv.flag = V_DEF;

   return ( -opt_logOverallPosteriorKernal(SMODEL_PS, &x0_sdv) );
}







//-----------------------------------------------------------------------
// Conjugate gradient method I minimization package.
//-----------------------------------------------------------------------
struct TSpackage_congrad1_tag *CreateTSpackage_congrad1(void)
{
   //===
   struct TSpackage_congrad1_tag *package_congrad1_ps = tzMalloc(1, struct TSpackage_congrad1_tag);

   package_congrad1_ps->crit = CRIT_CONGRAD1;
   package_congrad1_ps->itmax = ITMAX_CONGRAD1;

   return (package_congrad1_ps);
}
//---
struct TSpackage_congrad1_tag *DestroyTSpackage_congrad1(struct TSpackage_congrad1_tag *package_congrad1_ps)
{
   if (package_congrad1_ps)
   {
      //===
      free(package_congrad1_ps);
      return ((struct TSpackage_congrad1_tag *)NULL);
   }
   else  return (package_congrad1_ps);
}


/**
static void imslconlin_gradcd(int n, double *x, double *g) {
   //Outputs:
   //  g: the gradient n-by-1 g (no need to be initialized).
   //Inputs:
   //  x: the vector point at which the gradient is evaluated.  No change in the end although will be added or substracted by dh during the function (but in the end the original value will be put back).
   //  n: the dimension of g or x.
   int _i;

   // printf("\n=== Entering the gradient function. ===\n");
   // fflush(stdout);
   GLB_DISPLAY = 0;   //This guarantees so that the objective function printouts will not show when the gradient is computed.
   gradcd_gen(g, x, n, congrad_ObjFuncForTVBVAR, (double *)NULL, congrad_ObjFuncForTVBVAR(x, n));
   //=== Prints the intermediate gradient to a file.
   fptr_startingpoint_grad = tzFopen(filename_spgrad,"w");
   fprintf(fptr_startingpoint_grad, "--------Numerical gradient---------\n");
   for (_i=0; _i<n; _i++)  fprintf(fptr_startingpoint_grad, " %0.16e ", g[_i]);
   tzFclose(fptr_startingpoint_grad);
   GLB_DISPLAY = 1;
   // printf("\n=== Leaving the gradient function. ===\n");
   // fflush(stdout);
}
/**/



/**
   //=== Conjugate gradient method, which works too slowly but is reliable.  Thus, it is used to finish it up.
   congradmin_SetPrintFile(filename_spvec);
   frprmn(x0_dv->v, x0_dv->n, &niter, &fret, congrad_ObjFuncForTVBVAR, gradcd_gen, &ftol, &itmax, (double *)NULL, (int *)NULL, (double *)NULL);


void frprmn(double p[], int n, int *iter, double *fret,
            double (*func)(double [], int), void (*dfunc)(double [], double [], int, double (*func)(double [], int), double *, double),
         double *ftol_p, int *itmax_p, double *tol_brent_p, int *itmax_brent_p, double *grdh_p) {
   //Outputs:
   //  p[0, ..., n-1]:  the location of the minimum if it converges, which replaces the starting value.
   //  iter:  pointer to the number of iterations that were performed.
   //  fret:  pointer to the minimum value of the function.
   //Inputs:
   //  p[0, ..., n-1]:  a starting point for the minimization.
   //  n:  the dimension of p.
   //  ftol_p:  pointer to the convergence tolerance on the objective function value. Default: 1.0e-4 if NULL.
   //  itmax_p:    pointer to the maximum number of iterations in the main minimization program frprmn().  Default: 2000 if NULL.
   //  tol_brent_p:  pointer to the convergence tolerance for the line minimization in brent().  Default: 2.0e-4 if NULL.
   //  itmax_brent_p:  pointer to the maximum number of iterations for the line minimization in brent().  Default: 100 if NULL.
   //  grdh:  pointer to the user's specified step size for a numerical gradient.  If NULL, dfunc() (i.e., gradcd_gen()) will select grdh automatically.
   //  func():  the objective function.
   //  dfunc(): the gradient function computing the numerical gradient.  In the form of gradcd_gen() in cstz.c.


//------- For csminwel only. -------
typedef struct TSetc_congrad1_tag {
   //=== Optional input arguments, often NOT used, so we set to NULL at this point.
   double **args; //k-by-q.
   int *dims;   //k-by-1;
   int _k;

   //=== Mandatory input arguments.
   TSdmatrix *Hx_dm;  //n-by-n inverse Hessian.  Output as well, when csminwel is done.
   double crit;   //Overall convergence criterion for the function value.
   int itmax;  //Maximum number of iterations.
//   double grdh;  //Step size for the numerical gradient if no analytical gradient is available.

   //=== Some reported input arguments.
   double ini_h_csminwel;
   int indxnumgrad_csminwel;
   double gradstps_csminwel;


   //=== Output arguments.
   int badg;    //If (badg==0), analytical gradient is used; otherwise, numerical gradient will be produced.
   int niter;   //Number of iterations taken by csminwel.
   int fcount;   //Number of function evaluations used by csminwel.
   int retcode;  //Return code for the terminating condition.
                // 0, normal step (converged). 1, zero gradient (converged).
                // 4,2, back and forth adjustment of stepsize didn't finish.
                // 3, smallest stepsize still improves too slow. 5, largest step still improves too fast.
                // 6, no improvement found.
} TSetc_csminwel;


/**/



/**
//------- For IMSL multivariate linearly-constrained minimizaiton package only. -------
typedef struct TSetc_imslconlin_tag {
   //=== Non-trivial constraint arguments, whose arrays will point to the constraints specified in the specific project minpack->etc_project_ps.
   int nvars; //Total number of free parameters for the optimaization.
   int neqs;  //Number of linear equality constrains.  Must be no greater than ncons.
   int ncons; //Total number of linear equality and non-equality constrains (excluding simple bounds).
   double *lh_coefs_pd;  //ncons*nvars-by-1. Left-hand coefficients in the linear constrains (excluding simple bounds).
                         //lh_coefs_pd stacks the neqs rows first, followed by the inequality constraints.
                         //Set to NULL if ncons=0;
   double *rh_constraints_pd;  //ncons-by-1.  Right-hand constrains in the equality and non-equality constrains (excluding simple bounds).
                              //Set to NULL if ncons=0;

   //=== Trivial or simple bounds.
   double *lowbounds_pd;  //nvars-by-1.  Simple lower bounds.  If a component is unbounded, choose a very negative large value (e.g., -BIGREALNUMBER).
   double *upperbounds_pd;  //nvars-by-1.  Simple upper bounds.  If a component is unbounded, choose a very positive large value (e.g., BIGREALNUMBER).

   //=== Other inputs.
   int itmax;     //Maximum number of iterations.
   double crit;   //Overall convergence criterion on the first-order conditions.
} TSetc_imslconlin;


//-----------------------------------------------------------------------
// Linearly-constrained IMSL minimization package.
//-----------------------------------------------------------------------
static struct TSetc_imslconlin_tag *CreateTSetc_imslconlin(struct TSminpack_tag *minpack_psconst int nvars, const int neqs, const int ncons)
{
   //===
   struct TSetc_imslconlin *etc_imslconlin_ps = tzMalloc(1, struct TSetc_imslconlin_tag);

   if (neqs.ncons)  fn_DisplayErrors("optpackage.c/CreateTSetc_imslconlin: make sure # of equality constraints no greater than total # of constraints");


   return (etc_imslconlin_ps);
}

/**/
