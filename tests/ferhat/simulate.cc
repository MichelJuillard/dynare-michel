  ////////////////////////////////////////////////////////////////////////
  //                           simulate.c                               //
  //              simulate file designed for GNU GCC C++ compiler       //
  //                 use GCC_COMPILER option in MODEL command           //
  ////////////////////////////////////////////////////////////////////////

//#define PRINT_OUT
//#define PRINT_OUT_p
#include <stack>
#include <math.h>
#include <iostream>
#include <fstream>
#include "pctimer_h.hh"
#include "mex.h"


/*#include "simulate.hh"*/

int Per_y_, Per_u_, it_, nb_row_x, u_size, y_size, x_size, y_kmin, y_kmax, y_decal;
int periods, maxit_;
double *params, *u, *y, *x, *r, *g1, *g2;
double slowc, solve_tolf, max_res, res1, res2;
bool cvg;
pctimer_t t0, t1;

typedef struct t_table_y
{
  int index, nb;
  int *u_index, *y_index;
};

typedef struct t_table_u
{
  t_table_u* pNext;
  unsigned char type;
  int index;
  int op1, op2;
};


std::ifstream SaveCode;




void
read_file_table_u(t_table_u **table_u, t_table_u **F_table_u, t_table_u **i_table_u, t_table_u **F_i_table_u, int *nb_table_u, bool i_to_do, bool shifting, int *nb_add_u_count)
{
  char type;
  int i;
  SaveCode.read(reinterpret_cast<char *>(nb_table_u), sizeof(*nb_table_u));
#ifdef PRINT_OUT
  mexPrintf("->*nb_table_u=%d\n", *nb_table_u);
#endif
  *table_u = (t_table_u*)mxMalloc(sizeof(t_table_u) - 2 * sizeof(int));
  *F_table_u = *table_u;
   if(i_to_do)
    {
#ifdef PRINT_OUT
      mexPrintf("=>i_table\n");
#endif
      (*i_table_u) = (t_table_u*)mxMalloc(sizeof(t_table_u) - 2 * sizeof(int));
      (*F_i_table_u) = (*i_table_u);
    }
  for(i = 0;i < *nb_table_u;i++)
    {
      SaveCode.read(reinterpret_cast<char *>(&type), sizeof(type));
      switch (type)
        {
        case 3:
        case 7:
          (*table_u)->pNext = (t_table_u*)mxMalloc(sizeof(t_table_u) - sizeof(int));
          (*table_u) = (*table_u)->pNext;
          (*table_u)->type = type;
          SaveCode.read(reinterpret_cast<char *>(&(*table_u)->index), sizeof((*table_u)->index));
          SaveCode.read(reinterpret_cast<char *>(&(*table_u)->op1), sizeof((*table_u)->op1));
          if(shifting)
            {
              (*table_u)->index -= y_kmin * u_size;
              (*table_u)->op1 -= y_kmin * u_size;
            }
#ifdef PRINT_OUT

          if((*table_u)->type == 3)
            mexPrintf("u[%d]=-1/u[%d]\n", (*table_u)->index, (*table_u)->op1);
          else
            mexPrintf("u[%d]*=u[%d]\n", (*table_u)->index, (*table_u)->op1);
#endif
          if(i_to_do)
            {
              (*i_table_u)->pNext = (t_table_u*)mxMalloc(sizeof(t_table_u) - sizeof(int));
              (*i_table_u) = (*i_table_u)->pNext;
              (*i_table_u)->type = type;
              SaveCode.read(reinterpret_cast<char *>(&(*i_table_u)->index), sizeof((*i_table_u)->index));
              SaveCode.read(reinterpret_cast<char *>(&(*i_table_u)->op1), sizeof((*i_table_u)->op1));
#ifdef FIXE
              (*i_table_u)->index = u_size;
              (*i_table_u)->op1 = u_size;
#endif
#ifdef PRINT_OUT
              if((*i_table_u)->type == 3)
                mexPrintf("i u[%d]=1/(1-u[%d])\n", (*i_table_u)->index, (*i_table_u)->op1);
              else
                mexPrintf("i u[%d]*=u[%d]\n", (*i_table_u)->index, (*i_table_u)->op1);
#endif
            }
          break;
        case 1:
        case 2:
        case 6:
          (*table_u)->pNext = (t_table_u*)mxMalloc(sizeof(t_table_u));
          (*table_u) = (*table_u)->pNext;
          (*table_u)->type = type;
          SaveCode.read(reinterpret_cast<char *>(&(*table_u)->index), sizeof((*table_u)->index));
          SaveCode.read(reinterpret_cast<char *>(&(*table_u)->op1), sizeof((*table_u)->op1));
          SaveCode.read(reinterpret_cast<char *>(&(*table_u)->op2), sizeof((*table_u)->op2));
          if(shifting)
            {
              (*table_u)->index -= y_kmin * u_size;
              (*table_u)->op1 -= y_kmin * u_size;
              (*table_u)->op2 -= y_kmin * u_size;
            }
          if((*table_u)->type == 1)
            {
#ifdef PRINT_OUT
              mexPrintf("u[%d]=u[%d]*u[%d]\n", (*table_u)->index, (*table_u)->op1, (*table_u)->op2);
#endif
              if(i_to_do)
                (*nb_add_u_count)++;
            }
#ifdef PRINT_OUT
          else if((*table_u)->type == 2)
            mexPrintf("u[%d]+=u[%d]*u[%d]\n", (*table_u)->index, (*table_u)->op1, (*table_u)->op2);
          else
            mexPrintf("u[%d]=1/(1-u[%d]*u[%d])\n", (*table_u)->index, (*table_u)->op1, (*table_u)->op2);
#endif
          if(i_to_do)
            {
              (*i_table_u)->pNext = (t_table_u*)mxMalloc(sizeof(t_table_u));
              (*i_table_u) = (*i_table_u)->pNext;
              (*i_table_u)->type = type;
              SaveCode.read(reinterpret_cast<char *>(&(*i_table_u)->index), sizeof((*i_table_u)->index));
              SaveCode.read(reinterpret_cast<char *>(&(*i_table_u)->op1), sizeof((*i_table_u)->op1));
              SaveCode.read(reinterpret_cast<char *>(&(*i_table_u)->op2), sizeof((*i_table_u)->op2));
#ifdef FIXE
              (*i_table_u)->index = u_size;
              (*i_table_u)->op1 = u_size;
              (*i_table_u)->op2 = u_size;
#endif
#ifdef PRINT_OUT
              if((*i_table_u)->type == 1)
                mexPrintf("i u[%d]=u[%d]*u[%d]\n", (*i_table_u)->index, (*i_table_u)->op1, (*i_table_u)->op2);
              else if((*i_table_u)->type == 2)
                mexPrintf("i u[%d]+=u[%d]*u[%d]\n", (*i_table_u)->index, (*i_table_u)->op1, (*i_table_u)->op2);
              else
                mexPrintf("i u[%d]=1/(1-u[%d]*u[%d])\n", (*i_table_u)->index, (*i_table_u)->op1, (*i_table_u)->op2);
#endif
            }
          break;
        case 5:
          (*table_u)->pNext = (t_table_u*)mxMalloc(sizeof(t_table_u) - 2 * sizeof(int));
          (*table_u) = (*table_u)->pNext;
          (*table_u)->type = type;
          SaveCode.read(reinterpret_cast<char *>(&(*table_u)->index), sizeof((*table_u)->index));
          if(shifting)
            (*table_u)->index -= y_kmin * u_size;
#ifdef PRINT_OUT
          mexPrintf("push(u[%d])\n", (*table_u)->index);
#endif
          if(i_to_do)
            {
              (*i_table_u)->pNext = (t_table_u*)mxMalloc(sizeof(t_table_u) - 2 * sizeof(int));
              (*i_table_u) = (*i_table_u)->pNext;
              (*i_table_u)->type = type;
              SaveCode.read(reinterpret_cast<char *>(&(*i_table_u)->index), sizeof((*i_table_u)->index));
#ifdef FIXE
              (*i_table_u)->index = u_size;
#endif
#ifdef PRINT_OUT
              mexPrintf("i push(u[%d])\n", (*i_table_u)->index);
#endif
            }
          break;
        }
    }
}

void
read_file_table_y(t_table_y **table_y, t_table_y **i_table_y, int *nb_table_y, bool i_to_do, bool shifting)
{
  int i, k;
  SaveCode.read(reinterpret_cast<char *>(nb_table_y), sizeof(*nb_table_y));
#ifdef PRINT_OUT
  mexPrintf("nb_table_y=%d\n", *nb_table_y);
  mexPrintf("y_size=%d, u_size=%d, y_kmin=%d, y_kmax=%d\n", y_size, u_size, y_kmin, y_kmax);
#endif
  (*table_y) = (t_table_y*)mxMalloc((*nb_table_y) * sizeof(t_table_y));
  for(i = 0;i < *nb_table_y;i++)
    {
      SaveCode.read(reinterpret_cast<char *>(&((*table_y)[i].index)), sizeof((*table_y)[i].index));
      SaveCode.read(reinterpret_cast<char *>(&((*table_y)[i].nb)), sizeof((*table_y)[i].nb));
      if(shifting)
        (*table_y)[i].index -= y_kmin * y_size;
#ifdef PRINT_OUT
      mexPrintf("table_y[i].nb=%d\n", (*table_y)[i].nb);
      mexPrintf("y[%d]=", (*table_y)[i].index);
#endif
      (*table_y)[i].u_index = (int*)mxMalloc((*table_y)[i].nb * sizeof(int));
      (*table_y)[i].y_index = (int*)mxMalloc((*table_y)[i].nb * sizeof(int));
      for(k = 0;k < (*table_y)[i].nb;k++)
        {
          SaveCode.read(reinterpret_cast<char *>(&((*table_y)[i].u_index[k])), sizeof((*table_y)[i].u_index[k]));
          SaveCode.read(reinterpret_cast<char *>(&((*table_y)[i].y_index[k])), sizeof((*table_y)[i].y_index[k]));
          if(shifting)
            {
              (*table_y)[i].u_index[k] -= y_kmin * u_size;
              if(((*table_y)[i].y_index[k] > y_size*y_kmin) && ((*table_y)[i].y_index[k] < y_size*(2*y_kmin + y_kmax + 2)))
                {
                  (*table_y)[i].y_index[k] -= y_kmin * y_size;
                }
            }
#ifdef PRINT_OUT
          if(k < (*table_y)[i].nb - 1)
            mexPrintf("u[%d]*y[%d]+", (*table_y)[i].u_index[k], (*table_y)[i].y_index[k]);
          else
            mexPrintf("u[%d]*y[%d]\n", (*table_y)[i].u_index[k], (*table_y)[i].y_index[k]);
#endif
        }
    }
#ifdef PRINT_OUT
  mexPrintf("*nb_table_y=%d\n", *nb_table_y);
#endif
  if(i_to_do)
    {
      *i_table_y = (t_table_y*)mxMalloc((*nb_table_y) * sizeof(t_table_y));
      for(i = 0;i < *nb_table_y;i++)
        {
          SaveCode.read(reinterpret_cast<char *>(&((*i_table_y)[i].index)), sizeof((*i_table_y)[i].index));
          SaveCode.read(reinterpret_cast<char *>(&((*i_table_y)[i].nb)), sizeof((*i_table_y)[i].nb));
#ifdef PRINT_OUT
          mexPrintf("(*i_table_y)[i].nb=%d\n", (*i_table_y)[i].nb);
          mexPrintf("y_i[%d]=", (*i_table_y)[i].index);
#endif
          (*i_table_y)[i].u_index = (int*)mxMalloc((*i_table_y)[i].nb * sizeof(int));
          (*i_table_y)[i].y_index = (int*)mxMalloc((*i_table_y)[i].nb * sizeof(int));
          for(k = 0;k < (*i_table_y)[i].nb;k++)
            {
              SaveCode.read(reinterpret_cast<char *>(&((*i_table_y)[i].u_index[k])), sizeof((*i_table_y)[i].u_index[k]));
              SaveCode.read(reinterpret_cast<char *>(&((*i_table_y)[i].y_index[k])), sizeof((*i_table_y)[i].y_index[k]));
#ifdef PRINT_OUT
              if(k < (*i_table_y)[i].nb - 1)
                mexPrintf("u[%d]*y[%d]+", (*i_table_y)[i].u_index[k], (*i_table_y)[i].y_index[k]);
              else
                mexPrintf("u[%d]*y[%d]\n", (*i_table_y)[i].u_index[k], (*i_table_y)[i].y_index[k]);
#endif
            }
        }
    }
}


int i, j, k, nb_endo, u_count, u_count_init, iter;
int nb_prologue_table_u, nb_first_table_u, nb_middle_table_u, nb_last_table_u;
int nb_prologue_table_y, nb_first_table_y, nb_middle_table_y, nb_last_table_y;
int first_count_loop, middle_count_loop;
char type;
t_table_u *prologue_table_u, *first_table_u, *first_i_table_u, *middle_table_u, *middle_i_table_u, *last_table_u;
t_table_y *prologue_table_y, *first_table_y, *middle_table_y, *middle_i_table_y, *last_table_y;
t_table_u *F_prologue_table_u, *F_first_table_u, *F_first_i_table_u, *F_middle_table_u, *F_middle_i_table_u, *F_last_table_u;
std::string filename;


void
Read_file(std::string file_name, int periods, int u_size1, int y_size, int y_kmin, int y_kmax)
{
  int nb_add_u_count = 0;
  u_size = u_size1;
  filename=file_name;
#ifdef PRINT_OUT
  mexPrintf("%s\n", file_name.c_str());
#endif
  if(!SaveCode.is_open())
    {
#ifdef PRINT_OUT
      mexPrintf("file opened\n");
#endif
      SaveCode.open((file_name + ".bin").c_str(), std::ios::in | std::ios::binary);
      if(!SaveCode.is_open())
        {
          mexPrintf("Error : Can't open file \"%s\" for reading\n", (file_name + ".bin").c_str());
          mexErrMsgTxt("Exit from Dynare");
        }
#ifdef PRINT_OUT
      mexPrintf("done\n");
#endif
    }
  SaveCode.read(reinterpret_cast<char *>(&nb_endo), sizeof(nb_endo));
  SaveCode.read(reinterpret_cast<char *>(&u_count), sizeof(u_count));
  SaveCode.read(reinterpret_cast<char *>(&u_count_init), sizeof(u_count_init));
#ifdef PRINT_OUT
  mexPrintf("nb_endo=%d\n", nb_endo);
  mexPrintf("u_count=%d\n", u_count);
  mexPrintf("u_count_init=%d\n", u_count_init);
  //mexPrintf("first table_u\n");
#endif
  read_file_table_u(&first_table_u, &F_first_table_u, &first_i_table_u, &F_first_i_table_u, &nb_first_table_u, true, false, &nb_add_u_count);
#ifdef PRINT_OUT
  mexPrintf("nb_first_table_u=%d\n", nb_first_table_u);
#endif
//mexErrMsgTxt("Exit from Dynare");
#ifdef PRINT_OUT
  mexPrintf("prologue table_u\n");
#endif
  read_file_table_u(&prologue_table_u, &F_prologue_table_u, NULL, NULL, &nb_prologue_table_u, false, false, &nb_add_u_count);
#ifdef PRINT_OUT
  mexPrintf("nb_prologue_table_u=%d\n", nb_prologue_table_u);
#endif
  //mexErrMsgTxt("Exit from Dynare");
  SaveCode.read(reinterpret_cast<char *>(&middle_count_loop), sizeof(middle_count_loop));
#ifdef PRINT_OUT
  mexPrintf("middle_count_loop=%d\n",middle_count_loop);
#endif
  //mexErrMsgTxt("Exit from Dynare");
#ifdef PRINT_OUT
  mexPrintf("middle table_u\n");
#endif
  read_file_table_u(&middle_table_u, &F_middle_table_u, &middle_i_table_u, &F_middle_i_table_u, &nb_middle_table_u, true,  /*true*/false, &nb_add_u_count);
#ifdef PRINT_OUT
 mexPrintf("nb_middle_table_u=%d\n",nb_middle_table_u);
  //mexPrintf("last table_u\n");
#endif
  read_file_table_u(&last_table_u, &F_last_table_u, NULL, NULL, &nb_last_table_u, false, false, &nb_add_u_count);
#ifdef PRINT_OUT
  mexPrintf("->nb_last_table_u=%d\n", nb_last_table_u);
  mexPrintf("i=%d\n", i);
  mexPrintf("going to read prologue_table_y\n");
#endif
  read_file_table_y(&prologue_table_y, NULL, &nb_prologue_table_y, false, false);
#ifdef PRINT_OUT
  mexPrintf("nb_prologue_table_y=%d\n", nb_prologue_table_y);
  mexPrintf("going to read first_table_y\n");
#endif
  read_file_table_y(&first_table_y, NULL, &nb_first_table_y, false, false);
#ifdef PRINT_OUT
  mexPrintf("nb_first_table_y=%d\n", nb_first_table_y);
  mexPrintf("going to read middle_table_y\n");
#endif
  read_file_table_y(&middle_table_y, &middle_i_table_y, &nb_middle_table_y, true,  /*true*/false);
#ifdef PRINT_OUT
  mexPrintf("nb_middle_table_y=%d\n", nb_middle_table_y);
  mexPrintf("going to read last_table_y\n");
#endif
  read_file_table_y(&last_table_y, NULL, &nb_last_table_y, false, false);
#ifdef PRINT_OUT
  mexPrintf("nb_last_table_y=%d\n", nb_last_table_y);
  mexPrintf("->nb_last_table_y=%d\n", nb_last_table_y);
#endif
  if(nb_last_table_u > 0)
    {
#ifdef PRINT_OUT
      mexPrintf("y_size=%d, periods=%d, y_kmin=%d, y_kmax=%d\n", y_size, periods, y_kmin, y_kmax);
      mexPrintf("u=mxMalloc(%d)\n", u_count + 1);
#endif
      u = (double*)mxMalloc((u_count + 1) * sizeof(double));
    }
  else
    {
#ifdef PRINT_OUT
      mexPrintf("u_size=%d, y_size=%d, periods=%d, y_kmin=%d, y_kmax=%d, u_count=%d, nb_add_u_count=%d\n", u_size, y_size, periods, y_kmin, y_kmax, u_count, nb_add_u_count);
      mexPrintf("u=mxMalloc(%d)\n", u_count + (periods + y_kmin + y_kmax)* /*(u_count-u_size*(periods+y_kmin+y_kmax))*/nb_add_u_count);
#endif
      u = (double*)mxMalloc((u_count + (periods + y_kmin + y_kmax)* /*(u_count-u_size*(periods+y_kmin+y_kmax)*/nb_add_u_count) * sizeof(double));
      memset(u, 0, (u_count + (periods + y_kmin + y_kmax)* /*(u_count-u_size*(periods+y_kmin+y_kmax)*/nb_add_u_count)*sizeof(double));
    }
  if(u == NULL)
    {
      mexPrintf("memory exhausted\n");
      mexErrMsgTxt("Exit from Dynare");
    }
  // mexErrMsgTxt("Exit from Dynare");
}


std::stack <double> Stack;

void
simulate(int blck, int y_size, int it_, int y_kmin, int y_kmax)
{
  int i, j, k, l, m, m1, nop;
  int period = it_ * y_size, s_middle_count_loop = 0 ;
  pctimer_t t1 = pctimer();
  double uu, yy;
  char tmp_s[150];
#ifdef PRINT_OUT
  for(j = 0;j < it_ -y_kmin;j++)
    {
      for(i = 0;i < u_size;i++)
        {
          mexPrintf("u[%d]=%f ", j*u_size + i, u[j*u_size + i]);
        }
      mexPrintf("\n");
    }
#endif
  if(nb_first_table_u > 0)
    {
      first_count_loop = it_;
      s_middle_count_loop = it_ -y_kmin - middle_count_loop + 1;
//#ifdef PRINT_OUT
      mexPrintf("----------------------------------------------------------------------\n");
      mexPrintf("      Simulate     itération° %d     \n",iter+1);
      mexPrintf("      max. error=%.10e       \n",max_res);
      mexPrintf("      sqr. error=%.10e       \n",res2);
      mexPrintf("      abs. error=%.10e       \n",res1);
      mexPrintf("----------------------------------------------------------------------\n");
//#endif
    }
  nop = 0;
  for(j = 0 ;j < first_count_loop - y_kmin;j++)
    {
      first_table_u = F_first_table_u->pNext;
      first_i_table_u = F_first_i_table_u->pNext;
      for(i = 0;i < nb_first_table_u;i++)
        {
          switch (first_table_u->type)
            {
            case 1:
              u[first_table_u->index + j*first_i_table_u->index] = u[first_table_u->op1 + j * first_i_table_u->op1] * u[first_table_u->op2 + j * first_i_table_u->op2];
#ifdef PRINT_OUT
              mexPrintf("u[%d]=u[%d]*u[%d]=%f\n", first_table_u->index + j*first_i_table_u->index , first_table_u->op1 + j*first_i_table_u->op1, first_table_u->op2 + j*first_i_table_u->op2, u[first_table_u->index + j*first_i_table_u->index]);
#endif
              break;
            case 2:
              u[first_table_u->index + j*first_i_table_u->index] += u[first_table_u->op1 + j * first_i_table_u->op1] * u[first_table_u->op2 + j * first_i_table_u->op2];
#ifdef PRINT_OUT
              mexPrintf("u[%d]+=u[%d]*u[%d]=%f\n" , first_table_u->index + j*first_i_table_u->index, first_table_u->op1 + j*first_i_table_u->op1, first_table_u->op2 + j*first_i_table_u->op2, u[first_table_u->index + j*first_i_table_u->index]);
#endif
              break;
            case 3:
              u[first_table_u->index + j*first_i_table_u->index] = 1 / ( -u[first_table_u->op1 + j * first_i_table_u->op1]);
#ifdef PRINT_OUT
              mexPrintf("u[%d]=1/(-u[%d])=%f\n", first_table_u->index + j*first_i_table_u->index, first_table_u->op1 + j*first_i_table_u->op1, u[first_table_u->index + j*first_i_table_u->index]);
#endif
              break;
            case 5:
              Stack.push(u[first_table_u->index + j*first_i_table_u->index]);
#ifdef PRINT_OUT
              mexPrintf("push(u[%d])\n", first_table_u->index + j*first_i_table_u->index);
#endif
              break;
            case 6:
              u[first_table_u->index + j*first_i_table_u->index] = 1 / (1 - u[first_table_u->op1 + j * first_i_table_u->op1] * u[first_table_u->op2 + j * first_i_table_u->op2]);
#ifdef PRINT_OUT
              mexPrintf("u[%d]=1/(1-u[%d]*u[%d])=%f\n", first_table_u->index + j*first_i_table_u->index, first_table_u->op1 + j*first_i_table_u->op1, first_table_u->op2 + j*first_i_table_u->op2, u[first_table_u->index + j*first_i_table_u->index]);
#endif
              break;
            case 7:
              u[first_table_u->index + j*first_i_table_u->index] *= u[first_table_u->op1 + j * first_i_table_u->op1];
#ifdef PRINT_OUT
              mexPrintf("u[%d]*=u[%d]=%f\n", first_table_u->index + j*first_i_table_u->index, first_table_u->op1 + j*first_i_table_u->op1, u[first_table_u->index + j*first_i_table_u->index]);
#endif
              break;
            }
          if(isnan(u[first_table_u->index+ j*first_i_table_u->index]) || isinf(u[first_table_u->index+ j*first_i_table_u->index]))
           {
             mexPrintf("Error during the computation of u[%d] at time %d (in first_table_u) (operation type %d)",first_table_u->index,j,int(first_table_u->type));
             filename+=" stopped";
             mexErrMsgTxt(filename.c_str());
           }
          first_table_u = first_table_u->pNext;
          first_i_table_u = first_i_table_u->pNext;
          nop++;
        }
    }
#ifdef PRINT_OUT_p
  mexPrintf("prologue\n");
#endif
  //int nb_prologue_push=0;
  prologue_table_u = F_prologue_table_u->pNext;
  for(i = 0;i < nb_prologue_table_u;i++)
    {
      switch (prologue_table_u->type)
        {
        case 1:
          u[prologue_table_u->index ] = u[prologue_table_u->op1 ] * u[prologue_table_u->op2 ];
#ifdef PRINT_OUT_p
          mexPrintf("u[%d]=u[%d]*u[%d]=%f\n", prologue_table_u->index , prologue_table_u->op1 , prologue_table_u->op2 , u[prologue_table_u->index ]);
#endif
          break;
        case 2:
          u[prologue_table_u->index ] += u[prologue_table_u->op1 ] * u[prologue_table_u->op2 ];
#ifdef PRINT_OUT_p
          mexPrintf("u[%d]+=u[%d]*u[%d]=%f\n" , prologue_table_u->index , prologue_table_u->op1 , prologue_table_u->op2 , u[prologue_table_u->index ]);
#endif
          break;
        case 3:
          u[prologue_table_u->index ] = 1 / ( -u[prologue_table_u->op1 ]);
#ifdef PRINT_OUT_p
          mexPrintf("u[%d]=1/(-u[%d])=%f\n", prologue_table_u->index, prologue_table_u->op1, u[prologue_table_u->index]);
#endif
          break;
        case 5:
          //nb_prologue_push++;
          Stack.push(u[prologue_table_u->index]);
#ifdef PRINT_OUT_p
          mexPrintf("push(u[%d])\n", prologue_table_u->index );
#endif
          break;
        case 6:
          u[prologue_table_u->index ] = 1 / (1 - u[prologue_table_u->op1] * u[prologue_table_u->op2]);
#ifdef PRINT_OUT_p
          mexPrintf("u[%d]=1/(1-u[%d]*u[%d])=%f\n", prologue_table_u->index, prologue_table_u->op1, prologue_table_u->op2, u[prologue_table_u->index]);
#endif
          break;
        case 7:
          u[prologue_table_u->index] *= u[prologue_table_u->op1];
#ifdef PRINT_OUT_p
          mexPrintf("u[%d]*=u[%d]=%f\n", prologue_table_u->index, prologue_table_u->op1, u[prologue_table_u->index]);
#endif
          break;
        }
      if(isnan(u[prologue_table_u->index]) || isinf(u[prologue_table_u->index]))
        {
          mexPrintf("Error during the computation of u[%d] (in prologue_table_u)",prologue_table_u->index);
          filename+=" stopped";
          mexErrMsgTxt(filename.c_str());
        }
      prologue_table_u = prologue_table_u->pNext;
      nop++;
    }
#ifdef PRINT_OUT
  mexPrintf("middle_u (s_middle_count_loop=%d\n", s_middle_count_loop);
#endif
  //int nb_middle_push=0;
  for(j = 0;j < s_middle_count_loop - y_kmin;j++)
    {
      //cout << "j=" << j << "\n";
#ifdef PRINT_OUT
      mexPrintf("-----------------------------------------------------------------\n");
#endif
      middle_table_u = F_middle_table_u->pNext;
      middle_i_table_u = F_middle_i_table_u->pNext;
      for(i = 0;i < nb_middle_table_u;i++)
        {
          switch (middle_table_u->type)
            {
            case 1:
              u[middle_table_u->index + j*middle_i_table_u->index] = u[middle_table_u->op1 + j * middle_i_table_u->op1] * u[middle_table_u->op2 + j * middle_i_table_u->op2];
#ifdef PRINT_OUT
              mexPrintf("u[%d+%d*%d=%d]=u[%d]*u[%d]=%f\n", middle_table_u->index, j, middle_i_table_u->index, middle_table_u->index + j*middle_i_table_u->index, middle_table_u->op1 + j*middle_i_table_u->op1, middle_table_u->op2 + j*middle_i_table_u->op2, u[middle_table_u->index + j*middle_i_table_u->index]);
#endif
              break;
            case 2:
              u[middle_table_u->index + j*middle_i_table_u->index] += u[middle_table_u->op1 + j * middle_i_table_u->op1] * u[middle_table_u->op2 + j * middle_i_table_u->op2];
#ifdef PRINT_OUT
              mexPrintf("u[%d+%d*%d=%d]+=u[%d]*u[%d]=%f\n" , middle_table_u->index, j, middle_i_table_u->index , middle_table_u->index + j*middle_i_table_u->index , middle_table_u->op1 + j*middle_i_table_u->op1, middle_table_u->op2 + j*middle_i_table_u->op2, u[middle_table_u->index + j*middle_i_table_u->index]);
#endif
              break;
            case 3:
              u[middle_table_u->index + middle_i_table_u->index] = 1 / ( -u[middle_table_u->op1 + j * middle_i_table_u->op1]);
#ifdef PRINT_OUT
              mexPrintf("u[%d+%d*%d=%d]=1/(-u[%d])=%f\n", middle_table_u->index, j, middle_i_table_u->index, middle_table_u->index + j*middle_i_table_u->index, middle_table_u->op1 + j*middle_i_table_u->op1, u[middle_table_u->index + j*middle_i_table_u->index]);
#endif
              break;
            case 5:
#ifdef PRINT_OUT
              mexPrintf("push(u[%d+%d*%d=%d])\n", middle_table_u->index, j, middle_i_table_u->index, middle_table_u->index + j*middle_i_table_u->index);
#endif
              //nb_middle_push++;
              Stack.push(u[middle_table_u->index + j*middle_i_table_u->index]);
              break;
            case 6:
              u[middle_table_u->index + j*middle_i_table_u->index] = 1 / (1 - u[middle_table_u->op1 + j * middle_i_table_u->op1] * u[middle_table_u->op2 + j * middle_i_table_u->op2]);
#ifdef PRINT_OUT
              mexPrintf("u[%d+%d*%d=%d]=1/(1-u[%d]*u[%d])=%f\n", middle_table_u->index, j, middle_i_table_u->index, middle_table_u->index + j*middle_i_table_u->index, middle_table_u->op1 + j*middle_i_table_u->op1, middle_table_u->op2 + j*middle_i_table_u->op2, u[middle_table_u->index + j*middle_i_table_u->index]);
#endif
              break;
            case 7:
              u[middle_table_u->index + j*middle_i_table_u->index] *= u[middle_table_u->op1 + j * middle_i_table_u->op1];
#ifdef PRINT_OUT
              mexPrintf("u[%d+%d*%d=%d]*=u[%d]=%f\n", middle_table_u->index, j, middle_i_table_u->index, middle_table_u->index + j*middle_i_table_u->index, middle_table_u->op1 + j*middle_i_table_u->op1, u[middle_table_u->index + j*middle_i_table_u->index]);
#endif
              break;
            }
          if(isnan(u[middle_table_u->index+ j*middle_i_table_u->index]) || isinf(u[middle_table_u->index+ j*middle_i_table_u->index]))
           {
             mexPrintf("Error during the computation of u[%d] at time %d (in middle_table_u)",middle_table_u->index,j);
             filename+=" stopped";
             mexErrMsgTxt(filename.c_str());
           }
          middle_table_u = middle_table_u->pNext;
          middle_i_table_u = middle_i_table_u->pNext;
          nop++;
        }
    }
#ifdef PRINT_OUT
  mexPrintf("last_u\n");
#endif
  last_table_u = F_last_table_u->pNext;
  for(i = 0;i < nb_last_table_u ;i++)
    {
      switch (last_table_u->type)
        {
        case 1:
          u[last_table_u->index] = u[last_table_u->op1] * u[last_table_u->op2];
#ifdef PRINT_OUT
          mexPrintf("u[%d]=u[%d]*u[%d]=%f\n", last_table_u->index, last_table_u->op1, last_table_u->op2, u[last_table_u->index]);
#endif
          break;
        case 2:
          u[last_table_u->index] += u[last_table_u->op1] * u[last_table_u->op2];
#ifdef PRINT_OUT
          mexPrintf("u[%d]+=u[%d]*u[%d]=%f\n", last_table_u->index, last_table_u->op1, last_table_u->op2, u[last_table_u->index]);
#endif
          break;
        case 3:
          u[last_table_u->index] = 1 / ( -u[last_table_u->op1]);
#ifdef PRINT_OUT
          mexPrintf("u[%d]=1/(-u[%d])=%f\n", last_table_u->index, last_table_u->op1, u[last_table_u->index]);
#endif
          break;
        case 5:
          Stack.push(u[last_table_u->index]);
#ifdef PRINT_OUT
          mexPrintf("push(u[%d])\n", last_table_u->index);
#endif
          break;
        case 6:
          u[last_table_u->index] = 1 / (1 - u[last_table_u->op1] * u[last_table_u->op2]);
#ifdef PRINT_OUT
          mexPrintf("u[%d]=1/(1-u[%d]*u[%d])=%f\n", last_table_u->index, last_table_u->op1, last_table_u->op2, u[last_table_u->index]);
#endif
          break;
        case 7:
          u[last_table_u->index] *= u[last_table_u->op1];
#ifdef PRINT_OUT
          mexPrintf("u[%d]*=u[%d]=%f\n", last_table_u->index, last_table_u->op1, u[last_table_u->index]);
#endif
          break;
        }
      if(isnan(u[last_table_u->index]) || isinf(u[last_table_u->index]))
        {
          mexPrintf("Error during the computation of u[%d] (in last_table_u)",last_table_u->index);
          filename+=" stopped";
          mexErrMsgTxt(filename.c_str());
        }
      last_table_u = last_table_u->pNext;
      nop++;
    }
  for(i = nb_last_table_y - 1;i >= 0;i--)
    {
      k = last_table_y[i].index;
      yy = 0;
      /*y[period + k] = 0;*/
//#ifdef PRINT_OUT
      mexPrintf("it_=%d\n", it_);
      mexPrintf("->y[it_*y_size+%d]=y[%d]=", k, it_*y_size + k);
//#endif
      for(j = last_table_y[i].nb - 1;j >= 0;j--)
        {
          uu = Stack.top();
          Stack.pop();
          m = last_table_y[i].y_index[j];
#ifdef PRINT_OUT
          if(j > 0)
            {
              if(m >= 0)
                mexPrintf("u[%d](%f)*y[%d]+", last_table_y[i].u_index[j], uu, last_table_y[i].y_index[j]);
              else
                mexPrintf("u[%d](%f)+", last_table_y[i].u_index[j], uu);
            }
          else
            {
              if(m >= 0)
                mexPrintf("u[%d](%f)*y[%d]", last_table_y[i].u_index[j], uu, last_table_y[i].y_index[j]);
              else
                mexPrintf("u[%d](%f)", last_table_y[i].u_index[j], uu);
            }
#endif
          if(m >= 0)
            yy/*y[period + k]*/ += uu * y[period + m];
          else
            yy/*y[period + k]*/ += uu;
        }
       if(isnan(yy) || isinf(yy))
         {
           mexPrintf("Error during the computation of y[%d] at time %d (in last_table_u)",k,period);
           filename+=" stopped";
           mexErrMsgTxt(filename.c_str());
         }
       /*if(((k-73) % y_size)==0)
         mexPrintf("y[it_*y_size +73]=%f \n",yy);*/
       y[period + k] += slowc*(yy - y[period + k]);
//#ifdef PRINT_OUT
      mexPrintf("=%f\n", y[period + k]);
//#endif
      nop++;
    }
  int nb_middle_pop=0;
  for(j = s_middle_count_loop - y_kmin - 1+y_decal;j >= y_decal;j--)
    {
      for(i = nb_middle_table_y - 1;i >= 0;i--)
        {
          k = middle_table_y[i].index + j * middle_i_table_y[i].index;
          yy = 0;
//#ifdef PRINT_OUT
          mexPrintf("(0)y[%d]=", k );
//#endif
          for(l = middle_table_y[i].nb - 1;l >= 0;l--)
            {
              uu = Stack.top();
              //mexPrintf("{");
              Stack.pop();
              //mexPrintf("}");
              //nb_middle_pop++;
              m = middle_table_y[i].y_index[l] + j * middle_i_table_y[i].y_index[l];
//#ifdef PRINT_OUT
              if(m >= 0)
                {
                  m1 = middle_table_y[i].u_index[l] + j * middle_i_table_y[i].u_index[l];
                  if(l > 0)
                    mexPrintf("u[%d](%f)*y[%d](%f)+", m1, uu, m, y[m]);
                  else
                    mexPrintf("u[%d](%f)*y[%d](%f)", m1, uu, m, y[m]);
                }
              else
                {
                  m1 = middle_table_y[i].u_index[l] + j * middle_i_table_y[i].u_index[l];
                  if(l > 0)
                    mexPrintf("u[%d](%f)*y[%d](%f)+", m1, uu, m, 1.0);
                  else
                    mexPrintf("u[%d](%f)*y[%d](%f)", m1, uu, m, 1.0);
                }
//#endif
              if(m >= 0)
                yy += uu * y[m];
              else
                yy += uu;
            }
           //mexPrintf("y[%d]=%f\n",k,yy);
           if(isnan(yy) || isinf(yy))
             {
               mexPrintf("Error during the computation of y[%d] at time %d (in middle_table_u)",middle_table_y[i].index % y_size,j+middle_table_y[i].index / y_size);
               filename+=" stopped";
               mexErrMsgTxt(filename.c_str());
             }
           /*if(((k-73) % y_size)==0)
             mexPrintf("y[it_*y_size +73]=%f \n",yy);*/
           y[k] += slowc*(yy - y[k]);
//#ifdef PRINT_OUT
          mexPrintf("=%f\n", y[k]);
//#endif
          nop++;
        }
    }
  //mexPrintf("nb_middle_push=%d et nb_middle_pop=%d\n",nb_middle_push,nb_middle_pop);
  //int nb_prologue_pop=0;
  for(i = nb_prologue_table_y - 1;i >= 0;i--)
    {
      k = prologue_table_y[i].index;
      yy = 0;
//#ifdef PRINT_OUT_p
      mexPrintf("(1)y[%d]=", k+y_decal*y_size);
//#endif
      for(j = prologue_table_y[i].nb - 1;j >= 0;j--)
        {
          //nb_prologue_pop++;
          uu = Stack.top();
          Stack.pop();
//#ifdef PRINT_OUT_p
          if(prologue_table_y[i].y_index[j] >= 0)
            {
              if(j > 0)
                mexPrintf("u[%d](%f)*y[%d](%f)+", prologue_table_y[i].u_index[j], uu, prologue_table_y[i].y_index[j], y[prologue_table_y[i].y_index[j]]);
              else
                mexPrintf("u[%d](%f)*y[%d](%f)", prologue_table_y[i].u_index[j], uu, prologue_table_y[i].y_index[j], y[prologue_table_y[i].y_index[j]]);
            }
          else
            {
              if(j > 0)
                mexPrintf("u[%d](%f)*y[%d](%f)+", prologue_table_y[i].u_index[j], uu, prologue_table_y[i].y_index[j], 1.0);
              else
                mexPrintf("u[%d](%f)*y[%d](%f)", prologue_table_y[i].u_index[j], uu, prologue_table_y[i].y_index[j], 1.0);
            }
//#endif
          if(prologue_table_y[i].y_index[j] >= 0)
            yy += uu * y[prologue_table_y[i].y_index[j]+y_decal*y_size];
          else
            yy += uu;
        }
       if(isnan(yy) || isinf(yy))
         {
           mexPrintf("Error during the computation of y[%d] at time %d (in prologue_table_u)",k % y_size, k / y_size);
           filename+=" stopped";
           mexErrMsgTxt(filename.c_str());
         }
       /*if(((k-73) % y_size)==0)
         mexPrintf("y[it_*y_size +73]=%f \n",yy);*/
       y[k+y_decal*y_size] += slowc*(yy - y[k+y_decal*y_size]);
//#ifdef PRINT_OUT_p
      mexPrintf("=%f\n", y[k+y_decal*y_size]);
//#endif
      nop++;
    }
  //mexPrintf("nb_prologue_push=%d et nb_prologue_pop=%d\n",nb_prologue_push,nb_prologue_pop);
  for(i = nb_first_table_y - 1;i >= 0;i--)
    {
      k = first_table_y[i].index;
      yy = 0;
//#ifdef PRINT_OUT_p
      mexPrintf("(2)y[%d]=", k);
//#endif
      for(j = first_table_y[i].nb - 1;j >= 0;j--)
        {
          uu = Stack.top();
          Stack.pop();
#ifdef PRINT_OUT_p
          if(j > 0)
            mexPrintf("u[%d](%f)*y[%d](%f)+", first_table_y[i].u_index[j], uu, first_table_y[i].y_index[j], y[first_table_y[i].y_index[j]]);
          else
            mexPrintf("u[%d](%f)*y[%d](%f)", first_table_y[i].u_index[j], uu, first_table_y[i].y_index[j], y[first_table_y[i].y_index[j]]);
#endif
          if(m >= 0)
            yy += uu * y[first_table_y[i].y_index[j]];
          else
            yy += uu;
        }
      if(isnan(yy) || isinf(yy))
         {
           mexPrintf("Error during the computation of y[%d] (in first_table_u)",k);
           filename+=" stopped";
           mexErrMsgTxt(filename.c_str());
         }
      /*if(((k-73) % y_size)==0)
         mexPrintf("y[it_*y_size +73]=%f \n",yy);*/
      y[k] += slowc*(yy - y[k]);
//#ifdef PRINT_OUT_p
      mexPrintf("=%f\n", y[k]);
//#endif
      nop++;
    }
  pctimer_t t2 = pctimer();
  if(nb_first_table_u > 0)
    {
      mexPrintf("(**%f milliseconds u_count : %d  nop : %d **)\n", 1000*(t2 - t1), u_count, nop);
      mexEvalString("drawnow;");
    }
  /*if((nb_last_table_u>0)&&(it_>y_kmin))
    mexErrMsgTxt("Exit from Dynare");*/
}

