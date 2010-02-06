/*
 * Copyright (C) 2007-2009 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

//#define _GLIBCXX_USE_C99_FENV_TR1 1
//#include <cfenv>

#include <cstring>
#include <ctime>
#include <sstream>
#include "SparseMatrix.hh"

#ifdef _MSC_VER
unsigned long _nan[2] = { 0xffffffff, 0x7fffffff };
double NAN = *((double *) _nan);
#endif

SparseMatrix::SparseMatrix()
{
  pivotva = NULL;
  g_save_op = NULL;
  g_nop_all = 0;
  mem_mngr.init_Mem();
  symbolic = true;
  alt_symbolic = false;
  alt_symbolic_count = 0;
  max_u = 0;
  min_u = 0x7FFFFFFF;
  res1a = 9.0e60;
  tbreak_g = 0;
  start_compare = 0;
  restart = 0;
  IM_i.clear();
}

int
SparseMatrix::NRow(int r)
{
  return NbNZRow[r];
}

int
SparseMatrix::NCol(int c)
{
  return NbNZCol[c];
}

int
SparseMatrix::At_Row(int r, NonZeroElem **first)
{
  (*first) = FNZE_R[r];
  return NbNZRow[r];
}

int
SparseMatrix::Union_Row(int row1, int row2)
{
  NonZeroElem *first1, *first2;
  int n1 = At_Row(row1, &first1);
  int n2 = At_Row(row2, &first2);
  int i1 = 0, i2 = 0, nb_elem = 0;
  while (i1 < n1 && i2 < n2)
    {
      if (first1->c_index == first2->c_index)
        {
          nb_elem++;
          i1++;
          i2++;
          first1 = first1->NZE_R_N;
          first2 = first2->NZE_R_N;
        }
      else if (first1->c_index < first2->c_index)
        {
          nb_elem++;
          i1++;
          first1 = first1->NZE_R_N;
        }
      else
        {
          nb_elem++;
          i2++;
          first2 = first2->NZE_R_N;
        }
    }
  return nb_elem;
}

int
SparseMatrix::At_Pos(int r, int c, NonZeroElem **first)
{
  (*first) = FNZE_R[r];
  while ((*first)->c_index != c)
    (*first) = (*first)->NZE_R_N;
  return NbNZRow[r];
}

int
SparseMatrix::At_Col(int c, NonZeroElem **first)
{
  (*first) = FNZE_C[c];
  return NbNZCol[c];
}

int
SparseMatrix::At_Col(int c, int lag, NonZeroElem **first)
{
  (*first) = FNZE_C[c];
  int i = 0;
  while ((*first)->lag_index != lag && (*first))
    (*first) = (*first)->NZE_C_N;
  if ((*first))
    {
      NonZeroElem *firsta = (*first);
      if (!firsta->NZE_C_N)
        i++;
      else
        {
          while (firsta->lag_index == lag && firsta->NZE_C_N)
            {
              firsta = firsta->NZE_C_N;
              i++;
            }
          if (firsta->lag_index == lag) i++;
        }
    }
  return i;
}

void
SparseMatrix::Delete(const int r, const int c)
{
  NonZeroElem *first = FNZE_R[r], *firsta = NULL;

  while (first->c_index != c)
    {
      firsta = first;
      first = first->NZE_R_N;
    }
  if (firsta != NULL)
    firsta->NZE_R_N = first->NZE_R_N;
  if (first == FNZE_R[r])
    FNZE_R[r] = first->NZE_R_N;
  NbNZRow[r]--;

  first = FNZE_C[c];
  firsta = NULL;
  while (first->r_index != r)
    {
      firsta = first;
      first = first->NZE_C_N;
    }

  if (firsta != NULL)
    firsta->NZE_C_N = first->NZE_C_N;
  if (first == FNZE_C[c])
    FNZE_C[c] = first->NZE_C_N;

  u_liste.push_back(first->u_index);
  mem_mngr.mxFree_NZE(first);
  NbNZCol[c]--;
}

void
SparseMatrix::Print(int Size, int *b)
{
  int a, i, j, k, l;
  mexPrintf("   ");
  for (k = 0; k < Size*periods; k++)
    mexPrintf("%-2d ", k);
  mexPrintf("    |    ");
  for (k = 0; k < Size*periods; k++)
    mexPrintf("%8d", k);
  mexPrintf("\n");
  for (i = 0; i < Size*periods; i++)
    {
      NonZeroElem *first = FNZE_R[i];
      j = NbNZRow[i];
      mexPrintf("%-2d ", i);
      a = 0;
      for (k = 0; k < j; k++)
        {
          for (l = 0; l < (first->c_index-a); l++)
            mexPrintf("   ");
          mexPrintf("%-2d ", first->u_index);
          a = first->c_index+1;
          first = first->NZE_R_N;
        }
      for (k = a; k < Size*periods; k++)
        mexPrintf("   ");
      mexPrintf("%-2d ", b[i]);

      first = FNZE_R[i];
      j = NbNZRow[i];
      mexPrintf(" | %-2d ", i);
      a = 0;
      for (k = 0; k < j; k++)
        {
          for (l = 0; l < (first->c_index-a); l++)
            mexPrintf("        ");
          mexPrintf("%8.4f", double (u[first->u_index]));
          a = first->c_index+1;
          first = first->NZE_R_N;
        }
      for (k = a; k < Size*periods; k++)
        mexPrintf("        ");
      mexPrintf("%8.4f", double (u[b[i]]));
      mexPrintf("\n");
    }
}

void
SparseMatrix::Insert(const int r, const int c, const int u_index, const int lag_index)
{
  NonZeroElem *firstn, *first, *firsta, *a;
  firstn = mem_mngr.mxMalloc_NZE();
  first = FNZE_R[r];
  firsta = NULL;
  while (first->c_index < c && (a = first->NZE_R_N))
    {
      firsta = first;
      first = a;
    }
  firstn->u_index = u_index;
  firstn->r_index = r;
  firstn->c_index = c;
  firstn->lag_index = lag_index;
  if (first->c_index > c)
    {
      if (first == FNZE_R[r])
        FNZE_R[r] = firstn;
      if (firsta != NULL)
        firsta->NZE_R_N = firstn;
      firstn->NZE_R_N = first;
    }
  else
    {
      first->NZE_R_N = firstn;
      firstn->NZE_R_N = NULL;
    }
  NbNZRow[r]++;
  first = FNZE_C[c];
  firsta = NULL;
  while (first->r_index < r && (a = first->NZE_C_N))
    {
      firsta = first;
      first = a;
    }
  if (first->r_index > r)
    {
      if (first == FNZE_C[c])
        FNZE_C[c] = firstn;
      if (firsta != NULL)
        firsta->NZE_C_N = firstn;
      firstn->NZE_C_N = first;
    }
  else
    {
      first->NZE_C_N = firstn;
      firstn->NZE_C_N = NULL;
    }

  NbNZCol[c]++;
}

void
SparseMatrix::Read_SparseMatrix(string file_name, const int Size, int periods, int y_kmin, int y_kmax, bool steady_state, bool two_boundaries)
{
  unsigned int eq, var;
  int i, j, lag;
  filename = file_name;
  mem_mngr.fixe_file_name(file_name);
  if (!SaveCode.is_open())
    {
      if (steady_state)
        SaveCode.open((file_name + "_static.bin").c_str(), ios::in | ios::binary);
      else
        SaveCode.open((file_name + "_dynamic.bin").c_str(), ios::in | ios::binary);
      if (!SaveCode.is_open())
        {
          if (steady_state)
            mexPrintf("Error : Can't open file \"%s\" for reading\n", (file_name + "_static.bin").c_str());
          else
            mexPrintf("Error : Can't open file \"%s\" for reading\n", (file_name + "_dynamic.bin").c_str());
          mexEvalString("st=fclose('all');clear all;");
          mexErrMsgTxt("Exit from Dynare");
        }
    }
  IM_i.clear();
  if (two_boundaries)
    {
      for (i = 0; i < u_count_init-Size; i++)
        {
          SaveCode.read(reinterpret_cast<char *>(&eq), sizeof(eq));
          SaveCode.read(reinterpret_cast<char *>(&var), sizeof(var));
          SaveCode.read(reinterpret_cast<char *>(&lag), sizeof(lag));
          SaveCode.read(reinterpret_cast<char *>(&j), sizeof(j));
          IM_i[make_pair(make_pair(eq, var), lag)] = j;
        }
      for (j = 0; j < Size; j++)
        IM_i[make_pair(make_pair(j, Size*(periods+y_kmax)), 0)] = j;
    }
  else
    {
      for (i = 0; i < u_count_init; i++)
        {
          SaveCode.read(reinterpret_cast<char *>(&eq), sizeof(eq));
          SaveCode.read(reinterpret_cast<char *>(&var), sizeof(var));
          SaveCode.read(reinterpret_cast<char *>(&lag), sizeof(lag));
          SaveCode.read(reinterpret_cast<char *>(&j), sizeof(j));
          IM_i[make_pair(make_pair(eq, var), lag)] = j;
        }
    }
  index_vara = (int *) mxMalloc(Size*(periods+y_kmin+y_kmax)*sizeof(int));
  for (j = 0; j < Size; j++)
    SaveCode.read(reinterpret_cast<char *>(&index_vara[j]), sizeof(*index_vara));
  if (periods+y_kmin+y_kmax > 1)
    for (i = 1; i < periods+y_kmin+y_kmax; i++)
      for (j = 0; j < Size; j++)
        index_vara[j+Size*i] = index_vara[j+Size*(i-1)]+y_size;
  index_equa = (int *) mxMalloc(Size*sizeof(int));
  for (j = 0; j < Size; j++)
    SaveCode.read(reinterpret_cast<char *>(&index_equa[j]), sizeof(*index_equa));
}

void
SparseMatrix::Simple_Init(int it_, int y_kmin, int y_kmax, int Size, map<pair<pair<int, int>, int>, int> &IM)
{
  int i, eq, var, lag;
  map<pair<pair<int, int>, int>, int>::iterator it4;
  NonZeroElem *first;
  pivot = (int *) mxMalloc(Size*sizeof(int));
  pivot_save = (int *) mxMalloc(Size*sizeof(int));
  pivotk = (int *) mxMalloc(Size*sizeof(int));
  pivotv = (double *) mxMalloc(Size*sizeof(double));
  pivotva = (double *) mxMalloc(Size*sizeof(double));
  b = (int *) mxMalloc(Size*sizeof(int));
  line_done = (bool *) mxMalloc(Size*sizeof(bool));

  mem_mngr.init_CHUNK_BLCK_SIZE(u_count);
  g_save_op = NULL;
  g_nop_all = 0;
  i = Size*sizeof(NonZeroElem *);
  FNZE_R = (NonZeroElem **) mxMalloc(i);
  FNZE_C = (NonZeroElem **) mxMalloc(i);
  NonZeroElem **temp_NZE_R = (NonZeroElem **) mxMalloc(i);
  NonZeroElem **temp_NZE_C = (NonZeroElem **) mxMalloc(i);
  i = Size*sizeof(int);
  NbNZRow = (int *) mxMalloc(i);
  NbNZCol = (int *) mxMalloc(i);
  it4 = IM.begin();
  eq = -1;
  //#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
  for (i = 0; i < Size; i++)
    {
      line_done[i] = 0;
      FNZE_C[i] = 0;
      FNZE_R[i] = 0;
      temp_NZE_C[i] = 0;
      temp_NZE_R[i] = 0;
      NbNZRow[i] = 0;
      NbNZCol[i] = 0;
    }
  int u_count1 = Size;
  while (it4 != IM.end())
    {
      var = it4->first.first.second;
      eq = it4->first.first.first;
      lag = it4->first.second;
      if (lag == 0)   /*Build the index for sparse matrix containing the jacobian : u*/
        {
          NbNZRow[eq]++;
          NbNZCol[var]++;
          first = mem_mngr.mxMalloc_NZE();
          first->NZE_C_N = NULL;
          first->NZE_R_N = NULL;
          first->u_index = u_count1;
          first->r_index = eq;
          first->c_index = var;
          first->lag_index = lag;
          if (FNZE_R[eq] == NULL)
            FNZE_R[eq] = first;
          if (FNZE_C[var] == NULL)
            FNZE_C[var] = first;
          if (temp_NZE_R[eq] != NULL)
            temp_NZE_R[eq]->NZE_R_N = first;
          if (temp_NZE_C[var] != NULL)
            temp_NZE_C[var]->NZE_C_N = first;
          temp_NZE_R[eq] = first;
          temp_NZE_C[var] = first;
          u_count1++;
        }
      it4++;
    }
  //#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
  for (i = 0; i < Size; i++)
    b[i] = i;
  mxFree(temp_NZE_R);
  mxFree(temp_NZE_C);
  u_count = u_count1;
}

void
SparseMatrix::Init(int periods, int y_kmin, int y_kmax, int Size, map<pair<pair<int, int>, int>, int> &IM)
{
  int t, i, eq, var, lag, ti_y_kmin, ti_y_kmax;
  double tmp_b = 0.0;
  map<pair<pair<int, int>, int>, int>::iterator it4;
  NonZeroElem *first;
  pivot = (int *) mxMalloc(Size*periods*sizeof(int));
  pivot_save = (int *) mxMalloc(Size*periods*sizeof(int));
  pivotk = (int *) mxMalloc(Size*periods*sizeof(int));
  pivotv = (double *) mxMalloc(Size*periods*sizeof(double));
  pivotva = (double *) mxMalloc(Size*periods*sizeof(double));
  b = (int *) mxMalloc(Size*periods*sizeof(int));
  line_done = (bool *) mxMalloc(Size*periods*sizeof(bool));
  mem_mngr.init_CHUNK_BLCK_SIZE(u_count);
  g_save_op = NULL;
  g_nop_all = 0;
  i = (periods+y_kmax+1)*Size*sizeof(NonZeroElem *);
  FNZE_R = (NonZeroElem **) mxMalloc(i);
  FNZE_C = (NonZeroElem **) mxMalloc(i);
  NonZeroElem **temp_NZE_R = (NonZeroElem **) mxMalloc(i);
  NonZeroElem **temp_NZE_C = (NonZeroElem **) mxMalloc(i);
  i = (periods+y_kmax+1)*Size*sizeof(int);
  NbNZRow = (int *) mxMalloc(i);
  NbNZCol = (int *) mxMalloc(i);

  //#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
  for (i = 0; i < periods*Size; i++)
    {
      b[i] = 0;
      line_done[i] = 0;
    }
  //#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
  for (i = 0; i < (periods+y_kmax+1)*Size; i++)
    {
      FNZE_C[i] = 0;
      FNZE_R[i] = 0;
      temp_NZE_C[i] = NULL;
      temp_NZE_R[i] = NULL;
      NbNZRow[i] = 0;
      NbNZCol[i] = 0;
    }

  //#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) ordered private(it4, ti_y_kmin, ti_y_kmax, eq, var, lag) schedule(dynamic)
  for (t = 0; t < periods; t++)
    {
      ti_y_kmin = -min(t, y_kmin);
      ti_y_kmax = min(periods-(t+1), y_kmax);
      it4 = IM.begin();
      eq = -1;
      //#pragma omp ordered
      while (it4 != IM.end())
        {
          var = it4->first.first.second;
          if (eq != it4->first.first.first+Size*t)
            tmp_b = 0;
          eq = it4->first.first.first+Size*t;
          lag = it4->first.second;
          if (var < (periods+y_kmax)*Size)
            {
              lag = it4->first.second;
              if (lag <= ti_y_kmax && lag >= ti_y_kmin)   /*Build the index for sparse matrix containing the jacobian : u*/
                {
                  var += Size*t;
                  NbNZRow[eq]++;
                  NbNZCol[var]++;
                  first = mem_mngr.mxMalloc_NZE();
                  first->NZE_C_N = NULL;
                  first->NZE_R_N = NULL;
                  first->u_index = it4->second+u_count_init*t;
                  first->r_index = eq;
                  first->c_index = var;
                  first->lag_index = lag;
                  if (FNZE_R[eq] == NULL)
                    FNZE_R[eq] = first;
                  if (FNZE_C[var] == NULL)
                    FNZE_C[var] = first;
                  if (temp_NZE_R[eq] != NULL)
                    temp_NZE_R[eq]->NZE_R_N = first;
                  if (temp_NZE_C[var] != NULL)
                    temp_NZE_C[var]->NZE_C_N = first;
                  temp_NZE_R[eq] = first;
                  temp_NZE_C[var] = first;
                }
              else       /*Build the additive terms ooutside the simulation periods related to the first lags and the last leads...*/
                {
                  if (lag < ti_y_kmin)
                    {
                      tmp_b += u[it4->second+u_count_init*t]*y[index_vara[var+Size*(y_kmin+t)]];
                    }
                  else
                    {
                      tmp_b += u[it4->second+u_count_init*t]*y[index_vara[var+Size*(y_kmin+t)]];

                    }
                }
            }
          else           /* ...and store it in the u vector*/
            {
              b[eq] = it4->second+u_count_init*t;
              u[b[eq]] += tmp_b;
              tmp_b = 0;
            }
          it4++;
        }
    }
  mxFree(temp_NZE_R);
  mxFree(temp_NZE_C);
}

void
SparseMatrix::ShortInit(int periods, int y_kmin, int y_kmax, int Size, map<pair<pair<int, int>, int>, int> &IM)
{
  int t, eq, var, lag, ti_y_kmin, ti_y_kmax;
  double tmp_b = 0.0;
  map<pair<pair<int, int>, int>, int>::iterator it4;
  //#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) ordered private(it4, ti_y_kmin, ti_y_kmax, eq, var, lag, tmp_b) schedule(dynamic)
  for (t = 0; t < periods; t++)
    {
      ti_y_kmin = -min(t, y_kmin);
      ti_y_kmax = min(periods-(t+1), y_kmax);
      it4 = IM.begin();
      eq = -1;
      while (it4 != IM.end())
        {
          var = it4->first.first.second;
          if (eq != it4->first.first.first+Size*t)
            tmp_b = 0;
          eq = it4->first.first.first+Size*t;
          if (var < (periods+y_kmax)*Size)
            {
              lag = it4->first.second;
              if (lag <= ti_y_kmax && lag >= ti_y_kmin)
                {
                  var += Size*t;
                }
              else
                {
                  tmp_b += u[it4->second+u_count_init*t]*y[index_vara[var+Size*(y_kmin+t)]];
                }
            }
          else
            {
              b[eq] = it4->second+u_count_init*t;
              u[b[eq]] += tmp_b;
            }
          it4++;
        }
    }
}

int
SparseMatrix::Get_u()
{
  if (!u_liste.empty())
    {
      int i = u_liste.back();
      u_liste.pop_back();
      return i;
    }
  else
    {
      if (u_count < u_count_alloc)
        {
          int i = u_count;
          u_count++;
          return i;
        }
      else
        {
          u_count_alloc += 5*u_count_alloc_save;
          u = (double *) mxRealloc(u, u_count_alloc*sizeof(double));
          if (!u)
            {
              mexPrintf("Error in Get_u: memory exhausted (realloc(%d))\n", u_count_alloc*sizeof(double));
              mexEvalString("st=fclose('all');clear all;");
              mexErrMsgTxt("Exit from Dynare");
            }
          int i = u_count;
          u_count++;
          return i;
        }
    }
}

void
SparseMatrix::Delete_u(int pos)
{
  u_liste.push_back(pos);
}

void
SparseMatrix::Clear_u()
{
  u_liste.clear();
}

void
SparseMatrix::Print_u()
{
  for (unsigned int i = 0; i < u_liste.size(); i++)
    mexPrintf("%d ", u_liste[i]);
}

void
SparseMatrix::End(int Size)
{
  mem_mngr.Free_All();
  mxFree(FNZE_R);
  mxFree(FNZE_C);
  mxFree(NbNZRow);
  mxFree(NbNZCol);
  mxFree(b);
  mxFree(line_done);
  mxFree(pivot);
  mxFree(pivot_save);
  mxFree(pivotk);
  mxFree(pivotv);
  mxFree(pivotva);
}

bool
SparseMatrix::compare(int *save_op, int *save_opa, int *save_opaa, int beg_t, int periods, long int nop4,  int Size)
{
  long int i, j, nop = nop4/2, t, k;
  double r = 0.0;
  bool OK = true;
  clock_t t001;
  t_save_op_s *save_op_s, *save_opa_s, *save_opaa_s;
  int *diff1, *diff2;
  t001 = clock();
  diff1 = (int *) mxMalloc(nop*sizeof(int));
  diff2 = (int *) mxMalloc(nop*sizeof(int));
  int max_save_ops_first = -1;
  j = k = i = 0;
  while (i < nop4 && OK)
    {
      save_op_s = (t_save_op_s *) &(save_op[i]);
      save_opa_s = (t_save_op_s *) &(save_opa[i]);
      save_opaa_s = (t_save_op_s *) &(save_opaa[i]);
      diff1[j] = save_op_s->first-save_opa_s->first;
      if (max_save_ops_first < save_op_s->first+diff1[j]*(periods-beg_t))
        {
          max_save_ops_first = save_op_s->first+diff1[j]*(periods-beg_t);
        }
      switch (save_op_s->operat)
        {
        case IFLD:
        case IFDIV:
          OK = (save_op_s->operat == save_opa_s->operat && save_opa_s->operat == save_opaa_s->operat
                && diff1[j] == (save_opa_s->first-save_opaa_s->first));
          i += 2;
          break;
        case IFLESS:
        case IFSUB:
          diff2[j] = save_op_s->second-save_opa_s->second;
          OK = (save_op_s->operat == save_opa_s->operat && save_opa_s->operat == save_opaa_s->operat
                && diff1[j] == (save_opa_s->first-save_opaa_s->first)
                && diff2[j] == (save_opa_s->second-save_opaa_s->second));
          i += 3;
          break;
        default:
          mexPrintf("unknown operator = %d ", save_op_s->operat);
          mexEvalString("st=fclose('all');clear all;");
          filename += " stopped";
          mexErrMsgTxt(filename.c_str());
          break;
        }
      j++;
    }
  // the same pivot for all remaining periods
  if (OK)
    {
      //#pragma omp parallel for  num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) ordered private(j) schedule(dynamic)
      for (i = beg_t; i < periods; i++)
        {
          for (j = 0; j < Size; j++)
            {
              ///#pragma omp ordered
              pivot[i*Size+j] = pivot[(i-1)*Size+j]+Size;
            }
        }
      if (max_save_ops_first >= u_count_alloc)
        {
          u_count_alloc += max_save_ops_first;
          u = (double *) mxRealloc(u, u_count_alloc*sizeof(double));
          if (!u)
            {
              mexPrintf("Error in Get_u: memory exhausted (realloc(%d))\n", u_count_alloc*sizeof(double));
              mexEvalString("st=fclose('all');clear all;");
              mexErrMsgTxt("Exit from Dynare");
            }
        }
      double *up;
      for (t = 1; t < periods-beg_t-y_kmax; t++)
        {
          i = j = 0;
          while (i < nop4)
            {
              save_op_s = (t_save_op_s *)(&(save_op[i]));
              up = &u[save_op_s->first+t*diff1[j]];
              switch (save_op_s->operat)
                {
                case IFLD:
                  r = *up;
                  i += 2;
                  break;
                case IFDIV:
                  *up /= r;
                  i += 2;
                  break;
                case IFSUB:
                  *up -= u[save_op_s->second+t*diff2[j]]*r;;
                  i += 3;
                  break;
                case IFLESS:
                  *up = -u[save_op_s->second+t*diff2[j]]*r;
                  i += 3;
                  break;
                }
              j++;
            }
        }
      int t1 = max(1, periods-beg_t-y_kmax);
      int periods_beg_t = periods-beg_t;
      for (t = t1; t < periods_beg_t; t++)
        {
          i = j = 0;
          while (i < nop4)
            {
              save_op_s = (t_save_op_s *)(&(save_op[i]));
              if (save_op_s->lag < (periods_beg_t-t))
                {
                  up = &u[save_op_s->first+t*diff1[j]];
                  switch (save_op_s->operat)
                    {
                    case IFLD:
                      r = *up;
                      i += 2;
                      break;
                    case IFDIV:
                      *up /= r;
                      i += 2;
                      break;
                    case IFSUB:
                      *up -= u[save_op_s->second+t*diff2[j]]*r;
                      i += 3;
                      break;
                    case IFLESS:
                      *up = -u[save_op_s->second+t*diff2[j]]*r;
                      i += 3;
                      break;
                    }
                }
              else
                {
                  switch (save_op_s->operat)
                    {
                    case IFLD:
                    case IFDIV:
                      i += 2;
                      break;
                    case IFSUB:
                    case IFLESS:
                      i += 3;
                      break;
                    }
                }
              j++;
            }
        }
    }
  mxFree(diff1);
  mxFree(diff2);
  return OK;
}

int
SparseMatrix::complete(int beg_t, int Size, int periods, int *b)
{
  long int i, j, k, nop, nopa, nop1, cal_y, nb_var, pos, t, ti, max_var, min_var;
  NonZeroElem *first;
  int *save_code;
  int *diff;
  double yy = 0.0, err;

  int size_of_save_code = (1+y_kmax)*Size*(Size+1+4)/2*4;
  save_code = (int *) mxMalloc(size_of_save_code*sizeof(int));
  int size_of_diff = (1+y_kmax)*Size*(Size+1+4);
  diff = (int *) mxMalloc(size_of_diff*sizeof(int));
  cal_y = y_size*y_kmin;

  i = (beg_t+1)*Size-1;
  nop = 0;
  for (j = i; j > i-Size; j--)
    {
      pos = pivot[j];
      nb_var = At_Row(pos, &first);
      first = first->NZE_R_N;
      nb_var--;
      save_code[nop] = IFLDZ;
      save_code[nop+1] = 0;
      save_code[nop+2] = 0;
      save_code[nop+3] = 0;
      if ((nop+3) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop+2, size_of_save_code);
      nop += 4;
      for (k = 0; k < nb_var; k++)
        {
          save_code[nop] = IFMUL;
          save_code[nop+1] = index_vara[first->c_index]+cal_y;
          save_code[nop+2] = first->u_index;
          save_code[nop+3] = first->lag_index;
          if ((nop+3) >= size_of_save_code)
            mexPrintf("out of save_code[%d] (bound=%d)\n", nop+2, size_of_save_code);
          nop += 4;
          first = first->NZE_R_N;
        }
      save_code[nop] = IFADD;
      save_code[nop+1] = b[pos];
      save_code[nop+2] = 0;
      save_code[nop+3] = 0;
      if ((nop+3) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop+2, size_of_save_code);
      nop += 4;
      save_code[nop] = IFSTP;
      save_code[nop+1] = index_vara[j]+y_size*y_kmin;
      save_code[nop+2] = 0;
      save_code[nop+3] = 0;
      if ((nop+2) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop+2, size_of_save_code);
      nop += 4;
    }
  i = beg_t*Size-1;
  nop1 = nopa = 0;
  for (j = i; j > i-Size; j--)
    {
      pos = pivot[j];
      nb_var = At_Row(pos, &first);
      first = first->NZE_R_N;
      nb_var--;
      diff[nopa] = 0;
      diff[nopa+1] = 0;
      nopa += 2;
      nop1 += 4;
      for (k = 0; k < nb_var; k++)
        {
          diff[nopa] = save_code[nop1+1]-(index_vara[first->c_index]+cal_y);
          diff[nopa+1] = save_code[nop1+2]-(first->u_index);
          if ((nop1+2) >= size_of_save_code)
            mexPrintf("out of save_code[%d] (bound=%d)\n", nop1+2, size_of_save_code);
          if ((nopa+1) >= size_of_diff)
            mexPrintf("out of diff[%d] (bound=%d)\n", nopa+2, size_of_diff);
          nopa += 2;
          nop1 += 4;
          first = first->NZE_R_N;
        }
      diff[nopa] = save_code[nop1+1]-(b[pos]);
      diff[nopa+1] = 0;
      if ((nop1+3) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop1+2, size_of_save_code);
      if ((nopa+1) >= size_of_diff)
        mexPrintf("out of diff[%d] (bound=%d)\n", nopa+2, size_of_diff);
      nopa += 2;
      nop1 += 4;
      diff[nopa] = save_code[nop1+1]-(index_vara[j]+y_size*y_kmin);
      diff[nopa+1] = 0;
      if ((nop1+4) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop1+2, size_of_save_code);
      if ((nopa+1) >= size_of_diff)
        mexPrintf("out of diff[%d] (bound=%d)\n", nopa+2, size_of_diff);
      nopa += 2;
      nop1 += 4;
    }
  max_var = (periods+y_kmin)*y_size;
  min_var = y_kmin*y_size;
  int k1 = 0;
  for (t = periods+y_kmin-1; t >= beg_t+y_kmin; t--)
    {
      j = 0;
      ti = t-y_kmin-beg_t;
      for (i = 0; i < nop; i += 4)
        {
          switch (save_code[i])
            {
            case IFLDZ:
              yy = 0;
              break;
            case IFMUL:
              k = save_code[i+1]+ti*diff[j];
              if (k < max_var && k > min_var)
                {
                  yy += y[k]*u[save_code[i+2]+ti*diff[j+1]];
                }
              break;
            case IFADD:
              yy = -(yy+u[save_code[i+1]+ti*diff[j]]);
              break;
            case IFSTP:
              k = save_code[i+1]+ti*diff[j];
              k1 = k;
              err = yy - y[k];
              y[k] += slowc*(err);
              break;
            }
          j += 2;
        }
    }
  mxFree(save_code);
  mxFree(diff);
  return (beg_t);
}

double
SparseMatrix::bksub(int tbreak, int last_period, int Size, double slowc_l)
{
  NonZeroElem *first;
  int i, j, k;
  double yy;
  for (i = 0; i < y_size*(periods+y_kmin); i++)
    y[i] = ya[i];
  if (symbolic && tbreak)
    last_period = complete(tbreak, Size, periods, b);
  else
    last_period = periods;
  for (int t = last_period+y_kmin-1; t >= y_kmin; t--)
    {
      int ti = (t-y_kmin)*Size;
      int cal = y_kmin*Size;
      int cal_y = y_size*y_kmin;
      for (i = ti-1; i >= ti-Size; i--)
        {
          j = i+cal;
          int pos = pivot[i+Size];
          int nb_var = At_Row(pos, &first);
          first = first->NZE_R_N;
          nb_var--;
          int eq = index_vara[j]+y_size;
          yy = 0;
          for (k = 0; k < nb_var; k++)
            {
              yy += y[index_vara[first->c_index]+cal_y]*u[first->u_index];
              first = first->NZE_R_N;
            }
          yy = -(yy+y[eq]+u[b[pos]]);
          direction[eq] = yy;
          y[eq] += slowc_l*yy;
        }
    }
  return res1;
}

double
SparseMatrix::simple_bksub(int it_, int Size, double slowc_l)
{
  int i, k;
  double yy;
  NonZeroElem *first;
  for (i = 0; i < y_size; i++)
    y[i+it_*y_size] = ya[i+it_*y_size];
  for (i = Size-1; i >= 0; i--)
    {
      int pos = pivot[i];
      int nb_var = At_Row(pos, &first);
      first = first->NZE_R_N;
      nb_var--;
      int eq = index_vara[i];
      yy = 0;
      for (k = 0; k < nb_var; k++)
        {
          yy += y[index_vara[first->c_index]+it_*y_size]*u[first->u_index];
          first = first->NZE_R_N;
        }
      yy = -(yy+y[eq+it_*y_size]+u[b[pos]]);
      direction[eq+it_*y_size] = yy;
      y[eq+it_*y_size] += slowc_l*yy;
    }
  return res1;
}

bool
SparseMatrix::simulate_NG(int blck, int y_size, int it_, int y_kmin, int y_kmax, int Size, bool print_it, bool cvg, int &iter, bool steady_state, int Block_number)
{
  int i, j, k;
  int pivj = 0, pivk = 0;
  double piv_abs;
  NonZeroElem *first, *firsta, *first_suba;
  double *piv_v;
  int *pivj_v, *pivk_v, *NR;
  int l, N_max;
  bool one;
  Clear_u();
  piv_v = (double *) mxMalloc(Size*sizeof(double));
  pivj_v = (int *) mxMalloc(Size*sizeof(int));
  pivk_v = (int *) mxMalloc(Size*sizeof(int));
  NR = (int *) mxMalloc(Size*sizeof(int));
  error_not_printed = true;
  u_count_alloc_save = u_count_alloc;
  if (isnan(res1) || isinf(res1) || (res2 > 12*g0 && iter>0))
    {
      if (iter == 0)
        {
          for (j = 0; j < y_size; j++)
            {
              bool select = false;
              for (int i = 0; i < Size; i++)
                if (j == index_vara[i])
                  {
                    select = true;
                    break;
                  }
              if (select)
                mexPrintf("-> variable %d at time %d = %f direction = %f\n", j+1, it_, y[j+it_*y_size], direction[j+it_*y_size]);
              else
                mexPrintf("   variable %d at time %d = %f direction = %f\n", j+1, it_, y[j+it_*y_size], direction[j+it_*y_size]);
            }
          mexPrintf("res1=%5.10f\n", res1);
          mexPrintf("The initial values of endogenous variables are too far from the solution.\n");
          mexPrintf("Change them!\n");
          mexEvalString("drawnow;");
          mxFree(piv_v);
          mxFree(pivj_v);
          mxFree(pivk_v);
          mxFree(NR);
          if (steady_state)
            return false;
          else
            {
              mexEvalString("st=fclose('all');clear all;");
              filename += " stopped";
              mexErrMsgTxt(filename.c_str());
            }
        }
      if (fabs(slowc_save) < 1e-8)
        {
          for (j = 0; j < y_size; j++)
            {
              bool select = false;
              for (int i = 0; i < Size; i++)
                if (j == index_vara[i])
                  {
                    select = true;
                    break;
                  }
              if (select)
                mexPrintf("-> variable %d at time %d = %f direction = %f\n", j+1, it_, y[j+it_*y_size], direction[j+it_*y_size]);
              else
                mexPrintf("   variable %d at time %d = %f direction = %f\n", j+1, it_, y[j+it_*y_size], direction[j+it_*y_size]);
            }
          mexPrintf("Dynare cannot improve the simulation in block %d at time %d (variable %d)\n", blck+1, it_+1, max_res_idx);
          mexEvalString("drawnow;");
          mxFree(piv_v);
          mxFree(pivj_v);
          mxFree(pivk_v);
          mxFree(NR);
          if (steady_state)
            return false;
          else
            {
              mexEvalString("st=fclose('all');clear all;");
              filename += " stopped";
              mexErrMsgTxt(filename.c_str());
            }
        }
      if(!(isnan(res1) || isinf(res1)) && !(isnan(g0) || isinf(g0)))
        {
          if (try_at_iteration == 0)
            {
              prev_slowc_save = slowc_save;
              slowc_save = max( - gp0 / (2 * (res2 - g0 - gp0)) , 0.1);
            }
          else
            {
              double t1 = res2 - gp0 * slowc_save - g0;
              double t2 = glambda2 - gp0 * prev_slowc_save - g0;
              double a = (1/(slowc_save * slowc_save) * t1 - 1/(prev_slowc_save * prev_slowc_save) * t2) / (slowc_save - prev_slowc_save);
              double b = (-prev_slowc_save/(slowc_save * slowc_save) * t1 + slowc_save/(prev_slowc_save * prev_slowc_save) * t2) / (slowc_save - prev_slowc_save);
              prev_slowc_save = slowc_save;
              slowc_save = max(min( -b + sqrt(b*b - 3 * a * gp0) / (3 * a), 0.5 * slowc_save), 0.1 * slowc_save);
            }
          glambda2 = res2;
          try_at_iteration ++;
        }
      else
        {
          prev_slowc_save = slowc_save;
          slowc_save /= 1.1;
        }
      if (print_it)
        mexPrintf("Error: Simulation diverging, trying to correct it using slowc=%f\n", slowc_save);
      for (i = 0; i < y_size; i++)
        y[i+it_*y_size] = ya[i+it_*y_size] + slowc_save*direction[i+it_*y_size];
      mxFree(piv_v);
      mxFree(pivj_v);
      mxFree(pivk_v);
      mxFree(NR);
      iter--;
      return true;
    }
  if (cvg)
    {
      mxFree(piv_v);
      mxFree(pivj_v);
      mxFree(pivk_v);
      mxFree(NR);
      return (true);
    }
  if (print_it )
    {
      mexPrintf("solwc=%f g0=%f res2=%f glambda2=%f\n",slowc_save,g0, res2, glambda2);
      mexPrintf("-----------------------------------\n");
      mexPrintf("      Simulate iteration no %d     \n", iter+1);
      mexPrintf("      max. error=%.10e       \n", double (max_res));
      mexPrintf("      sqr. error=%.10e       \n", double (res2));
      mexPrintf("      abs. error=%.10e       \n", double (res1));
      mexPrintf("-----------------------------------\n");
    }
  Simple_Init(it_, y_kmin, y_kmax, Size, IM_i);
  NonZeroElem **bc;
  bc = (NonZeroElem **) mxMalloc(Size*sizeof(*bc));
  for (i = 0; i < Size; i++)
    {
      /*finding the max-pivot*/
      double piv = piv_abs = 0;
      int nb_eq = At_Col(i, &first);
      l = 0; N_max = 0;
      one = false;
      piv_abs = 0;
      for (j = 0; j < nb_eq; j++)
        {
          if (!line_done[first->r_index])
            {
              k = first->u_index;
              int jj = first->r_index;
              int NRow_jj = NRow(jj);

              piv_v[l] = u[k];
              double piv_fabs = fabs(u[k]);
              pivj_v[l] = jj;
              pivk_v[l] = k;
              NR[l] = NRow_jj;
              if (NRow_jj == 1 && !one)
                {
                  one = true;
                  piv_abs = piv_fabs;
                  N_max = NRow_jj;
                }
              if (!one)
                {
                  if (piv_fabs > piv_abs)
                    piv_abs = piv_fabs;
                  if (NRow_jj > N_max)
                    N_max = NRow_jj;
                }
              else
                {
                  if (NRow_jj == 1)
                    {
                      if (piv_fabs > piv_abs)
                        piv_abs = piv_fabs;
                      if (NRow_jj > N_max)
                        N_max = NRow_jj;
                    }
                }
              l++;
            }
          first = first->NZE_C_N;
        }
      if (piv_abs < eps)
        {
          if (Block_number > 1)
            mexPrintf("Error: singular system in Simulate_NG in block %d\n", blck+1);
          else
            mexPrintf("Error: singular system in Simulate_NG\n");
          mexEvalString("drawnow;");
          mxFree(piv_v);
          mxFree(pivj_v);
          mxFree(pivk_v);
          mxFree(NR);
          mxFree(bc);
          if (steady_state)
            return false;
          else
            {
              mexEvalString("st=fclose('all');clear all;");
              filename += " stopped";
              mexErrMsgTxt(filename.c_str());
            }
        }
      double markovitz = 0, markovitz_max = -9e70;
      int NR_max = 0;
      if (!one)
        {
          for (j = 0; j < l; j++)
            {
              if (N_max > 0 && NR[j] > 0)
                {
                  if (fabs(piv_v[j]) > 0)
                    {
                      if (markowitz_c > 0)
                        markovitz = exp(log(fabs(piv_v[j])/piv_abs)-markowitz_c*log(double (NR[j])/double (N_max)));
                      else
                        markovitz = fabs(piv_v[j])/piv_abs;
                    }
                  else
                    markovitz = 0;
                }
              else
                markovitz = fabs(piv_v[j])/piv_abs;
              if (markovitz > markovitz_max)
                {
                  piv = piv_v[j];
                  pivj = pivj_v[j];   //Line number
                  pivk = pivk_v[j];   //positi
                  markovitz_max = markovitz;
                  NR_max = NR[j];
                }
            }
        }
      else
        {
          for (j = 0; j < l; j++)
            {
              if (N_max > 0 && NR[j] > 0)
                {
                  if (fabs(piv_v[j]) > 0)
                    {
                      if (markowitz_c > 0)
                        markovitz = exp(log(fabs(piv_v[j])/piv_abs)-markowitz_c*log(double (NR[j])/double (N_max)));
                      else
                        markovitz = fabs(piv_v[j])/piv_abs;
                    }
                  else
                    markovitz = 0;
                }
              else
                markovitz = fabs(piv_v[j])/piv_abs;
              if (NR[j] == 1)
                {
                  piv = piv_v[j];
                  pivj = pivj_v[j];   //Line number
                  pivk = pivk_v[j];   //positi
                  markovitz_max = markovitz;
                  NR_max = NR[j];
                }
            }
        }
      /*if (fabs(piv) < eps)
        mexPrintf("==> Error NR_max=%d, N_max=%d and piv=%f, piv_abs=%f, markovitz_max=%f\n",NR_max, N_max, piv, piv_abs, markovitz_max);
      if (NR_max == 0)
        mexPrintf("==> Error NR_max=0 and piv=%f, markovitz_max=%f\n",piv, markovitz_max);*/
      pivot[i] = pivj;
      pivotk[i] = pivk;
      pivotv[i] = piv;
      line_done[pivj] = true;

      /*divide all the non zeros elements of the line pivj by the max_pivot*/
      int nb_var = At_Row(pivj, &first);
      for (j = 0; j < nb_var; j++)
        {
          u[first->u_index] /= piv;
          first = first->NZE_R_N;
        }
      u[b[pivj]] /= piv;
      /*substract the elements on the non treated lines*/
      nb_eq = At_Col(i, &first);
      NonZeroElem *first_piva;
      int nb_var_piva = At_Row(pivj, &first_piva);
      int nb_eq_todo = 0;
      for (j = 0; j < nb_eq && first; j++)
        {
          if (!line_done[first->r_index])
            bc[nb_eq_todo++] = first;
          first = first->NZE_C_N;
        }
      //#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
      for (j = 0; j < nb_eq_todo; j++)
        {
          first = bc[j];
          int row = first->r_index;
          double first_elem = u[first->u_index];

          int nb_var_piv = nb_var_piva;
          NonZeroElem *first_piv = first_piva;
          NonZeroElem *first_sub;
          int nb_var_sub = At_Row(row, &first_sub);
          int l_sub = 0, l_piv = 0;
          int sub_c_index = first_sub->c_index, piv_c_index = first_piv->c_index;
          while (l_sub < nb_var_sub || l_piv < nb_var_piv)
            {
              if (l_sub < nb_var_sub && (sub_c_index < piv_c_index || l_piv >= nb_var_piv))
                {
                  first_sub = first_sub->NZE_R_N;
                  if (first_sub)
                    sub_c_index = first_sub->c_index;
                  else
                    sub_c_index = Size;
                  l_sub++;
                }
              else if (sub_c_index > piv_c_index || l_sub >= nb_var_sub)
                {
                  int tmp_u_count = Get_u();
                  Insert(row, first_piv->c_index, tmp_u_count, 0);
                  u[tmp_u_count] = -u[first_piv->u_index]*first_elem;
                  first_piv = first_piv->NZE_R_N;
                  if (first_piv)
                    piv_c_index = first_piv->c_index;
                  else
                    piv_c_index = Size;
                  l_piv++;
                }
              else
                {
                  if (i == sub_c_index)
                    {
                      firsta = first;
                      first_suba = first_sub->NZE_R_N;
                      Delete(first_sub->r_index, first_sub->c_index);
                      first = firsta->NZE_C_N;
                      first_sub = first_suba;
                      if (first_sub)
                        sub_c_index = first_sub->c_index;
                      else
                        sub_c_index = Size;
                      l_sub++;
                      first_piv = first_piv->NZE_R_N;
                      if (first_piv)
                        piv_c_index = first_piv->c_index;
                      else
                        piv_c_index = Size;
                      l_piv++;
                    }
                  else
                    {
                      u[first_sub->u_index] -= u[first_piv->u_index]*first_elem;
                      first_sub = first_sub->NZE_R_N;
                      if (first_sub)
                        sub_c_index = first_sub->c_index;
                      else
                        sub_c_index = Size;
                      l_sub++;
                      first_piv = first_piv->NZE_R_N;
                      if (first_piv)
                        piv_c_index = first_piv->c_index;
                      else
                        piv_c_index = Size;
                      l_piv++;
                    }
                }
            }
          u[b[row]] -= u[b[pivj]]*first_elem;
        }
    }
  double slowc_lbx = slowc, res1bx;
  for (i = 0; i < y_size; i++)
    ya[i+it_*y_size] = y[i+it_*y_size];
  slowc_save = slowc;
  res1bx = simple_bksub(it_, Size, slowc_lbx);
  End(Size);
  mxFree(piv_v);
  mxFree(pivj_v);
  mxFree(pivk_v);
  mxFree(NR);
  mxFree(bc);
  return true;
}

void
SparseMatrix::CheckIt(int y_size, int y_kmin, int y_kmax, int Size, int periods, int iter)
{
  const double epsilon = 1e-7;
  fstream SaveResult;
  ostringstream out;
  out << "Result" << iter;
  SaveResult.open(out.str().c_str(), ios::in);
  if (!SaveResult.is_open())
    {
      mexPrintf("Error : Can't open file \"%s\" for reading\n", "Result");
      mexEvalString("st=fclose('all');clear all;");
      mexErrMsgTxt("Exit from Dynare");
    }
  mexPrintf("Reading Result...");
  int row, col;
  SaveResult >> row;
  mexPrintf("row=%d\n", row);
  SaveResult >> col;
  mexPrintf("col=%d\n", col);
  double G1a;
  mexPrintf("Allocated\n");
  NonZeroElem *first;
  for (int j = 0; j < col; j++)
    {
      mexPrintf("j=%d ", j);
      int nb_equ = At_Col(j, &first);
      mexPrintf("nb_equ=%d\n", nb_equ);
      int line;
      if (first)
        line = first->r_index;
      else
        line = -9999999;
      for (int i = 0; i < row; i++)
        {
          SaveResult >> G1a;
          if (line == i)
            {
              if (abs(u[first->u_index]/G1a-1) > epsilon)
                mexPrintf("Problem at r=%d c=%d u[first->u_index]=%5.14f G1a[i][j]=%5.14f %f\n", i, j, u[first->u_index], G1a, u[first->u_index]/G1a-1);
              first = first->NZE_C_N;
              if (first)
                line = first->r_index;
              else
                line = -9999999;
            }
          else
            {
              if (G1a != 0.0)
                mexPrintf("Problem at r=%d c=%d G1a[i][j]=%f\n", i, j, G1a);
            }
        }
    }
  mexPrintf("G1a red done\n");
  SaveResult >> row;
  mexPrintf("row(2)=%d\n", row);
  double *B;
  B = (double *) mxMalloc(row*sizeof(double));
  for (int i = 0; i < row; i++)
    SaveResult >> B[i];
  SaveResult.close();
  mexPrintf("done\n");
  mexPrintf("Comparing...");
  for (int i = 0; i < row; i++)
    {
      if (abs(u[b[i]]+B[i]) > epsilon)
        mexPrintf("Problem at i=%d u[b[i]]=%f B[i]=%f\n", i, u[b[i]], B[i]);
    }
  mxFree(B);
}

void
SparseMatrix::Check_the_Solution(int periods, int y_kmin, int y_kmax, int Size, double *u, int *pivot, int *b)
{
  const double epsilon = 1e-10;
  Init(periods, y_kmin, y_kmax, Size, IM_i);
  NonZeroElem *first;
  int cal_y = y_kmin*Size;
  mexPrintf("     ");
  for (int i = 0; i < Size; i++)
    mexPrintf(" %8d", i);
  mexPrintf("\n");
  for (int t = y_kmin; t < periods+y_kmin; t++)
    {
      mexPrintf("t=%5d", t);
      for (int i = 0; i < Size; i++)
        mexPrintf(" %d %1.6f", t*y_size+index_vara[i], y[t*y_size+index_vara[i]]);
      mexPrintf("\n");
    }
  for (int i = 0; i < Size*periods; i++)
    {
      double res = 0;
      int pos = pivot[i];
      mexPrintf("pos[%d]=%d", i, pos);
      int nb_var = At_Row(pos, &first);
      mexPrintf(" nb_var=%d\n", nb_var);
      for (int j = 0; j < nb_var; j++)
        {
          mexPrintf("(y[%d]=%f)*(u[%d]=%f)(r=%d, c=%d)\n", index_vara[first->c_index]+cal_y, y[index_vara[first->c_index]+cal_y], first->u_index, u[first->u_index], first->r_index, first->c_index);
          res += y[index_vara[first->c_index]+cal_y]*u[first->u_index];
          first = first->NZE_R_N;
        }
      double tmp_ = res;
      res += u[b[pos]];
      if (abs(res) > epsilon)
        mexPrintf("Error for equation %d => res=%f y[%d]=%f u[b[%d]]=%f somme(y*u)=%f\n", pos, res, pos, y[index_vara[pos]], pos, u[b[pos]], tmp_);
    }
  filename += " stopped";
  mexErrMsgTxt(filename.c_str());
}

int
SparseMatrix::simulate_NG1(int blck, int y_size, int it_, int y_kmin, int y_kmax, int Size, int periods, bool print_it, bool cvg, int &iter, int minimal_solving_periods, int Block_number)
{
  /*Triangularisation at each period of a block using a simple gaussian Elimination*/
  t_save_op_s *save_op_s;
  bool record = false;
  int *save_op = NULL, *save_opa = NULL, *save_opaa = NULL;
  long int nop = 0, nopa = 0;
  int tbreak = 0, last_period = periods;
  int i, j, k;
  int pivj = 0, pivk = 0;
  int tmp_u_count, lag;
  NonZeroElem *first;
  double piv_abs;
  if (start_compare == 0)
    start_compare = y_kmin;
  u_count_alloc_save = u_count_alloc;
  clock_t t01;
  clock_t t1 = clock();
  nop1 = 0;
  error_not_printed = true;
  if (iter > 0)
    {
      mexPrintf("Sim : %f ms\n", (1000.0*(double (clock())-double (time00)))/double (CLOCKS_PER_SEC));
      mexEvalString("drawnow;");
      time00 = clock();
    }
  if (isnan(res1) || isinf(res1) || (res2 > 12*g0 && iter>0))
    {
      if (iter == 0)
        {
          for (j = 0; j < y_size; j++)
            {
              bool select = false;
              for (int i = 0; i < Size; i++)
                if (j == index_vara[i])
                  {
                    select = true;
                    break;
                  }
              if (select)
                mexPrintf("-> variable %d at time %d = %f direction = %f\n", j+1, it_, y[j+it_*y_size], direction[j+it_*y_size]);
              else
                mexPrintf("   variable %d at time %d = %f direction = %f\n", j+1, it_, y[j+it_*y_size], direction[j+it_*y_size]);
            }
          mexPrintf("res1=%5.10f\n", res1);
          mexPrintf("The initial values of endogenous variables are too far from the solution.\n");
          mexPrintf("Change them!\n");
          mexEvalString("drawnow;");
          filename += " stopped";
          mexErrMsgTxt(filename.c_str());
        }
      if (fabs(slowc_save) < 1e-8)
        {
          for (j = 0; j < y_size; j++)
            {
              bool select = false;
              for (int i = 0; i < Size; i++)
                if (j == index_vara[i])
                  {
                    select = true;
                    break;
                  }
              if (select)
                mexPrintf("-> variable %d at time %d = %f direction = %f\n", j+1, it_, y[j+it_*y_size], direction[j+it_*y_size]);
              else
                mexPrintf("   variable %d at time %d = %f direction = %f\n", j+1, it_, y[j+it_*y_size], direction[j+it_*y_size]);
            }
          mexPrintf("Dynare cannot improve the simulation in block %d at time %d (variable %d)\n", blck+1, it_+1, max_res_idx);
          mexEvalString("drawnow;");
          filename += " stopped";
          mexErrMsgTxt(filename.c_str());
        }
      if(!(isnan(res1) || isinf(res1)) && !(isnan(g0) || isinf(g0)))
        {
          if (try_at_iteration == 0)
            {
              prev_slowc_save = slowc_save;
              slowc_save = max( - gp0 / (2 * (res2 - g0 - gp0)) , 0.1);
            }
          else
            {
              double t1 = res2 - gp0 * slowc_save - g0;
              double t2 = glambda2 - gp0 * prev_slowc_save - g0;
              double a = (1/(slowc_save * slowc_save) * t1 - 1/(prev_slowc_save * prev_slowc_save) * t2) / (slowc_save - prev_slowc_save);
              double b = (-prev_slowc_save/(slowc_save * slowc_save) * t1 + slowc_save/(prev_slowc_save * prev_slowc_save) * t2) / (slowc_save - prev_slowc_save);
              prev_slowc_save = slowc_save;
              slowc_save = max(min( -b + sqrt(b*b - 3 * a * gp0) / (3 * a), 0.5 * slowc_save), 0.1 * slowc_save);
            }
          glambda2 = res2;
          try_at_iteration ++;
        }
      else
        {
          prev_slowc_save = slowc_save;
          slowc_save /= 1.1;
        }
      if (print_it)
        mexPrintf("Error: Simulation diverging, trying to correct it using slowc=%f\n", slowc_save);
      for (i = 0; i < y_size*(periods+y_kmin); i++)
        y[i] = ya[i]+slowc_save*direction[i];
      iter--;
      return 0;
    }
  /*if (isnan(res1) || isinf(res1))
    {
      if (iter == 0)
        {
          for (j = 0; j < y_size; j++)
            mexPrintf("variable %d at time %d = %f\n", j+1, it_, y[j+it_*y_size]);
          for (j = 0; j < Size; j++)
            mexPrintf("residual(%d)=%5.25f\n", j, u[j]);
          mexPrintf("res1=%5.25f\n", res1);
          mexPrintf("The initial values of endogenous variables are too far from the solution.\n");
          mexPrintf("Change them!\n");
          mexEvalString("drawnow;");
          mexEvalString("st=fclose('all');clear all;");
          filename += " stopped";
          mexErrMsgTxt(filename.c_str());
        }
      if (slowc_save < 1e-8)
        {
          mexPrintf("slowc_save=%g\n", slowc_save);
          for (j = 0; j < y_size; j++)
            mexPrintf("variable %d at time %d = %f direction = %f variable at last step = %f b = %f\n", j+1, it_+1, y[j+it_*y_size], direction[j+it_*y_size], ya[j+it_*y_size], u[pivot[j+it_*y_size]]);
          mexPrintf("Dynare cannot improve the simulation in block %d at time %d (variable %d max_res = %f, res1 = %f)\n", blck+1, it_+1, max_res_idx, max_res, res1);
          mexEvalString("drawnow;");
          mexEvalString("st=fclose('all');clear all;");
          filename += " stopped";
          mexErrMsgTxt(filename.c_str());
        }
      slowc_save /= 2;
      mexPrintf("Error: Simulation diverging, trying to correct it using slowc=%f\n", slowc_save);
      for (i = 0; i < y_size*(periods+y_kmin); i++)
        y[i] = ya[i]+slowc_save*direction[i];
      iter--;
      return (0);
    }*/
  u_count += u_count_init;
  if (alt_symbolic && alt_symbolic_count < alt_symbolic_count_max)
    {
      mexPrintf("Pivoting method will be applied only to the first periods.\n");
      alt_symbolic = false;
      symbolic = true;
      markowitz_c = markowitz_c_s;
      alt_symbolic_count++;
    }
  if (((res1/res1a-1) > -0.3) && symbolic && iter > 0)
    {
      if (restart > 2)
        {
          mexPrintf("Divergence or slowdown occured during simulation.\nIn the next iteration, pivoting method will be applied to all periods.\n");
          symbolic = false;
          alt_symbolic = true;
          markowitz_c_s = markowitz_c;
          markowitz_c = 0;
        }
      else
        {
          mexPrintf("Divergence or slowdown occured during simulation.\nIn the next iteration, pivoting method will be applied for a longer period.\n");
          start_compare = min(tbreak_g, periods);
          restart++;
        }
    }
  else
    {
      start_compare = max(y_kmin, minimal_solving_periods);
      restart = 0;
    }
  res1a = res1;

  if (print_it)
    {
      mexPrintf("-----------------------------------\n");
      mexPrintf("      Simulate iteration no %d     \n", iter+1);
      mexPrintf("      max. error=%.10e       \n", double (max_res));
      mexPrintf("      sqr. error=%.10e       \n", double (res2));
      mexPrintf("      abs. error=%.10e       \n", double (res1));
      mexPrintf("-----------------------------------\n");
    }
  if (cvg)
    {
      return (0);
    }
  else
    {
      Init(periods, y_kmin, y_kmax, Size, IM_i);
      double *piv_v;
      int *pivj_v, *pivk_v, *NR;
      piv_v = (double *) mxMalloc(Size*sizeof(double));
      pivj_v = (int *) mxMalloc(Size*sizeof(int));
      pivk_v = (int *) mxMalloc(Size*sizeof(int));
      NR = (int *) mxMalloc(Size*sizeof(int));
      for (int t = 0; t < periods; t++)
        {
          if (record && symbolic)
            {
              if (save_op)
                {
                  mxFree(save_op);
                  save_op = NULL;
                }
              save_op = (int *) mxMalloc(nop*sizeof(int));
              nopa = nop;
            }
          nop = 0;
          Clear_u();
          int ti = t*Size;
          for (i = ti; i < Size+ti; i++)
            {
              /*finding the max-pivot*/
              double piv = piv_abs = 0;
              int nb_eq = At_Col(i, 0, &first);
              if ((symbolic && t <= start_compare) || !symbolic)
                {
                  int l = 0, N_max = 0;
                  bool one = false;
                  piv_abs = 0;
                  for (j = 0; j < nb_eq; j++)
                    {
                      if (!line_done[first->r_index])
                        {
                          k = first->u_index;
                          int jj = first->r_index;
                          int NRow_jj = NRow(jj);
                          piv_v[l] = u[k];
                          double piv_fabs = fabs(u[k]);
                          pivj_v[l] = jj;
                          pivk_v[l] = k;
                          NR[l] = NRow_jj;
                          if (NRow_jj == 1 && !one)
                            {
                              one = true;
                              piv_abs = piv_fabs;
                              N_max = NRow_jj;
                            }
                          if (!one)
                            {
                              if (piv_fabs > piv_abs)
                                piv_abs = piv_fabs;
                              if (NRow_jj > N_max)
                                N_max = NRow_jj;
                            }
                          else
                            {
                              if (NRow_jj == 1)
                                {
                                  if (piv_fabs > piv_abs)
                                    piv_abs = piv_fabs;
                                  if (NRow_jj > N_max)
                                    N_max = NRow_jj;
                                }
                            }
                          l++;
                        }
                      first = first->NZE_C_N;
                    }
                  double markovitz = 0, markovitz_max = -9e70;
                  int NR_max = 0;
                  if (!one)
                    {
                      for (j = 0; j < l; j++)
                        {
                          if (N_max > 0 && NR[j] > 0)
                            {
                              if (fabs(piv_v[j]) > 0)
                                {
                                  if (markowitz_c > 0)
                                    markovitz = exp(log(fabs(piv_v[j])/piv_abs)-markowitz_c*log(double (NR[j])/double (N_max)));
                                  else
                                    markovitz = fabs(piv_v[j])/piv_abs;
                                }
                              else
                                markovitz = 0;
                            }
                          else
                            markovitz = fabs(piv_v[j])/piv_abs;
                          if (markovitz > markovitz_max)
                            {
                              piv = piv_v[j];
                              pivj = pivj_v[j];   //Line number
                              pivk = pivk_v[j];   //positi
                              markovitz_max = markovitz;
                              NR_max = NR[j];
                            }
                        }
                    }
                  else
                    {
                      for (j = 0; j < l; j++)
                        {
                          if (N_max > 0 && NR[j] > 0)
                            {
                              if (fabs(piv_v[j]) > 0)
                                {
                                  if (markowitz_c > 0)
                                    markovitz = exp(log(fabs(piv_v[j])/piv_abs)-markowitz_c*log(double (NR[j])/double (N_max)));
                                  else
                                    markovitz = fabs(piv_v[j])/piv_abs;
                                }
                              else
                                markovitz = 0;
                            }
                          else
                            markovitz = fabs(piv_v[j])/piv_abs;
                          if (NR[j] == 1)
                            {
                              piv = piv_v[j];
                              pivj = pivj_v[j];   //Line number
                              pivk = pivk_v[j];   //positi
                              markovitz_max = markovitz;
                              NR_max = NR[j];
                            }
                        }
                    }
                  if (fabs(piv) < eps)
                    mexPrintf("==> Error NR_max=%d, N_max=%d and piv=%f, piv_abs=%f, markovitz_max=%f\n",NR_max, N_max, piv, piv_abs, markovitz_max);
                  if (NR_max == 0)
                    mexPrintf("==> Error NR_max=0 and piv=%f, markovitz_max=%f\n",piv, markovitz_max);
                  pivot[i] = pivj;
                  pivot_save[i] = pivj;
                  pivotk[i] = pivk;
                  pivotv[i] = piv;
                }
              else
                {
                  pivj = pivot[i-Size]+Size;
                  pivot[i] = pivj;
                  At_Pos(pivj, i, &first);
                  pivk = first->u_index;
                  piv = u[pivk];
                  piv_abs = fabs(piv);
                }
              line_done[pivj] = true;
              if (symbolic)
                {
                  if (record)
                    {
                      if (nop+1 >= nopa)
                        {
                          nopa = int (1.5*nopa);
                          save_op = (int *) mxRealloc(save_op, nopa*sizeof(int));
                        }
                      save_op_s = (t_save_op_s *)(&(save_op[nop]));
                      save_op_s->operat = IFLD;
                      save_op_s->first = pivk;
                      save_op_s->lag = 0;
                    }
                  nop += 2;
                }
              if (piv_abs < eps)
                {
                  if (Block_number>1)
                    mexPrintf("Error: singular system in Simulate_NG1 in block %d\n", blck);
                  else
                    mexPrintf("Error: singular system in Simulate_NG1\n");
                  mexEvalString("drawnow;");
                  filename += " stopped";
                  mexErrMsgTxt(filename.c_str());
                }
              /*divide all the non zeros elements of the line pivj by the max_pivot*/
              int nb_var = At_Row(pivj, &first);
              NonZeroElem **bb;
              bb = (NonZeroElem **) mxMalloc(nb_var*sizeof(first));
              for (j = 0; j < nb_var; j++)
                {
                  bb[j] = first;
                  first = first->NZE_R_N;
                }

              for (j = 0; j < nb_var; j++)
                {
                  first = bb[j];
                  u[first->u_index] /= piv;
                  if (symbolic)
                    {
                      if (record)
                        {
                          if (nop+j*2+1 >= nopa)
                            {
                              nopa = int (1.5*nopa);
                              save_op = (int *) mxRealloc(save_op, nopa*sizeof(int));
                            }
                          save_op_s = (t_save_op_s *)(&(save_op[nop+j*2]));
                          save_op_s->operat = IFDIV;
                          save_op_s->first = first->u_index;
                          save_op_s->lag = first->lag_index;
                        }
                    }
                }
              mxFree(bb);
              nop += nb_var*2;
              u[b[pivj]] /= piv;
              if (symbolic)
                {
                  if (record)
                    {
                      if (nop+1 >= nopa)
                        {
                          nopa = int (1.5*nopa);
                          save_op = (int *) mxRealloc(save_op, nopa*sizeof(int));
                        }
                      save_op_s = (t_save_op_s *)(&(save_op[nop]));
                      save_op_s->operat = IFDIV;
                      save_op_s->first = b[pivj];
                      save_op_s->lag = 0;
                    }
                  nop += 2;
                }
              /*substract the elements on the non treated lines*/
              nb_eq = At_Col(i, &first);
              NonZeroElem *first_piva;
              int nb_var_piva = At_Row(pivj, &first_piva);

              NonZeroElem **bc;
              bc = (NonZeroElem **) mxMalloc(nb_eq*sizeof(first));
              int nb_eq_todo = 0;
              for (j = 0; j < nb_eq && first; j++)
                {
                  if (!line_done[first->r_index])
                    bc[nb_eq_todo++] = first;
                  first = first->NZE_C_N;
                }
              //#pragma omp parallel for num_threads(2) shared(nb_var_piva, first_piva, nopa, nop, save_op, record)
              for (j = 0; j < nb_eq_todo; j++)
                {
                  t_save_op_s *save_op_s_l;
                  first = bc[j];
                  int row = first->r_index;
                  double first_elem = u[first->u_index];
                  if (symbolic)
                    {
                      if (record)
                        {
                          if (nop+2 >= nopa)
                            {
                              //#pragma omp critical
                              {
                                nopa = int (1.5*nopa);
                                save_op = (int *) mxRealloc(save_op, nopa*sizeof(int));
                              }
                            }
                          save_op_s_l = (t_save_op_s *)(&(save_op[nop]));
                          save_op_s_l->operat = IFLD;
                          save_op_s_l->first = first->u_index;
                          save_op_s_l->lag = abs(first->lag_index);
                        }
                      nop += 2;
                    }

                  int nb_var_piv = nb_var_piva;
                  NonZeroElem *first_piv = first_piva;
                  NonZeroElem *first_sub;
                  int nb_var_sub = At_Row(row, &first_sub);
                  int l_sub = 0;
                  int l_piv = 0;
                  int sub_c_index = first_sub->c_index;
                  int piv_c_index = first_piv->c_index;
                  int tmp_lag = first_sub->lag_index;
                  while (l_sub < nb_var_sub || l_piv < nb_var_piv)
                    {
                      if (l_sub < nb_var_sub && (sub_c_index < piv_c_index || l_piv >= nb_var_piv))
                        {
                          //There is no nonzero element at row pivot for this column=> Nothing to do for the current element got to next column
                          first_sub = first_sub->NZE_R_N;
                          if (first_sub)
                            sub_c_index = first_sub->c_index;
                          else
                            sub_c_index = Size*periods;
                          l_sub++;
                        }
                      else if (sub_c_index > piv_c_index || l_sub >= nb_var_sub)
                        {
                          // There is an nonzero element at row pivot but not at the current row=> insert a negative element in the current row
                          tmp_u_count = Get_u();
                          lag = first_piv->c_index/Size-row/Size;
                          //#pragma omp critical
                          {
                            Insert(row, first_piv->c_index, tmp_u_count, lag);
                          }
                          u[tmp_u_count] = -u[first_piv->u_index]*first_elem;
                          if (symbolic)
                            {
                              if (record)
                                {
                                  if (nop+2 >= nopa)
                                    {
                                      //#pragma omp critical
                                      {
                                        nopa = int (1.5*nopa);
                                        save_op = (int *) mxRealloc(save_op, nopa*sizeof(int));
                                      }
                                    }
                                  save_op_s_l = (t_save_op_s *)(&(save_op[nop]));
                                  save_op_s_l->operat = IFLESS;
                                  save_op_s_l->first = tmp_u_count;
                                  save_op_s_l->second = first_piv->u_index;
                                  save_op_s_l->lag = max(first_piv->lag_index, abs(tmp_lag));
                                }
                              nop += 3;
                            }
                          first_piv = first_piv->NZE_R_N;
                          if (first_piv)
                            piv_c_index = first_piv->c_index;
                          else
                            piv_c_index = Size*periods;
                          l_piv++;
                        }
                      else /*first_sub->c_index==first_piv->c_index*/
                        {
                          if (i == sub_c_index)
                            {

                              //#pragma omp barrier
                              //#pragma omp single
                              //#pragma omp critical
                              {
                                NonZeroElem *firsta = first;
                                NonZeroElem *first_suba = first_sub->NZE_R_N;
                                Delete(first_sub->r_index, first_sub->c_index);
                                first = firsta->NZE_C_N;
                                first_sub = first_suba;
                              }

                              if (first_sub)
                                sub_c_index = first_sub->c_index;
                              else
                                sub_c_index = Size*periods;
                              l_sub++;
                              first_piv = first_piv->NZE_R_N;
                              if (first_piv)
                                piv_c_index = first_piv->c_index;
                              else
                                piv_c_index = Size*periods;
                              l_piv++;
                            }
                          else
                            {
                              u[first_sub->u_index] -= u[first_piv->u_index]*first_elem;
                              if (symbolic)
                                {
                                  if (record)
                                    {
                                      if (nop+3 >= nopa)
                                        {
                                          //#pragma omp critical
                                          {
                                            nopa = int (1.5*nopa);
                                            save_op = (int *) mxRealloc(save_op, nopa*sizeof(int));
                                          }
                                        }
                                      save_op_s_l = (t_save_op_s *)(&(save_op[nop]));
                                      save_op_s_l->operat = IFSUB;
                                      save_op_s_l->first = first_sub->u_index;
                                      save_op_s_l->second = first_piv->u_index;
                                      save_op_s_l->lag = max(abs(tmp_lag), first_piv->lag_index);
                                    }
                                  nop += 3;
                                }
                              first_sub = first_sub->NZE_R_N;
                              if (first_sub)
                                sub_c_index = first_sub->c_index;
                              else
                                sub_c_index = Size*periods;
                              l_sub++;
                              first_piv = first_piv->NZE_R_N;
                              if (first_piv)
                                piv_c_index = first_piv->c_index;
                              else
                                piv_c_index = Size*periods;
                              l_piv++;
                            }
                        }
                    }
                  u[b[row]] -= u[b[pivj]]*first_elem;

                  if (symbolic)
                    {
                      if (record)
                        {
                          if (nop+3 >= nopa)
                            {
                              //#pragma omp critical
                              {
                                nopa = int (1.5*nopa);
                                save_op = (int *) mxRealloc(save_op, nopa*sizeof(int));
                              }
                            }
                          save_op_s_l = (t_save_op_s *)(&(save_op[nop]));
                          save_op_s_l->operat = IFSUB;
                          save_op_s_l->first = b[row];
                          save_op_s_l->second = b[pivj];
                          save_op_s_l->lag = abs(tmp_lag);
                        }
                      nop += 3;
                    }
                }
              mxFree(bc);
            }
          if (symbolic)
            {
              if (record && (nop == nop1))
                {
                  if (save_opa && save_opaa)
                    {
                      if (compare(save_op, save_opa, save_opaa, t, periods, nop, Size))
                        {
                          tbreak = t;
                          tbreak_g = tbreak;
                          break;
                        }
                    }
                  if (save_opa)
                    {
                      if (save_opaa)
                        {
                          mxFree(save_opaa);
                          save_opaa = NULL;
                        }
                      save_opaa = (int *) mxMalloc(nop1*sizeof(int));
                      memcpy(save_opaa, save_opa, nop1*sizeof(int));
                    }
                  if (save_opa)
                    {
                      mxFree(save_opa);
                      save_opa = NULL;
                    }
                  save_opa = (int *) mxMalloc(nop*sizeof(int));
                  memcpy(save_opa, save_op, nop*sizeof(int));
                }
              else
                {
                  if (nop == nop1)
                    record = true;
                  else
                    {
                      record = false;
                      if (save_opa)
                        {
                          mxFree(save_opa);
                          save_opa = NULL;
                        }
                      if (save_opaa)
                        {
                          mxFree(save_opaa);
                          save_opaa = NULL;
                        }
                    }
                }
              nop2 = nop1;
              nop1 = nop;
            }
        }
      mxFree(piv_v);
      mxFree(pivj_v);
      mxFree(pivk_v);
      mxFree(NR);
    }
  nop_all += nop;
  if (symbolic)
    {
      if (save_op)
        mxFree(save_op);
      if (save_opa)
        mxFree(save_opa);
      if (save_opaa)
        mxFree(save_opaa);
    }

  /*The backward substitution*/
  double slowc_lbx = slowc, res1bx;
  for (i = 0; i < y_size*(periods+y_kmin); i++)
    ya[i] = y[i];
  slowc_save = slowc;
  res1bx = bksub(tbreak, last_period, Size, slowc_lbx);
  t01 = clock();
  End(Size);
  if (print_it)
    {
      clock_t t2 = clock();
      mexPrintf("(** %f milliseconds **)\n", 1000.0*(double (t2) - double (t1))/double (CLOCKS_PER_SEC));
      mexEvalString("drawnow;");
    }

  time00 = clock();
  if (tbreak_g == 0)
    tbreak_g = periods;

  return (0);
}

void
SparseMatrix::fixe_u(double **u, int u_count_int, int max_lag_plus_max_lead_plus_1)
{
  u_count = u_count_int * periods;
  u_count_alloc = 2*u_count;
  (*u) = (double *) mxMalloc(u_count_alloc*sizeof(double));
  memset((*u), 0, u_count_alloc*sizeof(double));
  u_count_init = max_lag_plus_max_lead_plus_1;
}

