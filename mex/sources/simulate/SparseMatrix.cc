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

#include <cstring>
#include <sstream>
#include "SparseMatrix.hh"

SparseMatrix::SparseMatrix()
{
  pivotva=NULL;
  g_save_op=NULL;
  g_nop_all=0;
  mem_mngr.init_Mem();
  symbolic=true;
  alt_symbolic=false;
  alt_symbolic_count=0;
  max_u=0;
  min_u=0x7FFFFFFF;
  res1a=9.0e60;
  tbreak_g=0;
  start_compare=0;
  restart = 0;
}


//


int SparseMatrix::NRow(int r)
{
  return NbNZRow[r];
}

int SparseMatrix::NCol(int c)
{
  return NbNZCol[c];
}


int SparseMatrix::At_Row(int r, NonZeroElem **first)
{
  (*first)=FNZE_R[r];
  return NbNZRow[r];
}

int
SparseMatrix::Union_Row(int row1, int row2)
{
  NonZeroElem *first1, *first2;
  int n1=At_Row(row1, &first1);
  int n2=At_Row(row2, &first2);
  int i1=0, i2=0, nb_elem=0;
  while (i1<n1 && i2<n2)
    {
      if (first1->c_index==first2->c_index)
        {
          nb_elem++;
          i1++;
          i2++;
          first1=first1->NZE_R_N;
          first2=first2->NZE_R_N;
        }
      else if (first1->c_index<first2->c_index)
        {
          nb_elem++;
          i1++;
          first1=first1->NZE_R_N;
        }
      else
        {
          nb_elem++;
          i2++;
          first2=first2->NZE_R_N;
        }
    }
  return nb_elem;
}

int
SparseMatrix::At_Pos(int r, int c, NonZeroElem **first)
{
  (*first)=FNZE_R[r];
  while ((*first)->c_index!=c /*&& (*first)->NZE_R_N*/)
    {
#ifdef PRINT_OUT
      mexPrintf("looking not CRS [%d, %d]\n",(*first)->r_index,(*first)->c_index);
#endif
      (*first)=(*first)->NZE_R_N;
    }
  /*if ((*first)->c_index!=c)
    mexPrintf("-----------------------  cannot find M[%d, %d]\n",r,c);*/
  return NbNZRow[r];
}


int SparseMatrix::At_Col(int c, NonZeroElem **first)
{
  (*first)=FNZE_C[c];
  return NbNZCol[c];
}

int SparseMatrix::At_Col(int c, int lag, NonZeroElem **first)
{
  (*first)=FNZE_C[c];
  int i=0;
  while ((*first)->lag_index!=lag && (*first))
    {
#ifdef PRINT_OUT
      mexPrintf("first->lag_index(%d) != %d\n",(*first)->lag_index,lag);
#endif
      (*first)=(*first)->NZE_C_N;
    }
  if ((*first))
    {
#ifdef PRINT_OUT
      mexPrintf("first=%x\n",(*first));
#endif
      NonZeroElem* firsta=(*first);
      if (!firsta->NZE_C_N)
        i++;
      else
        {
          while (firsta->lag_index==lag && firsta->NZE_C_N)
            {
#ifdef PRINT_OUT
              mexPrintf("firsta->lag_index(%d) == %d, eq=%d, var=%d\n",firsta->lag_index,lag, firsta->r_index, firsta->c_index);
#endif
              firsta=firsta->NZE_C_N;
              i++;
            }
          if (firsta->lag_index==lag) i++;
        }
    }
#ifdef PRINT_OUT
  mexPrintf("i=%d\n",i);
#endif
  return i;
}

#ifdef PROFILER
double tdelete1=0, tdelete2=0, tdelete21=0, tdelete22=0, tdelete221=0, tdelete222=0, tcompare=0;
#endif

void SparseMatrix::Delete(const int r,const int c, const int Size)
{
	//mexPrintf("Delete r=%d c=%d\n",r,c);
  NonZeroElem *first=FNZE_R[r], *firsta=NULL;
#ifdef PROFILER
  clock_t td0, td1, td2;
  td0=clock();
#endif
  while (first->c_index!=c)
    {
      firsta=first;
      first=first->NZE_R_N;
    }
#ifdef PRINT_OUT
      mexPrintf("CRS [%d, %d]=c(%d)\n",first->r_index,first->c_index,c);
      mexEvalString("drawnow;");
#endif
  if (firsta!=NULL)
    firsta->NZE_R_N=first->NZE_R_N;
  if (first==FNZE_R[r])
    FNZE_R[r]=first->NZE_R_N;
  NbNZRow[r]--;
#ifdef PROFILER
  tdelete1+=clock()-td0;
  td0=clock();
  td1=clock();
#endif
  first=FNZE_C[c];
  firsta=NULL;
  while (first->r_index!=r)
    {
      firsta=first;
      first=first->NZE_C_N;
    }
#ifdef PRINT_OUT
  mexPrintf("CSS [%d, %d]=r(%d)\n",first->r_index,first->c_index,r);
  mexEvalString("drawnow;");
#endif
#ifdef PROFILER
  tdelete21+=clock()-td1;
  td1=clock();
#endif
  if (firsta!=NULL)
    firsta->NZE_C_N=first->NZE_C_N;
  if (first==FNZE_C[c])
    FNZE_C[c]=first->NZE_C_N;
#ifdef PROFILER
  td2=clock();
#endif
  u_liste.push_back(first->u_index);
#ifdef PROFILER
  tdelete221+=clock()-td2;
  td2=clock();
#endif
#ifdef NEW_ALLOC
  mem_mngr.mxFree_NZE(first);
#else
  mxFree(first);
#endif
  NbNZCol[c]--;
#ifdef PROFILER
  tdelete222+=clock()-td2;
#endif
#ifdef PROFILER
  tdelete22+=clock()-td1;
  tdelete2+=clock()-td0;
#endif
  /*Check the deletition*/
  /*int nb_var=NbNZRow[r];
  first=FNZE_R[r];
  for(int j=0;j<nb_var;j++)
    {
    	if(!first)
    	  mexPrintf("Error in Delete (Row) r=%d and c=%d \n",r,c);
      first=first->NZE_R_N;
    }
	nb_var=NbNZCol[c];
  first=FNZE_C[c];
  for(int j=0;j<nb_var;j++)
    {
    	if(!first)
    	  mexPrintf("Error in Delete (Col) r=%d and c=%d \n",r,c);
      first=first->NZE_C_N;
    }*/
}


void SparseMatrix::Print(int Size, int *b)
{
  int a,i,j,k,l;
  mexPrintf("   ");
  for (k=0;k<Size*periods;k++)
    mexPrintf("%-2d ",k);
  mexPrintf("    |    ");
  for (k=0;k<Size*periods;k++)
    mexPrintf("%8d",k);
  mexPrintf("\n");
  for (i=0;i<Size*periods;i++)
    {
      NonZeroElem *first=FNZE_R[i];
      j=NbNZRow[i];
      mexPrintf("%-2d ",i);
      a=0;
      for (k=0;k<j;k++)
        {
          for (l=0;l<(first->c_index-a);l++)
            mexPrintf("   ");
          mexPrintf("%-2d ",first->u_index);
          a=first->c_index+1;
          first=first->NZE_R_N;
        }
      for (k=a;k<Size*periods;k++)
        mexPrintf("   ");
      mexPrintf("%-2d ",b[i]);

      first=FNZE_R[i];
      j=NbNZRow[i];
      mexPrintf(" | %-2d ",i);
      a=0;
      for (k=0;k<j;k++)
        {
          for (l=0;l<(first->c_index-a);l++)
            mexPrintf("        ");
          mexPrintf("%8.4f",double(u[first->u_index]));
          a=first->c_index+1;
          first=first->NZE_R_N;
        }
      for (k=a;k<Size*periods;k++)
        mexPrintf("        ");
      mexPrintf("%8.4f",double(u[b[i]]));
      mexPrintf("\n");
    }
}



void SparseMatrix::Insert(const int r, const int c, const int u_index, const int lag_index)
{
	//mexPrintf("Insert r=%d c=%d\n",r,c);
#ifdef PRINT_OUT
  mexPrintf("In Insert r=%d, c=%d, u=%d, lag=%d \n",r,c,u_index,lag_index);
#endif
  NonZeroElem *firstn, *first, *firsta;
  /*if (first)
    {*/
#ifdef NEW_ALLOC
  firstn=mem_mngr.mxMalloc_NZE();
#else
  firstn=(NonZeroElem*)mxMalloc(sizeof(NonZeroElem));
#endif
  first=FNZE_R[r];
  firsta=NULL;
#ifdef PRINT_OUT
  mexPrintf("first->c_index=%d, first->NZE_R_N=%x\n",first->c_index, first->NZE_R_N);
#endif
  while (first->c_index<c && first->NZE_R_N)
    {
      firsta=first;
#ifdef PRINT_OUT
      mexPrintf("drop first->c_index=%d c=%d\n",first->c_index,c);
#endif
      first=first->NZE_R_N;
    }
#ifdef PRINT_OUT
  mexPrintf("retain first->c_index=%d c=%d\n",first->c_index,c);
#endif
  firstn->u_index=u_index;
  firstn->r_index=r;
  firstn->c_index=c;
  firstn->lag_index=lag_index;
  if (first->c_index>c)
    {
      if (first==FNZE_R[r])
        FNZE_R[r]=firstn;
      if (firsta!=NULL)
        firsta->NZE_R_N=firstn;
      firstn->NZE_R_N=first;
    }
  else /*first.c_index<c*/
    {
    	/*if(first->c_index==c)
    	  mexPrintf("Error in Insert (r=%d, c=%d -Row-) already exist!!\n");*/
      first->NZE_R_N=firstn;
      firstn->NZE_R_N=NULL;
    }
  NbNZRow[r]++;
  first=FNZE_C[c];
  firsta=NULL;
  while (first->r_index<r && first->NZE_C_N)
    {
      firsta=first;
      first=first->NZE_C_N;
    }
  if (first->r_index>r)
    {
      if (first==FNZE_C[c])
        FNZE_C[c]=firstn;
      if (firsta!=NULL)
        firsta->NZE_C_N=firstn;
      firstn->NZE_C_N=first;
    }
  else /*first.r_index<r*/
    {
    	/*if(first->r_index==r)
    	  mexPrintf("Error in Insert (r=%d, c=%d -Col-) already exist!!\n");*/
      first->NZE_C_N=firstn;
      firstn->NZE_C_N=NULL;
    }
  NbNZCol[c]++;
  /*Check the insertion*/
  /*int nb_var=NbNZRow[r];
  first=FNZE_R[r];
  for(int j=0;j<nb_var;j++)
    {
    	if(!first)
    	  mexPrintf("Error in insert (Row) r=%d and c=%d \n",r,c);
      first=first->NZE_R_N;
    }
	nb_var=NbNZCol[c];
  first=FNZE_C[c];
  for(int j=0;j<nb_var;j++)
    {
    	if(!first)
    	  mexPrintf("Error in insert (Col) r=%d and c=%d \n",r,c);
      first=first->NZE_C_N;
    }*/
}

void SparseMatrix::Read_SparseMatrix(string file_name, int Size, int periods, int y_kmin, int y_kmax)
{
  int i,j,eq,var,lag;
  filename=file_name;
  mem_mngr.fixe_file_name(file_name);
  if (!SaveCode.is_open())
    {
#ifdef PRINT_OUT
      mexPrintf("file opened\n");
#endif
      SaveCode.open((file_name + ".bin").c_str(), std::ios::in | std::ios::binary);
      if (!SaveCode.is_open())
        {
          mexPrintf("Error : Can't open file \"%s\" for reading\n", (file_name + ".bin").c_str());
          mexEvalString("st=fclose('all');clear all;");
          mexErrMsgTxt("Exit from Dynare");
        }
#ifdef PRINT_OUT
      mexPrintf("done\n");
#endif
    }
  IM_i.clear();
  //mexPrintf("u_count_init=%d\n",u_count_init);
  for (i=0;i<u_count_init;i++)
    {
      SaveCode.read(reinterpret_cast<char *>(&eq), sizeof(eq));
      SaveCode.read(reinterpret_cast<char *>(&var), sizeof(var));
      SaveCode.read(reinterpret_cast<char *>(&lag), sizeof(lag));
      SaveCode.read(reinterpret_cast<char *>(&j), sizeof(j));
      //mexPrintf("eq=%d var=%d lag=%d j=%d\n",eq, var, lag, j);
      IM_i[std::make_pair(std::make_pair(eq, var), lag)] = j;
    }
#ifdef MEM_ALLOC_CHK
  mexPrintf("index_vara=(int*)mxMalloc(%d*sizeof(int))\n",Size*(periods+y_kmin+y_kmax));
#endif
  index_vara=(int*)mxMalloc(Size*(periods+y_kmin+y_kmax)*sizeof(int));
#ifdef MEM_ALLOC_CHK
  mexPrintf("ok\n");
#endif
  for (j=0;j<Size;j++)
    {
      SaveCode.read(reinterpret_cast<char *>(&index_vara[j]), sizeof(*index_vara));
    }
  if(periods+y_kmin+y_kmax>1)
    {
      for (i=1;i<periods+y_kmin+y_kmax;i++)
        {
          for (j=0;j<Size;j++)
           {
   #ifdef PRINT_OUT
              mexPrintf("index_vara[%d]=index_vara[%d]+y_size=",j+Size*i,j+Size*(i-1));
   #endif
             index_vara[j+Size*i]=index_vara[j+Size*(i-1)]+y_size;
   #ifdef PRINT_OUT
             mexPrintf("%d\n",index_vara[j+Size*i]);
  #endif
           }
       }
    }
  index_equa=(int*)mxMalloc(Size*sizeof(int));
  for(j=0;j<Size;j++)
    {
      SaveCode.read(reinterpret_cast<char *>(&index_equa[j]), sizeof(*index_equa));
    }
}



void SparseMatrix::Simple_Init(int it_, int y_kmin, int y_kmax, int Size, std::map<std::pair<std::pair<int, int> ,int>, int> &IM)
{
  int i, eq, var, lag;
  //double tmp_b=0.0;
  std::map<std::pair<std::pair<int, int> ,int>, int>::iterator it4;
  NonZeroElem* first;
  //mexPrintf("periods=%d, y_kmin=%d, y_kmax=%d, SizeInit=%d, IM.size()=%d\n",periods, y_kmin, y_kmax, Size, IM.size());
  pivot=(int*)mxMalloc(Size*sizeof(int));
  pivot_save=(int*)mxMalloc(Size*sizeof(int));
  pivotk=(int*)mxMalloc(Size*sizeof(int));
  pivotv=(double*)mxMalloc(Size*sizeof(double));
  pivotva=(double*)mxMalloc(Size*sizeof(double));
  b=(int*)mxMalloc(Size*sizeof(int));
  line_done=(bool*)mxMalloc(Size*sizeof(bool));
  //memset(line_done, 0, Size*sizeof(*line_done));

  mem_mngr.init_CHUNK_BLCK_SIZE(u_count);
  g_save_op=NULL;
  g_nop_all=0;
  i=Size*sizeof(NonZeroElem*);
  FNZE_R=(NonZeroElem**)mxMalloc(i);
  FNZE_C=(NonZeroElem**)mxMalloc(i);
  //memset(FNZE_R, 0, i);
  //memset(FNZE_C, 0, i);
  NonZeroElem** temp_NZE_R=(NonZeroElem**)mxMalloc(i);
  NonZeroElem** temp_NZE_C=(NonZeroElem**)mxMalloc(i);
  //memset(temp_NZE_R, 0, i);
  //memset(temp_NZE_C, 0, i);
  i=Size*sizeof(int);
  NbNZRow=(int*)mxMalloc(i);
  NbNZCol=(int*)mxMalloc(i);
  //memset(NbNZRow, 0, i);
  //memset(NbNZCol, 0, i);
  i=Size*sizeof(*b);
  //memset(b,0,i);
  it4=IM.begin();
  eq=-1;
  double tmp_b[Size];
  #pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
  for(i=0; i< Size;i++)
    {
      tmp_b[i]=0;//u[i];
      b[i]=0;
      line_done[i]=0;
      FNZE_C[i]=0;
      FNZE_R[i]=0;
      temp_NZE_C[i]=0;
      temp_NZE_R[i]=0;
      NbNZRow[i]=0;
      NbNZCol[i]=0;
    }
  int u_count1=Size;
  while (it4!=IM.end())
    {
      var=it4->first.first.second;
      /*if (eq!=it4->first.first.first)
        tmp_b=0;*/
      eq=it4->first.first.first;
      lag=it4->first.second;
      if (lag==0)   /*Build the index for sparse matrix containing the jacobian : u*/
        {
          //mexPrintf("Add eq=%d, var=%d, lag=%d at it_=%d u=%f\n",eq,var,lag, it_, u[u_count1]);
          //mexPrintf("    u_index=%d\n",/*it4->second+u_count_init*it_*/u_count1);
          NbNZRow[eq]++;
          NbNZCol[var]++;
#ifdef NEW_ALLOC
          first=mem_mngr.mxMalloc_NZE();
#else
          first=(NonZeroElem*)mxMalloc(sizeof(NonZeroElem));
#endif
          first->NZE_C_N=NULL;
          first->NZE_R_N=NULL;
          first->u_index=u_count1/*it4->second+u_count_init*it_*/;
          first->r_index=eq;
          first->c_index=var;
          first->lag_index=lag;
          //mexPrintf("  u[%d](%f)*y[%d](%f)=%f\n",u_count1, u[u_count1], index_vara[var]+it_*y_size, y[index_vara[var]+it_*y_size], u[u_count1]*y[index_vara[var]+it_*y_size]);
          tmp_b[eq] += u[u_count1]*y[index_vara[var]+it_*y_size];
          if (FNZE_R[eq]==NULL)
            {
              FNZE_R[eq]=first;
            }
          if (FNZE_C[var]==NULL)
            FNZE_C[var]=first;
          if (temp_NZE_R[eq]!=NULL)
            temp_NZE_R[eq]->NZE_R_N=first;
          if (temp_NZE_C[var]!=NULL)
            temp_NZE_C[var]->NZE_C_N=first;
          temp_NZE_R[eq]=first;
          temp_NZE_C[var]=first;
          u_count1++;
        }
      it4++;
    }
  #pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
  for(i=0;i<Size;i++)
    {
      b[i]=u_count1+i;
      u[b[i]]=-tmp_b[i];
    }
  //mexEvalString("Init");
  mxFree(temp_NZE_R);
  mxFree(temp_NZE_C);
  /*mexPrintf("end of Simple_Init\n");
  mexEvalString("drawnow;");*/
}


void SparseMatrix::Init(int periods, int y_kmin, int y_kmax, int Size, std::map<std::pair<std::pair<int, int> ,int>, int> &IM)
{
  int t,i, eq, var, lag, ti_y_kmin, ti_y_kmax;
  double tmp_b=0.0;
  std::map<std::pair<std::pair<int, int> ,int>, int>::iterator it4;
  NonZeroElem* first;
  //mexPrintf("periods=%d, y_kmin=%d, y_kmax=%d, SizeInit=%d, IM.size()=%d\n",periods, y_kmin, y_kmax, Size, IM.size());
#ifdef MEM_ALLOC_CHK
  mexPrintf("pivot=(int*)mxMalloc(%d*sizeof(int))\n",Size*periods);
#endif
  pivot=(int*)mxMalloc(Size*periods*sizeof(int));
  pivot_save=(int*)mxMalloc(Size*periods*sizeof(int));
#ifdef MEM_ALLOC_CHK
  mexPrintf("pivota=(int*)mxMalloc(%d*sizeof(int))\n",Size*periods);
#endif
  pivotk=(int*)mxMalloc(Size*periods*sizeof(int));
#ifdef MEM_ALLOC_CHK
  mexPrintf("pivotv=(double*)mxMalloc(%d*sizeof(double))\n",Size*periods);
#endif
  pivotv=(double*)mxMalloc(Size*periods*sizeof(double));
#ifdef MEM_ALLOC_CHK
  mexPrintf("pivotva=(double*)mxMalloc(%d*sizeof(double))\n",Size*periods);
#endif
  pivotva=(double*)mxMalloc(Size*periods*sizeof(double));
#ifdef MEM_ALLOC_CHK
  mexPrintf("b=(int*)mxMalloc(%d*sizeof(int))\n",Size*periods);
#endif
  b=(int*)mxMalloc(Size*periods*sizeof(int));
#ifdef MEM_ALLOC_CHK
  mexPrintf("line_done=(bool*)mxMalloc(%d*sizeof(bool))\n",Size*periods);
#endif
  line_done=(bool*)mxMalloc(Size*periods*sizeof(bool));
  //memset(line_done, 0, periods*Size*sizeof(*line_done));
  mem_mngr.init_CHUNK_BLCK_SIZE(u_count);
  g_save_op=NULL;
  g_nop_all=0;
#ifdef PRINT_OUT
  mexPrintf("sizeof(NonZeroElem)=%d sizeof(NonZeroElem*)=%d\n",sizeof(NonZeroElem),sizeof(NonZeroElem*));
#endif
  i=(periods+y_kmax+1)*Size*sizeof(NonZeroElem*);
#ifdef MEM_ALLOC_CHK
  mexPrintf("FNZE_R=(NonZeroElem**)mxMalloc(%d)\n",i);
#endif
  FNZE_R=(NonZeroElem**)mxMalloc(i);
#ifdef MEM_ALLOC_CHK
  mexPrintf("FNZE_C=(NonZeroElem**)mxMalloc(%d)\n",i);
#endif
  FNZE_C=(NonZeroElem**)mxMalloc(i);
  //memset(FNZE_R, 0, i);
  //memset(FNZE_C, 0, i);
#ifdef MEM_ALLOC_CHK
  mexPrintf("temp_NZE_R=(NonZeroElem**)(%d)\n",i);
#endif
  NonZeroElem** temp_NZE_R=(NonZeroElem**)mxMalloc(i);
#ifdef MEM_ALLOC_CHK
  mexPrintf("temp_NZE_R=(NonZeroElem**)(%d)\n",i);
#endif
  NonZeroElem** temp_NZE_C=(NonZeroElem**)mxMalloc(i);
  //memset(temp_NZE_R, 0, i);
  //memset(temp_NZE_C, 0, i);
  i=(periods+y_kmax+1)*Size*sizeof(int);
#ifdef MEM_ALLOC_CHK
  mexPrintf("NbNZRow=(int*)mxMalloc(%d)\n",i);
#endif
  NbNZRow=(int*)mxMalloc(i);
#ifdef MEM_ALLOC_CHK
  mexPrintf("NbNZCol=(int*)mxMalloc(%d)\n",i);
#endif
  NbNZCol=(int*)mxMalloc(i);
#ifdef MEM_ALLOC_CHK
  mexPrintf("ok\n");
#endif
  //memset(NbNZRow, 0, i);
  //memset(NbNZCol, 0, i);
  //i=periods*Size*sizeof(*b);
  //memset(b,0,i);

  #pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
  for(i=0; i< periods*Size;i++)
    {
      b[i]=0;
      line_done[i]=0;
    }
  #pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
  for(i=0; i< (periods+y_kmax+1)*Size;i++)
    {
      FNZE_C[i]=0;
      FNZE_R[i]=0;
      temp_NZE_C[i]=0;
      temp_NZE_R[i]=0;
      NbNZRow[i]=0;
      NbNZCol[i]=0;
    }

#ifdef PRINT_OUT
  mexPrintf("Now looping\n");
  mexEvalString("drawnow;");
#endif
  ///#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) ordered private(it4, ti_y_kmin, ti_y_kmax, eq, var, lag) schedule(dynamic)
  for (t=0;t<periods;t++)
    {
#ifdef PRINT_OUT
      mexPrintf("t=%d\n",t);
#endif
      ti_y_kmin=-min( t            , y_kmin);
      ti_y_kmax= min( periods-(t+1), y_kmax);
      it4=IM.begin();
      eq=-1;
      ///#pragma omp ordered
      while (it4!=IM.end())
        {
          var=it4->first.first.second;
          if (eq!=it4->first.first.first+Size*t)
            tmp_b=0;
          eq=it4->first.first.first+Size*t;
          lag=it4->first.second;
#ifdef PRINT_OUT
					mexPrintf("=) eq=%d var=%d lag=%d t=%d\n",eq,var, lag, t);
          mexPrintf("eq=%d, var=%d",eq,var);
          mexEvalString("drawnow;");
#endif
          if (var<(periods+y_kmax)*Size)
            {
              lag=it4->first.second;
#ifdef PRINT_OUT
              mexPrintf(", lag =%d, ti_y_kmin=%d, ti_y_kmax=%d ", lag, ti_y_kmin, ti_y_kmax);
#endif
              if (lag<=ti_y_kmax && lag>=ti_y_kmin)   /*Build the index for sparse matrix containing the jacobian : u*/
                {
                  var+=Size*t;
                  NbNZRow[eq]++;
                  NbNZCol[var]++;
#ifdef NEW_ALLOC
                  first=mem_mngr.mxMalloc_NZE();
#else
                  first=(NonZeroElem*)mxMalloc(sizeof(NonZeroElem));
#endif
								  //mexPrintf("=> eq=%d var=%d lag=%d u=%d\n",eq,var, lag, it4->second+u_count_init*t);
                  first->NZE_C_N=NULL;
                  first->NZE_R_N=NULL;
                  first->u_index=it4->second+u_count_init*t;
                  first->r_index=eq;
                  first->c_index=var;
                  first->lag_index=lag;
                  /*if(eq==0 && var==0)
                    mexPrintf("alloc FNZE_R[0]=%x\n",first);*/
                  if (FNZE_R[eq]==NULL)
                    FNZE_R[eq]=first;
                  if (FNZE_C[var]==NULL)
                    FNZE_C[var]=first;
                  if (temp_NZE_R[eq]!=NULL)
                    temp_NZE_R[eq]->NZE_R_N=first;
                  if (temp_NZE_C[var]!=NULL)
                    temp_NZE_C[var]->NZE_C_N=first;
                  temp_NZE_R[eq]=first;
                  temp_NZE_C[var]=first;
#ifdef PRINT_OUT
                  mexPrintf("=> ");
#endif
                }
              else       /*Build the additive terms ooutside the simulation periods related to the first lags and the last leads...*/
                {
                	if(lag<ti_y_kmin)
                	  {
#ifdef PRINT_OUT
                      mexPrintf("nn var=%d, Size=%d, t=%d, y_kmin=%d, y_kmax=%d\n", var, Size, t, y_kmin, y_kmax);
                      mexPrintf("   tmp_b+=u[%d]*y[index_var[%d]]\n", it4->second+u_count_init*t, var+Size*(y_kmin+t));
                      mexPrintf("   tmp_b+=u[%d](%f)*y[%d(%d)](f)\n", it4->second+u_count_init*t, u[it4->second+u_count_init*t], index_vara[var+Size*(y_kmin+t)],var+Size*(y_kmin+t)/*,y[index_vara[var+Size*(y_kmin+t)]]*/);
                      mexEvalString("drawnow;");
#endif
                      tmp_b+=u[it4->second+u_count_init*t]*y[index_vara[var+Size*(y_kmin+t)]];
                	  }
									else
									  {
#ifdef PRINT_OUT
									  	var -= Size;
                      mexPrintf("nn var=%d, Size=%d, t=%d, y_kmin=%d, y_kmax=%d\n", var, Size, t, y_kmin, y_kmax);
                      mexPrintf("   tmp_b+=u[%d]*y[index_var[%d]]\n", it4->second+u_count_init*t, var+Size*(y_kmin+t));
                      mexPrintf("   tmp_b+=u[%d](%f)*y[%d(%d)](f)\n", it4->second+u_count_init*t, u[it4->second+u_count_init*t], index_vara[var+Size*(y_kmin+t)],var+Size*(y_kmin+t)/*,y[index_vara[var+Size*(y_kmin+t)]]*/);
                      mexEvalString("drawnow;");
#endif
                      tmp_b+=u[it4->second+u_count_init*t]*y[index_vara[var+Size*(y_kmin+t)]];

									  }
                }
            }
          else           /* ...and store it in the u vector*/
            {
#ifdef PRINT_OUT
              mexPrintf("");
#endif
              b[eq]=it4->second+u_count_init*t;
              u[b[eq]]+=tmp_b;
              tmp_b = 0;
              //mexPrintf("u[%d]=%f corr=%f\n",b[eq],u[b[eq]],tmp_b);
#ifdef PRINT_OUT
              mexPrintf("=> u[b[%d]=%d]=%f\n", eq, b[eq], u[b[eq]]);
              mexEvalString("drawnow;");
#endif
            }
#ifdef PRINT_OUT
          mexPrintf(" u[%d] = %e\n",it4->second+u_count_init*t,double(u[it4->second+u_count_init*t]));
          mexEvalString("drawnow;");
#endif
          it4++;
        }
    }
#ifdef PRINT_OUT
  mexPrintf("end of Init\n");
  mexEvalString("drawnow;");
#endif
  mxFree(temp_NZE_R);
  mxFree(temp_NZE_C);
}

void SparseMatrix::ShortInit(int periods, int y_kmin, int y_kmax, int Size, std::map<std::pair<std::pair<int, int> ,int>, int> &IM)
{
  int t, eq, var, lag, ti_y_kmin, ti_y_kmax;
  double tmp_b=0.0;
  std::map<std::pair<std::pair<int, int> ,int>, int>::iterator it4;
  //#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) ordered private(it4, ti_y_kmin, ti_y_kmax, eq, var, lag, tmp_b) schedule(dynamic)
  for (t=0;t<periods;t++)
    {
#ifdef PRINT_OUT
      mexPrintf("t=%d\n",t);
#endif
      ti_y_kmin=-min( t            , y_kmin);
      ti_y_kmax= min( periods-(t+1), y_kmax);
      it4=IM.begin();
      eq=-1;
      while (it4!=IM.end())
        {
          var=it4->first.first.second;
          if (eq!=it4->first.first.first+Size*t)
            tmp_b=0;
          eq=it4->first.first.first+Size*t;
#ifdef PRINT_OUT
          mexPrintf("eq=%d, var=%d",eq,var);
#endif
          if (var<(periods+y_kmax)*Size)
            {
              lag=it4->first.second;
#ifdef PRINT_OUT
              mexPrintf(", lag =%d, ti_y_kmin=%d, ti_y_kmax=%d ", lag, ti_y_kmin, ti_y_kmax);
#endif
              if (lag<=ti_y_kmax && lag>=ti_y_kmin)
                {
                  var+=Size*t;
                }
              else
                {
#ifdef PRINT_OUT
                  mexPrintf("nn ");
                  mexPrintf("tmp_b+=u[%d]*y[index_var[%d]]\n",it4->second+u_count_init*t,var+Size*(y_kmin+t));
                  mexPrintf("tmp_b+=u[%d](%f)*y[%d(%d)](%f)",it4->second+u_count_init*t,u[it4->second+u_count_init*t], index_vara[var+Size*(y_kmin+t)],var+Size*(y_kmin+t),y[index_vara[var+Size*(y_kmin+t)]]);
#endif
                  tmp_b+=u[it4->second+u_count_init*t]*y[index_vara[var+Size*(y_kmin+t)]];
                }
            }
          else
            {
#ifdef PRINT_OUT
              mexPrintf("");
#endif
              b[eq]=it4->second+u_count_init*t;
              u[b[eq]]+=tmp_b;
              //mexPrintf("u[%d]=%f\n",b[eq],u[b[eq]]);
#ifdef PRINT_OUT
              mexPrintf("=> b[%d]=%f\n", eq, u[b[eq]]);
#endif
            }
#ifdef PRINT_OUT
          mexPrintf(" u[%d] = %e\n",it4->second+u_count_init*t,double(u[it4->second+u_count_init*t]));
#endif
          it4++;
        }
    }
}



int SparseMatrix::Get_u()
{
  if (!u_liste.empty())
    {
      int i=u_liste.back();
      u_liste.pop_back();
#ifdef PRINT_OUT
      mexPrintf("Get_u=%d\n",i);
#endif
      return i;
    }
  else
    {
      if (u_count<u_count_alloc)
        {
          int i=u_count;
          u_count++;
#ifdef PRINT_OUT
          mexPrintf("Get_u=%d\n",i);
#endif
          return i;
        }
      else
        {
          u_count_alloc+=5*u_count_alloc_save;
#ifdef MEM_ALLOC_CHK
          mexPrintf("u=(double*)mxRealloc(u,%d*sizeof(double))\n",u_count_alloc);
#endif
          u=(double*)mxRealloc(u,u_count_alloc*sizeof(double));
#ifdef MEM_ALLOC_CHK
          mexPrintf("ok\n");
#endif
          if (!u)
            {
              mexPrintf("Error in Get_u: memory exhausted (realloc(%d))\n",u_count_alloc*sizeof(double));
              mexEvalString("st=fclose('all');clear all;");
              mexErrMsgTxt("Exit from Dynare");
            }
          int i=u_count;
          u_count++;
          return i;
        }
    }
}

void SparseMatrix::Delete_u(int pos)
{
#ifdef PRINT_OUT
  mexPrintf("Delete_u=%d\n",pos);
#endif
  u_liste.push_back(pos);

}

void SparseMatrix::Clear_u()
{
  u_liste.clear();
}

void SparseMatrix::Print_u()
{
  for (unsigned int i=0;i<u_liste.size();i++)
    mexPrintf("%d ",u_liste[i]);
}

void SparseMatrix::End(int Size)
{
#ifdef NEW_ALLOC
  mem_mngr.Free_All();
#else
  for (int i=0;i<Size*periods;i++)
    {
      NonZeroElem *first=FNZE_R[i];
      while (!first)
        {
          NonZeroElem *firsta=first->NZE_R_N;
          mxFree(first);
          first=firsta;
        }
    }
#endif
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
SparseMatrix::compare( int *save_op, int *save_opa, int *save_opaa, int beg_t, int periods, long int nop4,  int Size
#ifdef PROFILER
, long int *ndiv, long int *nsub
#endif
)
{
  long int i,j,nop=nop4/2, t, index_d, k;
  double r=0.0;
  bool OK=true;
  t_save_op_s *save_op_s, *save_opa_s, *save_opaa_s;
  int *diff1, *diff2;
  diff1=(int*)mxMalloc(nop*sizeof(int));
  diff2=(int*)mxMalloc(nop*sizeof(int));
  int max_save_ops_first=-1;
  j=k=i=0;
  while (i<nop4 && OK)
    {
      save_op_s=(t_save_op_s*)&(save_op[i]);
      save_opa_s=(t_save_op_s*)&(save_opa[i]);
      save_opaa_s=(t_save_op_s*)&(save_opaa[i]);
      diff1[j]=save_op_s->first-save_opa_s->first;
      if(max_save_ops_first<save_op_s->first+diff1[j]*(periods-beg_t))
        {
          max_save_ops_first=save_op_s->first+diff1[j]*(periods-beg_t);
        }
      switch (save_op_s->operat)
        {
          case IFLD:
          case IFDIV:
            OK=(save_op_s->operat==save_opa_s->operat && save_opa_s->operat==save_opaa_s->operat
                && diff1[j]==(save_opa_s->first-save_opaa_s->first));
            i+=2;
            break;
          case IFLESS:
          case IFSUB:
            diff2[j]=save_op_s->second-save_opa_s->second;
            OK=(save_op_s->operat==save_opa_s->operat && save_opa_s->operat==save_opaa_s->operat
                && diff1[j]==(save_opa_s->first-save_opaa_s->first)
                && diff2[j]==(save_opa_s->second-save_opaa_s->second));
            i+=3;
            break;
          default:
            mexPrintf("unknown operator = %d ",save_op_s->operat);
            mexEvalString("st=fclose('all');clear all;");
            filename+=" stopped";
            mexErrMsgTxt(filename.c_str());
            break;
        }
      j++;
    }
  // the same pivot for all remaining periods
  if (OK)
    //#pragma omp parallel for  num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) ordered private(j) schedule(dynamic)
    for (i=beg_t;i<periods;i++)
      {
        for (j=0;j<Size;j++)
          {
            ///#pragma omp ordered
            pivot[i*Size+j]=pivot[(i-1)*Size+j]+Size;
          }
      }
  if (OK)
    {
      if (max_save_ops_first>=u_count_alloc)
        {
          u_count_alloc+=5*u_count_alloc_save;
          u=(double*)mxRealloc(u,u_count_alloc*sizeof(double));
          if (!u)
            {
              mexPrintf("Error in Get_u: memory exhausted (realloc(%d))\n",u_count_alloc*sizeof(double));
              mexEvalString("st=fclose('all');clear all;");
              mexErrMsgTxt("Exit from Dynare");
            }
          }
      for (t=1;t<periods-beg_t-y_kmax/*max(y_kmax,y_kmin)*/;t++)
        {
          i=j=0;
          while (i<nop4)
            {
              save_op_s=(t_save_op_s*)(&(save_op[i]));
              index_d=save_op_s->first+t*diff1[j];
              if (index_d>u_count_alloc)
                {
                  u_count_alloc+=2*u_count_alloc_save;
                  u=(double*)mxRealloc(u,u_count_alloc*sizeof(double));
                  if (!u)
                    {
                      mexPrintf("Error in Get_u: memory exhausted (realloc(%d))\n",u_count_alloc*sizeof(double));
                      mexEvalString("st=fclose('all');clear all;");
                      mexErrMsgTxt("Exit from Dynare");
                    }
                }
              switch (save_op_s->operat)
                {
                  case IFLD  :
                    r=u[index_d];
                    i+=2;
                    break;
                  case IFDIV :
                    u[index_d]/=r;
                    i+=2;
                    break;
                  case IFSUB :
                    u[index_d]-=u[save_op_s->second+t*diff2[j]]*r;
                    i+=3;
                    break;
                  case IFLESS:
                    u[index_d]=-u[save_op_s->second+t*diff2[j]]*r;
                    i+=3;
                    break;
                }
              j++;
            }
        }
      int t1=max(1,periods-beg_t-y_kmax);
      //#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) ordered private(t, i,j, save_op_s, index_d, r) schedule(dynamic)
      for (t=t1;t<periods-beg_t;t++)
        {
          i=j=0;
          //#pragma omp ordered
          while (i<nop4)
            {
              save_op_s=(t_save_op_s*)(&(save_op[i]));
              if (save_op_s->lag<((periods-beg_t)-t))
                {
                  index_d=save_op_s->first+t*diff1[j];
                  if (index_d>u_count_alloc)
                    {
                      u_count_alloc+=2*u_count_alloc_save;
                      u=(double*)mxRealloc(u,u_count_alloc*sizeof(double));
                      if (!u)
                        {
                          mexPrintf("Error in Get_u: memory exhausted (realloc(%d))\n",u_count_alloc*sizeof(double));
                          mexEvalString("st=fclose('all');clear all;");
                          mexErrMsgTxt("Exit from Dynare");
                        }
                    }
                  switch (save_op_s->operat)
                    {
                      case IFLD  :
                        r=u[index_d];
                        i+=2;
                        break;
                      case IFDIV :
                        u[index_d]/=r;
                        i+=2;
                        break;
                      case IFSUB :
                        u[index_d]-=u[save_op_s->second+t*diff2[j]]*r;
                        i+=3;
                        break;
                      case IFLESS:
                        u[index_d]=-u[save_op_s->second+t*diff2[j]]*r;
                        i+=3;
                        break;
                    }
                }
              else
                {
                  switch (save_op_s->operat)
                    {
                      case IFLD  :
                      case IFDIV :
                        i+=2;
                        break;
                      case IFSUB :
                      case IFLESS :
                        i+=3;
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

void
SparseMatrix::run_u_period1(int periods)
{
  double r=0;
  int index_d, t;
  for (t=0;t<periods;t++)
    {
      for (long int i=0;i<g_nop_all;)
        {
          index_d=g_save_op[i+1]+t*g_save_op[i+2];
          if (index_d>=u_count_alloc)
            {
              u_count_alloc+=5*u_count_alloc_save;
#ifdef MEM_ALLOC_CHK
              mexPrintf("u=(double*)mxRealloc(u,u_count_alloc*sizeof(double))\n",u_count_alloc);
#endif
              u=(double*)mxRealloc(u,u_count_alloc*sizeof(double));
#ifdef MEM_ALLOC_CHK
              mexPrintf("ok\n");
#endif
              if (!u)
                {
                  mexPrintf("Error in Get_u: memory exhausted (realloc(%d))\n",u_count_alloc*sizeof(double));
                  mexEvalString("st=fclose('all');clear all;");
                  mexErrMsgTxt("Exit from Dynare");
                }
            }
          switch (g_save_op[i])
            {
              case IFLD  :
                r=u[index_d];
#ifdef PRINT_u
                mexPrintf("FLD u[%d] (%f)\n",index_d,u[index_d]);
#endif
                break;
              case IFDIV :
                u[index_d]/=r;
#ifdef PRINT_u
                mexPrintf("FDIV u[%d](%f)/=r(%f)=(%f)\n",index_d,u[index_d],r,u[index_d]);
#endif
                break;
              case IFSUB :
                u[index_d]-=u[g_save_op[i+3]+t*g_save_op[i+4]]*r;
#ifdef PRINT_u
                mexPrintf("FSUB u[%d]-=u[%d](%f)*r(%f)=(%f) index1=%d index2=%d\n",index_d,g_save_op[i+3]+t*g_save_op[i+4],u[g_save_op[i+3]+t*g_save_op[i+4]],r,u[index_d],g_save_op[i+1],g_save_op[i+2] );
#endif
                break;
              case IFLESS:
                u[index_d]=-u[g_save_op[i+3]+t*g_save_op[i+4]]*r;
#ifdef PRINT_u
                mexPrintf("FLESS u[%d]=-u[%d](%f)*r(%f)=(%f) index1=%d index2=%d\n",index_d,g_save_op[i+3]+t*g_save_op[i+4],u[g_save_op[i+3]+t*g_save_op[i+4]],r,u[index_d],g_save_op[i+1],g_save_op[i+2] );
#endif
                break;
            }
          i+=5;
        }
    }
}




void
SparseMatrix::run_it(int nop_all,int *op_all)
{
  double r=0;
  int index_d;
  for (long int i=0;i<nop_all;)
    {
      index_d=op_all[i+1];
      if (index_d>=u_count_alloc)
        {
          u_count_alloc+=5*u_count_alloc_save;
#ifdef MEM_ALLOC_CHK
          mexPrintf("u=(double*)mxRealloc(u,u_count_alloc*sizeof(double))\n",u_count_alloc);
#endif
          u=(double*)mxRealloc(u,u_count_alloc*sizeof(double));
#ifdef MEM_ALLOC_CHK
          mexPrintf("ok\n");
#endif
          if (!u)
            {
              mexPrintf("Error in Get_u: memory exhausted (realloc(%d))\n",u_count_alloc*sizeof(double));
              mexEvalString("st=fclose('all');clear all;");
              mexErrMsgTxt("Exit from Dynare");
            }
        }
      switch (op_all[i])
        {
          case IFLD  :
            r=u[index_d];
#ifdef PRINT_u
            mexPrintf("FLD u[%d] (%f)\n",index_d,u[index_d]);
#endif
            i+=2;
            break;
          case IFDIV :
            u[index_d]/=r;
#ifdef PRINT_u
            mexPrintf("FDIV u[%d](%f)/=r(%f)=(%f)\n",index_d,u[index_d],r,u[index_d]);
#endif
            i+=2;
            break;
          case IFSUB :
            u[index_d]-=u[op_all[i+2]]*r;
#ifdef PRINT_u
            mexPrintf("FSUB u[%d]-=u[%d](%f)*r(%f)=(%f)\n",index_d,op_all[i+2],u[op_all[i+2]],r,u[index_d]);
#endif
            i+=3;
            break;
          case IFLESS:
            u[index_d]=-u[op_all[i+2]]*r;
#ifdef PRINT_u
            mexPrintf("FLESS u[%d]=-u[%d](%f)*r(%f)=(%f)\n",index_d,op_all[i+2],u[op_all[i+2]],r,u[index_d]);
#endif
            i+=3;
            break;
        }
    }
}


void
SparseMatrix::run_triangular(int nop_all,int *op_all)
{
  int j=0;
  //mexPrintf("begining of run_triangular nop_all=%d\n",nop_all);
  if (mem_mngr.swp_f)
    {
      bool OK=true;
      int* save_op;
      long int nop;
      while (OK)
        {
          mexPrintf("reading blck%d\n",j++);
          OK=mem_mngr.read_swp_f(&save_op,&nop);
          if (OK)
            {
              run_it(nop,save_op);
              mxFree(save_op);
            }
        }
    }
  run_it(nop_all,op_all);
}

int
SparseMatrix::complete(int beg_t, int Size, int periods, int *b)
{
  long int i, j, k, nop, nopa, nop1, cal_y, nb_var, pos, t, ti, max_var, min_var;
  NonZeroElem *first;
  int *save_code;
  int *diff;
  double yy=0.0, err;

  int size_of_save_code=(1+y_kmax)*Size*(Size+1+4)/2*4;
#ifdef MEM_ALLOC_CHK
  mexPrintf("save_code=(int*)mxMalloc(%d*sizeof(int))\n",size_of_save_code);
#endif
  save_code=(int*)mxMalloc(size_of_save_code*sizeof(int));
  int size_of_diff=(1+y_kmax)*Size*(Size+1+4);
#ifdef MEM_ALLOC_CHK
  mexPrintf("diff=(int*)mxMalloc(%d*sizeof(int))\n",size_of_diff);
#endif
  diff=(int*)mxMalloc(size_of_diff*sizeof(int));
#ifdef MEM_ALLOC_CHK
  mexPrintf("ok\n");
#endif
  cal_y=y_size*y_kmin;

  i=(beg_t+1)*Size-1;
  nop=0;
  for (j=i;j>i-Size;j--)
    {
      pos=pivot[j];
      nb_var=At_Row(pos,&first);
      first=first->NZE_R_N;
      nb_var--;
      save_code[nop]=IFLDZ;
      save_code[nop+1]=0;
      save_code[nop+2]=0;
      save_code[nop+3]=0;
      if ((nop+3)>=size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n",nop+2,size_of_save_code);
      nop+=4;
      for (k=0;k<nb_var;k++)
        {
          save_code[nop]=IFMUL;
          save_code[nop+1]=index_vara[first->c_index]+cal_y;
          save_code[nop+2]=first->u_index;
          save_code[nop+3]=first->lag_index;
          if ((nop+3)>=size_of_save_code)
            mexPrintf("out of save_code[%d] (bound=%d)\n",nop+2,size_of_save_code);
          nop+=4;
          first=first->NZE_R_N;
        }
      //yy=-(yy+u[b[pos]]);
      //mexPrintf("|u[%d]|\n",b[pos]);
      save_code[nop]=IFADD;
      save_code[nop+1]=b[pos];
      save_code[nop+2]=0;
      save_code[nop+3]=0;
      if ((nop+3)>=size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n",nop+2,size_of_save_code);
      nop+=4;
      save_code[nop]=IFSTP;
      save_code[nop+1]=index_vara[j]+y_size*y_kmin;
      save_code[nop+2]=0;
      save_code[nop+3]=0;
      if ((nop+2)>=size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n",nop+2,size_of_save_code);
      nop+=4;
    }
  i=beg_t*Size-1;
  nop1=nopa=0;
  for (j=i;j>i-Size;j--)
    {
      pos=pivot[j];
      nb_var=At_Row(pos,&first);
      first=first->NZE_R_N;
      nb_var--;
      diff[nopa]=0;
      diff[nopa+1]=0;
      nopa+=2;
      nop1+=4;
      for (k=0;k<nb_var;k++)
        {
          diff[nopa]=save_code[nop1+1]-(index_vara[first->c_index]+cal_y);
          diff[nopa+1]=save_code[nop1+2]-(first->u_index);
          if ((nop1+2)>=size_of_save_code)
            mexPrintf("out of save_code[%d] (bound=%d)\n",nop1+2,size_of_save_code);
          if ((nopa+1)>=size_of_diff)
            mexPrintf("out of diff[%d] (bound=%d)\n",nopa+2,size_of_diff);
          nopa+=2;
          nop1+=4;
          first=first->NZE_R_N;
        }
      diff[nopa]=save_code[nop1+1]-(b[pos]);
      diff[nopa+1]=0;
      if ((nop1+3)>=size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n",nop1+2,size_of_save_code);
      if ((nopa+1)>=size_of_diff)
        mexPrintf("out of diff[%d] (bound=%d)\n",nopa+2,size_of_diff);
      nopa+=2;
      nop1+=4;
      diff[nopa]=save_code[nop1+1]-(index_vara[j]+y_size*y_kmin);
      diff[nopa+1]=0;
      if ((nop1+4)>=size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n",nop1+2,size_of_save_code);
      if ((nopa+1)>=size_of_diff)
        mexPrintf("out of diff[%d] (bound=%d)\n",nopa+2,size_of_diff);
      nopa+=2;
      nop1+=4;
    }
  max_var=(periods+y_kmin)*y_size;
  min_var=y_kmin*y_size;
  int k1=0;
  for (t=periods+y_kmin-1;t>=beg_t+y_kmin;t--)
    {
      j=0;
      ti=t-y_kmin-beg_t;
      for (i=0;i<nop;i+=4)
        {
          switch (save_code[i])
            {
              case IFLDZ :
                yy=0;
                break;
              case IFMUL :
                k=save_code[i+1]+ti*diff[j];
                if (k<max_var && k>min_var)
                  {
                    yy+=y[k]*u[save_code[i+2]+ti*diff[j+1]];
#ifdef PRINT_OUT_y1
                    mexPrintf("y[%d](%f)*u[%d](%f)+",k, double(y[k]), save_code[i+2]+ti*diff[j+1], double(u[save_code[i+2]+ti*diff[j+1]]));
#endif
                  }
                break;
              case IFADD :
                yy=-(yy+u[save_code[i+1]+ti*diff[j]]);
#ifdef PRINT_OUT_y1
                mexPrintf("|u[%d](%f)|",save_code[i+1]+ti*diff[j],double(u[save_code[i+1]+ti*diff[j]]));
#endif
                break;
              case IFSTP :
                k=save_code[i+1]+ti*diff[j];
                k1=k;
                err = yy - y[k];
                y[k] += slowc*(err);
#ifdef PRINT_OUT_y1
                mexPrintf("=y[%d]=%f  diff[%d]=%d save_code[%d]=%d ti=%d\n",save_code[i+1]+ti*diff[j],y[k],j,diff[j],i+1,save_code[i+1],ti);
#endif
                break;
            }
          j+=2;
        }
    }
  mxFree(save_code);
  mxFree(diff);
  return(beg_t);
}

void
SparseMatrix::close_swp_file()
{
  mem_mngr.close_swp_f();
}










double
SparseMatrix::bksub( int tbreak, int last_period, int Size, double slowc_l
#ifdef PROFILER
, /*NonZeroElem *first,*/ long int *nmul
#endif
)
{
  NonZeroElem *first;
  int i, j, k;
  double yy;
  res1 = res2 = max_res = 0;
  for (i=0;i<y_size*(periods+y_kmin);i++)
    y[i]=ya[i];
  if (symbolic && tbreak)
    last_period=complete(tbreak, Size, periods, b);
  else
    last_period=periods;
  for (int t=last_period+y_kmin-1;t>=y_kmin;t--)
    {
      int ti=(t-y_kmin)*Size;
#ifdef PRINT_OUT
      mexPrintf("t=%d ti=%d\n",t,ti);
#endif
      int cal=y_kmin*Size;
      int cal_y=y_size*y_kmin;
      for (i=ti-1;i>=ti-Size;i--)
        {
          j=i+cal;
#ifdef PRINT_OUT_y
          mexPrintf("t=%d, ti=%d j=%d i+Size=%d\n",t,ti,j,i+Size);
#endif
          int pos=pivot[/*j*/i+Size];
#ifdef PRINT_OUT_y
          mexPrintf("i-ti+Size=%d pos=%d j=%d\n",i-ti+Size,pos,j);
#endif
          int nb_var=At_Row(pos,&first);
          first=first->NZE_R_N;
          nb_var--;
          int eq=index_vara[j]+y_size;
#ifdef PRINT_OUT_y1
          mexPrintf("y[index_vara[%d]=%d]=",j,index_vara[j]+y_size);
#endif
          yy=0;
          for (k=0;k<nb_var;k++)
            {
              yy+=y[index_vara[first->c_index]+cal_y]*u[first->u_index];
#ifdef PROFILER
              (*nmul)++;
#endif
              first=first->NZE_R_N;
            }
#ifdef PRINT_OUT_y1
          mexPrintf("|u[%d](%f)|",b[pos],double(u[b[pos]]));
#endif
          yy=-(yy+y[eq]+u[b[pos]]);
          direction[eq]=yy;
          //mexPrintf("direction[%d] = %f\n",eq,yy);
          y[eq] += slowc_l*yy;
#ifdef PRINT_OUT_y1
          mexPrintf("=%f (%f)\n",double(yy),double(y[eq]));
#endif
        }
    }
  return res1;
}


double
SparseMatrix::simple_bksub(int it_, int Size, double slowc_l)
{
  int i,k;
  double yy;
  NonZeroElem *first;
  res1 = res2 = max_res = 0;
  //mexPrintf("simple_bksub\n");
  for (i=0;i<y_size;i++)
    y[i+it_*y_size]=ya[i+it_*y_size];
  //mexPrintf("setp 1\n");
  //mexPrintf("Size=%d\n",Size);
  for(i=Size-1;i>=0;i--)
    {
      //mexPrintf("i=%d\n",i);
      int pos=pivot[i];
      //mexPrintf("pos=%d\n",pos);
      int nb_var=At_Row(pos,&first);
      //mexPrintf("nb_var=%d\n",nb_var);
      first=first->NZE_R_N;
      nb_var--;
      //mexPrintf("i=%d\n",i);
      int eq=index_vara[i];
      //mexPrintf("eq=%d\n",eq);
      yy = 0;
      for (k=0;k<nb_var;k++)
        {
          //mexPrintf("y[index_vara[%d]=%d]=%f\n",first->c_index, index_vara[first->c_index], y[index_vara[first->c_index]+it_*y_size]);
          //mexPrintf("u[first->u_index=%d]=%f\n",first->u_index, u[first->u_index]);
          yy+=y[index_vara[first->c_index]+it_*y_size]*u[first->u_index];
          first=first->NZE_R_N;
        }
      yy=-(yy+y[eq+it_*y_size]+u[b[pos]]);
      //mexPrintf("yy=%f\n",yy);
      direction[eq+it_*y_size]=yy;
      y[eq+it_*y_size] += slowc_l*yy;
    }
  return res1;
}


int
SparseMatrix::simulate_NG(int blck, int y_size, int it_, int y_kmin, int y_kmax, int Size, bool print_it, bool cvg, int &iter)
{
  int i, j, k;
  int pivj=0, pivk=0;
  double piv_abs, first_elem;
  NonZeroElem *first, *firsta, *first_sub, *first_piv, *first_suba;
#ifdef MARKOVITZ
  double piv_v[Size];
  int pivj_v[Size], pivk_v[Size], NR[Size], l, N_max;
  bool one;
#endif
#ifdef PROFILER
  //long int ndiv=0, nsub=0, ncomp=0, nmul=0;
  //double tinsert=0, tdelete=0, tpivot=0, tbigloop=0;
  //clock_t td1;
  int nbpivot=0, nbpivot_it=0;
  //int nbdiv=0, nbless=0, nbRealloc=0, insert=0;
#endif
  /*mexPrintf("begining\n");
  mexEvalString("drawnow;");*/
  if (cvg)
    return(0);
  /*mexPrintf("begining after cvg Size=%d\n", Size);
  mexEvalString("drawnow;");*/
  Simple_Init(it_, y_kmin, y_kmax, Size, IM_i);
  /*mexPrintf("begining after Simple_Init\n");
  mexEvalString("drawnow;");*/
  if (isnan(res1) || isinf(res1))
    {
      if (slowc_save<1e-8)
        {
          mexPrintf("slowc_save=%g\n", slowc_save);
          mexPrintf("Dynare cannot improve the simulation in block %d at time %d\n", blck, it_);
          mexEvalString("drawnow;");
          mexEvalString("st=fclose('all');clear all;");
          filename+=" stopped";
          mexErrMsgTxt(filename.c_str());
        }
      slowc_save/=2;
      mexPrintf("Error: Simulation diverging, trying to correct it using slowc=%f\n",slowc_save);
      for (i=0;i<y_size;i++)
        y[i+it_*y_size]=ya[i+it_*y_size]+slowc_save*direction[i+it_*y_size];
      iter--;
      return(0);
    }
  //mexPrintf("begining after nan\n");
  //mexEvalString("drawnow;");
  //mexPrintf("before the main loop y_size=%d\n",y_size);
  for (i=0;i<Size;i++)
    {
      //mexPrintf("i=%d\n",i);
      //mexEvalString("drawnow;");
      /*finding the max-pivot*/
      double piv=piv_abs=0;
      //int nb_eq=At_Col(i, 0, &first);
      int nb_eq=At_Col(i, &first);
      //mexPrintf("nb_eq=%d\n",nb_eq);
      //mexEvalString("drawnow;");
#ifdef MARKOVITZ
      l=0; N_max=0;
      one=false;
      piv_abs=0;
#endif
      for (j=0;j<nb_eq/*Size*/;j++)
        {
          //mexPrintf("j=%d \n",j);
          //mexEvalString("drawnow;");
          //mexPrintf("first->r_index=%d \n",first->r_index);
          //mexPrintf("line_done[%d]=%d \n",first->r_index,line_done[first->r_index]);
          //mexEvalString("drawnow;");
          if (!line_done[first->r_index])
            {
              //mexPrintf("first->u_index=%d \n",first->u_index);
              //mexEvalString("drawnow;");
              k=first->u_index;
              int jj=first->r_index;
              int NRow_jj=NRow(jj);
#ifdef PROFILER
              nbpivot++;
              nbpivot_it++;
#endif

#ifdef MARKOVITZ
              piv_v[l]=u[k];
              //mexPrintf("piv_v[%d]=%f\n",l, piv_v[l]);
              //mexEvalString("drawnow;");
              double piv_fabs=fabs(u[k]);
              pivj_v[l]=jj;
              pivk_v[l]=k;
              NR[l]=NRow_jj;
              if (NRow_jj==1 && !one)
                {
                  one=true;
                  piv_abs=piv_fabs;
                  N_max=NRow_jj;
                }
              if (!one)
                {
                  if (piv_fabs>piv_abs)
                    piv_abs=piv_fabs;
                  if (NRow_jj>N_max)
                    N_max=NRow_jj;
                }
              else
                {
                  if (NRow_jj==1)
                    {
                      if (piv_fabs>piv_abs)
                        piv_abs=piv_fabs;
                      if (NRow_jj>N_max)
                        N_max=NRow_jj;
                    }
                }
              l++;
#else
              if (piv_abs<fabs(u[k])||NRow_jj==1)
                {
                  piv=u[k];
                  piv_abs=fabs(piv);
                  pivj=jj;   //Line number
                  pivk=k;   //position in u
                  if (NRow_jj==1)
                  break;
                }
#endif
            }
          first=first->NZE_C_N;
        }
      /*mexPrintf("pivot found piv_v[0]=%f pivk_v[0]=%d l=%d one=%d\n", piv_v[0], pivk_v[0], l, one);
      mexEvalString("drawnow;");*/
#ifdef MARKOVITZ
      double markovitz=0, markovitz_max=-9e70;
      if (!one)
        {
          for (j=0;j<l;j++)
            {
              markovitz=exp(log(fabs(piv_v[j])/piv_abs)-markowitz_c*log(double(NR[j])/double(N_max)));
              if (markovitz>markovitz_max)
                {
                   piv=piv_v[j];
                   pivj=pivj_v[j];   //Line number
                   pivk=pivk_v[j];   //positi
                   markovitz_max=markovitz;
                }
            }
        }
      else
        {
          for (j=0;j<l;j++)
            {
              markovitz=exp(log(fabs(piv_v[j])/piv_abs)-markowitz_c*log(NR[j]/N_max));
              //mexPrintf("j=%d markovitz=%f\n",j, markovitz);
              if (markovitz>markovitz_max && NR[j]==1)
                {
                  piv=piv_v[j];
                  pivj=pivj_v[j];   //Line number
                  pivk=pivk_v[j];   //positi
                  markovitz_max=markovitz;
                  //mexPrintf("stored\n");
                }
            }
        }
#endif
      //mexPrintf("OK i=%d\n", i);
      pivot[i]=pivj;
      //mexPrintf("pivot[%d]=%d\n",i,pivot[i]);
      pivotk[i]=pivk;
      pivotv[i]=piv;
      line_done[pivj]=true;
      if (piv_abs<eps)
        {
          mexPrintf("Error: singular system in Simulate_NG\n");
          mexEvalString("drawnow;");
          mexEvalString("st=fclose('all');clear all;");
          filename+=" stopped";
          mexErrMsgTxt(filename.c_str());
        }
      //mexPrintf("piv=%f\n",piv);
      /*divide all the non zeros elements of the line pivj by the max_pivot*/
      int nb_var=At_Row(pivj,&first);
      for (j=0;j<nb_var;j++)
        {
          u[first->u_index]/=piv;
          first=first->NZE_R_N;
        }
      u[b[pivj]]/=piv;
      /*substract the elements on the non treated lines*/
      nb_eq=At_Col(i,&first);
      NonZeroElem *first_piva;
      int nb_var_piva=At_Row(pivj,&first_piva);
      //#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
      for (j=0;j<Size /*and first*/;j++)
        {
        	if(first)
        	{
          int row=first->r_index;
          //mexPrintf("j=%d row=%d line_done[row]=%d\n",j, row, line_done[row]);
          if (!line_done[row])
            {
              first_elem=u[first->u_index];
              int nb_var_piv=nb_var_piva;
              first_piv=first_piva;
              int nb_var_sub=At_Row(row,&first_sub);
              //mexPrintf("nb_var_sub=%d nb_var_piv=%d\n",nb_var_sub,nb_var_piv);
              int l_sub=0, l_piv=0;
              int sub_c_index=first_sub->c_index, piv_c_index=first_piv->c_index;
              //int tmp_lag=first_sub->lag_index;
              while (l_sub<nb_var_sub || l_piv<nb_var_piv)
                {
                  if (l_sub<nb_var_sub && (sub_c_index<piv_c_index || l_piv>=nb_var_piv))
                    {
                      first_sub=first_sub->NZE_R_N;
                      if (first_sub)
                        sub_c_index=first_sub->c_index;
                      else
                        sub_c_index=Size/* *periods*/;
                      l_sub++;
                    }
                  else if (sub_c_index>piv_c_index || l_sub>=nb_var_sub)
                    {
                      int tmp_u_count=Get_u();
                      //int lag=first_piv->c_index/Size-row/Size;
                      Insert(row,first_piv->c_index,tmp_u_count,/*lag*/0);
                      u[tmp_u_count]=-u[first_piv->u_index]*first_elem;
                      first_piv=first_piv->NZE_R_N;
                      if (first_piv)
                        piv_c_index=first_piv->c_index;
                      else
                        piv_c_index=Size/* *periods*/;
                      l_piv++;
                    }
                  else /*first_sub->c_index==first_piv->c_index*/
                    {
                      if (i==sub_c_index)
                        {
                          firsta=first;
                          first_suba=first_sub->NZE_R_N;
                          Delete(first_sub->r_index,first_sub->c_index, Size);
                          first=firsta->NZE_C_N;
                          first_sub=first_suba;
                          if (first_sub)
                            sub_c_index=first_sub->c_index;
                          else
                            sub_c_index=Size/* *periods*/;
                          l_sub++;
                          first_piv=first_piv->NZE_R_N;
                          if (first_piv)
                            piv_c_index=first_piv->c_index;
                          else
                            piv_c_index=Size;
                          l_piv++;
                        }
                      else
                        {
                          u[first_sub->u_index]-=u[first_piv->u_index]*first_elem;
                          first_sub=first_sub->NZE_R_N;
                          if (first_sub)
                            sub_c_index=first_sub->c_index;
                          else
                            sub_c_index=Size /* *periods*/;
                          l_sub++;
                          first_piv=first_piv->NZE_R_N;
                          if (first_piv)
                            piv_c_index=first_piv->c_index;
                          else
                            piv_c_index=Size/* *periods*/;
                          l_piv++;
                        }
                    }
                }
              u[b[row]]-=u[b[pivj]]*first_elem;
            }
          else
            first=first->NZE_C_N;
            /*first=first->NZE_R_N;*/
        }
        }
    }
  /*mexPrintf("before bcksub\n");
  mexEvalString("drawnow;");*/
  double slowc_lbx=slowc, res1bx;
  //mexPrintf("before bksub it_=%d\n",it_);
  for (i=0;i<y_size;i++)
    ya[i+it_*y_size]=y[i+it_*y_size];
  slowc_save=slowc;
  //res1bx=bksub( NULL, 1, y_size, slowc_lbx);
  res1bx=simple_bksub(it_,Size,slowc_lbx);
  //mexPrintf("End of simulate_NG\n");
  End(Size);
  return(0);
}




void
SparseMatrix::CheckIt(int y_size, int y_kmin, int y_kmax, int Size, int periods, int iter)
{
	const double epsilon=1e-7;
	fstream SaveResult;
  ostringstream out;
  out << "Result" << iter;
	SaveResult.open(out.str().c_str(), std::ios::in );
  if (!SaveResult.is_open())
    {
      mexPrintf("Error : Can't open file \"%s\" for reading\n", "Result");
      mexEvalString("st=fclose('all');clear all;");
      mexErrMsgTxt("Exit from Dynare");
    }
	mexPrintf("Reading Result...");
	int row, col;
	SaveResult >> row;
	mexPrintf("row=%d\n",row);
	SaveResult >> col;
	mexPrintf("col=%d\n",col);
	//double G1a[row][col];
	double G1a;
	mexPrintf("Allocated\n");
	NonZeroElem *first;
	for(int j=0; j< col; j++)
	  {
	  	mexPrintf("j=%d ",j);
      int nb_equ=At_Col(j,&first);
      mexPrintf("nb_equ=%d\n",nb_equ);
      int line;
      if(first)
      	line = first->r_index;
			else
			  line = -9999999;
	    for(int i=0; i< row; i++)
	      {
	        SaveResult >> G1a;
	        if(line == i)
            {
            	if(abs(u[first->u_index]/G1a-1)>epsilon)
            	  mexPrintf("Problem at r=%d c=%d u[first->u_index]=%5.14f G1a[i][j]=%5.14f %f\n",i,j,u[first->u_index],G1a, u[first->u_index]/G1a-1);
				  	  first=first->NZE_C_N;
					  	if(first)
						  	 line = first->r_index;
							else
  							 line = -9999999;
						}
		  		else
			  		{
				  	  if(G1a!=0.0)
  				  	  mexPrintf("Problem at r=%d c=%d G1a[i][j]=%f\n",i,j,G1a);
            }
	      }
	  }
	mexPrintf("G1a red done\n");
	SaveResult >> row;
	mexPrintf("row(2)=%d\n",row);
	double B[row];
	for(int i=0; i< row; i++)
	  SaveResult >> B[i];
	SaveResult.close();
  mexPrintf("done\n");
  mexPrintf("Comparing...");
  /*NonZeroElem *first;
  for(int i=0;i<row;i++)
    {
    	mexPrintf("i=%d ",i);
      int nb_var=At_Row(i,&first);
      mexPrintf("nb_var=%d\n",nb_var);
      int column;
      if(first)
      	column = first->c_index;
			else
			  column = -9999999;
      for(int j=0;j<col;j++)
        {
          if(column == j)
            {
            	if(abs(u[first->u_index]-G1a[i][j])>epsilon)
            	  mexPrintf("Problem at r=%d c=%d u[first->u_index]=%f G1a[i][j]=%f\n",i,j,u[first->u_index],G1a[i][j]);
				  	  first=first->NZE_R_N;
					  	if(first)
						  	 column = first->c_index;
							else
  							 column = -9999999;
						}
		  		else
			  		{
				  	  if(G1a[i][j]!=0)
  				  	  mexPrintf("Problem at r=%d c=%d G1a[i][j]=%f\n",i,j,G1a[i][j]);
            }
        }
    }*/
	for(int i=0; i<row; i++)
	  {
	  	if(abs(u[b[i]]+B[i])>epsilon)
	  	  mexPrintf("Problem at i=%d u[b[i]]=%f B[i]=%f\n",i,u[b[i]],B[i]);
	  }
}


void
SparseMatrix::Check_the_Solution(int periods, int y_kmin, int y_kmax, int Size, double *u, int *pivot, int* b)
{
	const double epsilon=1e-10;
	//std::map<std::pair<std::pair<int, int> ,int>, int> IM_i;
	Init(periods, y_kmin, y_kmax, Size, IM_i);
	NonZeroElem *first;
	int cal_y = y_kmin*Size;
	mexPrintf("     ");
	for(int i=0; i<Size; i++)
	  mexPrintf(" %8d",i);
	mexPrintf("\n");
	for(int t=y_kmin; t<periods+y_kmin; t++)
	  {
	  	mexPrintf("t=%5d",t);
      for(int i=0; i<Size; i++)
	  	   mexPrintf(" %d %1.6f",t*y_size+index_vara[i], y[t*y_size+index_vara[i]]);
			mexPrintf("\n");
	  }
	for(int i=0;i<Size*periods;i++)
	  {
	  	double res=0;
	  	int pos = pivot[i];
	  	mexPrintf("pos[%d]=%d",i,pos);
	  	int nb_var = At_Row(pos, &first);
	  	mexPrintf(" nb_var=%d\n",nb_var);
	  	for(int j=0;j<nb_var; j++)
	  	  {
	  	  	mexPrintf("(y[%d]=%f)*(u[%d]=%f)(r=%d, c=%d)\n",index_vara[first->c_index]+cal_y, y[index_vara[first->c_index]+cal_y], first->u_index, u[first->u_index], first->r_index, first->c_index);
	  	  	res += y[index_vara[first->c_index]+cal_y]*u[first->u_index];
	  	  	first=first->NZE_R_N;
	  	  }
			double tmp_=res;
			res += u[b[pos]];
			if(abs(res)>epsilon)
			  mexPrintf("Error for equation %d => res=%f y[%d]=%f u[b[%d]]=%f somme(y*u)=%f\n",pos,res,pos,y[index_vara[pos]], pos, u[b[pos]], tmp_);
	  }
	filename+=" stopped";
	mexErrMsgTxt(filename.c_str());
}





int
SparseMatrix::simulate_NG1(int blck, int y_size, int it_, int y_kmin, int y_kmax, int Size, int periods, bool print_it, bool cvg, int &iter)
{
  /*Triangularisation at each period of a block using a simple gaussian Elimination*/
  t_save_op_s *save_op_s;
  bool record=false;
  int *save_op=NULL, *save_opa=NULL, *save_opaa=NULL;
  long int nop=0, nopa=0;
  int tbreak=0, last_period=periods;
  int i,j,k;
  int pivj=0, pivk=0;
  int row, nb_var_piv, nb_var_sub, l_sub, sub_c_index, tmp_lag, l_piv, piv_c_index, tmp_u_count, lag;
  NonZeroElem *first, *firsta, *first_sub, *first_piv, *first_suba;
  double piv_abs, first_elem;
  if(start_compare==0)
    start_compare=y_kmin;;
#ifdef RECORD_ALL
  int save_u_count=u_count;
#endif
  u_count_alloc_save=u_count_alloc;
#ifdef PROFILER
  long int ndiv=0, nsub=0, ncomp=0, nmul=0;
  double tinsert=0, tdelete=0, tpivot=0, tbigloop=0;
  clock_t td1;
  int nbpivot=0, nbdiv=0, nbless=0, nbpivot_it=0, nbRealloc=0, ninsert=0;
#endif
  //pctimer_t t01;
  clock_t t01;
  //pctimer_t t1 = pctimer();
  clock_t t1 = clock();
#ifdef PROFILER
  tdelete1=0; tdelete2=0; tdelete21=0; tdelete22=0; tdelete221=0; tdelete222=0;
#endif
  nop1 = 0;
  if (iter>0)
    {
      mexPrintf("Sim : %f ms\n",(1000.0*(double(clock())-double(time00)))/double(CLOCKS_PER_SEC));
      mexEvalString("drawnow;");
    }
#ifdef MEMORY_LEAKS
  mexEvalString("feature('memstats');");
#endif
  if (isnan(res1) || isinf(res1))
    {
    	if(iter==0)
    	  {
    	  	for(j=0;j<y_size; j++)
            mexPrintf("variable %d at time %d = %f, %f\n",j+1, it_, y[j+it_*y_size], y[j+(it_+1)*y_size]);
    	  	mexPrintf("The initial values of endogenous variables are too far from the solution.\n");
    	  	mexPrintf("Change them!\n");
    	  	mexEvalString("drawnow;");
          mexEvalString("st=fclose('all');clear all;");
          filename+=" stopped";
          mexErrMsgTxt(filename.c_str());
    	  }
      if (slowc_save<1e-8)
        {
          mexPrintf("slowc_save=%g\n", slowc_save);
          for(j=0;j<y_size; j++)
            mexPrintf("variable %d at time %d = %f, %f\n",j+1, it_, y[j+it_*y_size], y[j+(it_+1)*y_size]);
          mexPrintf("Dynare cannot improve the simulation in block %d at time %d (variable %d)\n", blck+1, it_+1, max_res_idx);
          mexEvalString("drawnow;");
          mexEvalString("st=fclose('all');clear all;");
          filename+=" stopped";
          mexErrMsgTxt(filename.c_str());
#ifdef DEBUG_EX
          exit(-1);
#endif
        }
      slowc_save/=2;
      mexPrintf("Error: Simulation diverging, trying to correct it using slowc=%f\n",slowc_save);
      for (i=0;i<y_size*(periods+y_kmin);i++)
        y[i]=ya[i]+slowc_save*direction[i];
      iter--;
      return(0);
    }
  u_count+=u_count_init;
  if (alt_symbolic && alt_symbolic_count<alt_symbolic_count_max)
    {
      mexPrintf("Pivoting method will be applied only to the first periods.\n");
      alt_symbolic=false;
      symbolic=true;
      markowitz_c=markowitz_c_s;
      alt_symbolic_count++;
    }
  if (((res1/res1a-1)>-0.3) && symbolic && iter>0)
    {
    	if(restart>2)
          {
            mexPrintf("Divergence or slowdown occured during simulation.\nIn the next iteration, pivoting method will be applied to all periods.\n");
            symbolic=false;
            alt_symbolic=true;
            markowitz_c_s=markowitz_c;
            markowitz_c=0;
          }
      else
        {
          mexPrintf("Divergence or slowdown occured during simulation.\nIn the next iteration, pivoting method will be applied for a longer period.\n");
          start_compare=min(tbreak_g,periods);
          restart++;
        }
    }
	else
	  {
	    start_compare=y_kmin;
	    restart = 0;
	  }
  res1a=res1;



  if(print_it)
    {
      mexPrintf("-----------------------------------\n");
      mexPrintf("      Simulate     iteration° %d     \n",iter+1);
      mexPrintf("      max. error=%.10e       \n",double(max_res));
      mexPrintf("      sqr. error=%.10e       \n",double(res2));
      mexPrintf("      abs. error=%.10e       \n",double(res1));
      mexPrintf("-----------------------------------\n");
    }
	//Print(Size, b);
  if (cvg)
    {
      /*mexPrintf("End of simulate_NG1\n");
      mexEvalString("drawnow;");*/
      return(0);
    }
#ifdef PRINT_OUT
  mexPrintf("Size=%d y_size=%d y_kmin=%d y_kmax=%d u_count=%d u_count_alloc=%d periods=%d\n",Size,y_size,y_kmin,y_kmax,u_count,u_count_alloc,periods);
  mexEvalString("drawnow;");
#endif
#ifdef PROFILER
  clock_t t00 = clock();
#endif
#ifdef WRITE_u
  fstream toto;
  int i_toto=0;
  if (!symbolic)
    {
      toto.open("compare_good.txt", std::ios::out);
    }
#endif
#ifdef PROFILER
  t01=clock();
#endif
#ifdef PROFILER
  mexPrintf("initialization time=%f ms\n",1000.0*(double(t01)-double(t00))/double(CLOCKS_PER_SEC));
  mexEvalString("drawnow;");
#endif

#ifdef PRINT_OUT
  mexPrintf("sizeof(NonZeroElem)=%d\n",sizeof(NonZeroElem));
  /*for (i=0;i<Size*periods;i++)
    mexPrintf("b[%d]=%f\n",i,double(b[i]));*/
#endif
#ifdef RECORD_ALL
  if (record_all && !save_op_all)
    {
      nopa_all=(Size*periods)*Size*periods*4;  /*Initial guess on the total number of operations*/
#ifdef MEM_ALLOC_CHK
      mexPrintf("record_all save_op_all=(int*)mxMalloc(%d*sizeof(int))\n",nopa_all);
#endif
      save_op_all=(int*)mxMalloc(nopa_all*sizeof(int));
#ifdef MEM_ALLOC_CHK
      mexPrintf("ok\n");
#endif
      nop_all=0;
    }
#endif
  /*Add the first and the last values of endogenous to the exogenous */
#ifdef PROFILER
  t00 = clock();
#endif
#ifdef RECORD_ALL
  if (record_all && nop_all)
    {
//#ifdef PRINT_OUT
      mexPrintf("ShortInit\n");
      mexEvalString("drawnow;");
//#endif
      ShortInit(periods, y_kmin, y_kmax, Size, IM_i);
//#ifdef PRINT_OUT
      mexPrintf("run_triangular\n");
      mexEvalString("drawnow;");
//#endif
      run_triangular(nop_all,save_op_all);
//#ifdef PRINT_OUT
      mexPrintf("OK\n");
      mexEvalString("drawnow;");
//#endif
    }
  else
#endif
    if (g_nop_all>0)
      {
#ifdef PRINT_OUT
        mexPrintf("run_triangular\n");
        mexEvalString("drawnow;");
#endif
        run_u_period1(periods);
#ifdef PRINT_OUT
        mexPrintf("done\n");
        mexEvalString("drawnow;");
#endif
      }
    else
      {
#ifdef PRINT_OUT
        mexPrintf("Init\n");
        mexEvalString("drawnow;");
#endif
        Init(periods, y_kmin, y_kmax, Size, IM_i);
	      /*ua = (double*)mxMalloc(u_count_init * periods*sizeof(double));
	      for(i=0; i< u_count_init * periods;i++)
	        ua[i] = u[i];*/
#ifdef PRINT_OUT
        mexPrintf("done\n");
        mexEvalString("drawnow;");
#endif
        //Print(Size, b);

        //CheckIt(y_size, y_kmin, y_kmax, Size, periods, iter);

        //mexErrMsgTxt("Exit from Dynare");
        /*for(int i=0; i<row; i++)
	        {
	        	u[b[i]] = - u[b[i]];
	        }*/

        for (int t=0;t<periods;t++)
          {
          	//mexPrintf("t=%d periods=%d\n",t,periods);
#ifdef WRITE_u
            if (!symbolic && ((periods-t)<=y_kmax))
              {
                toto << "t=" << t << endl;
              }
#endif
            if (record && symbolic)
              {
                //nop*=8;
                if (save_op);
                {
                  mxFree(save_op);
                  save_op=NULL;
                }
#ifdef MEM_ALLOC_CHK
                mexPrintf("save_op=(int*)mxMalloc(%d*sizeof(int))\n",nop);
#endif
#ifdef N_MX_ALLOC
                save_op=malloc_std(nop);
#else
                save_op=(int*)mxMalloc(nop*sizeof(int));
#endif
                nopa=nop;
#ifdef MEM_ALLOC_CHK
                mexPrintf("ok\n");
#endif
              }
            nop=0;
#ifdef PRINT_OUT
            mexPrintf("---------------------------------------------------------\n");
            mexPrintf("t=%d\n",t);
            mexPrintf("---------------------------------------------------------\n");
#endif
            Clear_u();
            int ti=t*Size;
#ifdef PROFILER
            if (t<=start_compare)
              nbpivot_it=0;
#endif
            for (i=ti;i<Size+ti;i++)
              {
                /*finding the max-pivot*/
#ifdef PRINT_OUT
                Print(Size,b);
                mexPrintf("*************************************\n");
                mexPrintf("Finding the Pivot at column i=%d\n",i);
#endif
#ifdef PROFILER
                clock_t td0=clock();
#endif
                double piv=piv_abs=0;
                int nb_eq=At_Col(i, 0, &first);
#ifdef PRINT_OUT
                mexPrintf("nb_eq=%d\n",nb_eq);
#endif
                //mexPrintf("symbolic=%d t=%d start_compare=%d\n",symbolic, t, start_compare);
                if ((symbolic && t<=start_compare) || !symbolic)
                  {
#ifdef MARKOVITZ
                    double piv_v[Size];
                    int pivj_v[Size], pivk_v[Size], NR[Size], l=0, N_max=0;
                    bool one=false;
                    piv_abs=0;
#endif
                    for (j=0;j<nb_eq;j++)
                      {
                        //mexPrintf("j=%d\n",j);
#ifdef PRINT_OUT
                        mexPrintf("first=%x\n",first);
                        mexPrintf("examine col %d row %d with lag 0 line_done=%d \n",i,first->r_index,line_done[first->r_index]);
#endif
                        if (!line_done[first->r_index])
                          {
                            k=first->u_index;
#ifdef PRINT_OUT
                            mexPrintf("u[%d]=%f fabs(u[%d])=%f\n",k,double(u[k]),k,double(fabs(u[k])));
#endif
                            int jj=first->r_index;
                            int NRow_jj=NRow(jj);
#ifdef PROFILER
                            nbpivot++;
                            nbpivot_it++;
#endif

#ifdef MARKOVITZ
                            piv_v[l]=u[k];
                            double piv_fabs=fabs(u[k]);
                            pivj_v[l]=jj;
                            pivk_v[l]=k;
                            NR[l]=NRow_jj;
                            if (NRow_jj==1 && !one)
                              {
                                one=true;
                                piv_abs=piv_fabs;
                                N_max=NRow_jj;
                              }
                            if (!one)
                              {
                                if (piv_fabs>piv_abs)
                                  piv_abs=piv_fabs;
                                if (NRow_jj>N_max)
                                  N_max=NRow_jj;
                              }
                            else
                              {
                                if (NRow_jj==1)
                                  {
                                    if (piv_fabs>piv_abs)
                                      piv_abs=piv_fabs;
                                    if (NRow_jj>N_max)
                                      N_max=NRow_jj;
                                  }
                              }
                            l++;
#else
                            if (piv_abs<fabs(u[k])||NRow_jj==1)
                              {
                                piv=u[k];
                                piv_abs=fabs(piv);
                                pivj=jj;   //Line number
                                pivk=k;   //position in u
                                if (NRow_jj==1)
                                  break;
                              }
#endif
                          }
                        first=first->NZE_C_N;
                      }
#ifdef MARKOVITZ
                    double markovitz=0, markovitz_max=-9e70;
                    if (!one)
                      {
                        for (j=0;j<l;j++)
                          {
                            markovitz=exp(log(fabs(piv_v[j])/piv_abs)-markowitz_c*log(double(NR[j])/double(N_max)));
                            if (markovitz>markovitz_max)
                              {
                                piv=piv_v[j];
                                pivj=pivj_v[j];   //Line number
                                pivk=pivk_v[j];   //positi
                                markovitz_max=markovitz;
                              }
                          }
                      }
                    else
                      {
                        for (j=0;j<l;j++)
                          {
                            markovitz=exp(log(fabs(piv_v[j])/piv_abs)-markowitz_c*log(NR[j]/N_max));
                            if (markovitz>markovitz_max && NR[j]==1)
                              {
                                piv=piv_v[j];
                                pivj=pivj_v[j];   //Line number
                                pivk=pivk_v[j];   //positi
                                markovitz_max=markovitz;
                              }
                          }
                      }
#endif
#ifdef PROFILER
                    tpivot+=clock()-td0;
#endif

#ifdef PRINT_OUT
                    mexPrintf("Thats the pivot: %d with value %f in u[%d] \n",pivj,double(piv),pivk);
                    mexPrintf("______________________________________________\n");
                    mexPrintf("pivot[%d]=%d\n",i,pivj);
                    mexEvalString("drawnow;");
#endif
										if(iter>0 && t>start_compare)
										  {
										  	if(pivot_save[i-Size]+Size!=pivj)
										  	  mexPrintf("At t=%d at line i=%d pivj=%d and pivot_save[i-Size]+Size=%d\n",t,i,pivj, pivot_save[i-Size]+Size);
										  }
                    pivot[i]=pivj;
                    pivot_save[i]=pivj;
                    pivotk[i]=pivk;
                    pivotv[i]=piv;
                  }
                else
                  {
                    pivj=pivot[i-Size]+Size;
                    pivot[i]=pivj;
                    At_Pos(pivj, i, &first);
                    pivk=first->u_index;
                    piv=u[pivk];
                    piv_abs=fabs(piv);
                  }
                line_done[pivj]=true;
#ifdef PRINT_u
                mexPrintf("FLD u[%d] (%f=%f)   |",pivk,u[pivk],piv);
                Print_u();mexPrintf("\n");
                mexEvalString("drawnow;");
#endif
                if (symbolic)
                  {
                    if (record)
                      {
                        if (nop+1>=nopa)
                          {
                            nopa=int(1.5*nopa);
#ifdef MEM_ALLOC_CHK
                            mexPrintf("save_op=(int*)mxRealloc(save_op,%d*sizeof(int))\n",nopa);
#endif
#ifdef PROFILER
                            nbRealloc++;
#endif
#ifdef N_MX_ALLOC
                            save_op=realloc_std(save_op, nopa);
#else
                            save_op=(int*)mxRealloc(save_op,nopa*sizeof(int));
#endif
#ifdef MEM_ALLOC_CHK
                            mexPrintf("ok\n");
#endif
                          }
                        save_op_s=(t_save_op_s*)(&(save_op[nop]));
                        save_op_s->operat=IFLD;
                        save_op_s->first=pivk;
                        save_op_s->lag=0;
                      }
                    nop+=2;
                  }
#ifdef RECORD_ALL
                else if (record_all)
                  {
                    if (nop_all+1>=nopa_all)
                      chk_avail_mem(&save_op_all,&nop_all,&nopa_all,1,t);
                    save_op_all[nop_all]=IFLD;
                    save_op_all[nop_all+1]=pivk;
                    nop_all+=2;
                  }
#endif
								/*mexPrintf("piv_abs=%f\n",piv_abs);
								mexEvalString("drawnow;");*/
                if (piv_abs<eps)
                  {
                    mexPrintf("Error: singular system in Simulate_NG1\n");
                    mexEvalString("drawnow;");
                    mexEvalString("st=fclose('all');clear all;");
                    filename+=" stopped";
                    mexErrMsgTxt(filename.c_str());
                  }
                /*divide all the non zeros elements of the line pivj by the max_pivot*/
                int nb_var=At_Row(pivj,&first);
                NonZeroElem* bb[nb_var];
                /*mexPrintf("nb_var=%d\n",nb_var);
                mexEvalString("drawnow;");*/
                for(j=0;j<nb_var;j++)
                  {
                    bb[j]=first;
                    /*mexPrintf("j=%d",j);
                    mexPrintf(" first->NZE_R_N=%x\n",first->NZE_R_N);*/
                    first=first->NZE_R_N;
                  }

#ifdef PRINT_OUT
                mexPrintf("nb_var=%d\n",nb_var);
#endif
                for (j=0;j<nb_var;j++)
                  {
                    first=bb[j];
#ifdef PRINT_OUT
                    mexPrintf("j=%d ",j);
                    mexPrintf("first=%x ",first);
                    mexPrintf("dividing at lag %d [%d, %d] u[%d]\n",first->lag_index, first->r_index, first->c_index, first->u_index);
#endif
                    u[first->u_index]/=piv;
#ifdef PROFILER
                    nbdiv++;
#endif
#ifdef WRITE_u
                    if ((periods-t)<=y_kmax)
                      {
                        toto << i_toto << " u[" /*<< first->u_index*/ << "]/=" << piv << "=" << u[first->u_index] << endl;
                        i_toto++;
                      }
#endif
#ifdef PRINT_u
                    mexPrintf("FDIV u[%d](%f)/=piv(%f)=(%f)   |",first->u_index,u[first->u_index],piv,u[first->u_index]);
                    Print_u();mexPrintf("\n");
#endif
                    if (symbolic)
                      {
                        if (record)
                          {
                            if (nop+j*2+1>=nopa)
                              {
                                nopa=int(1.5*nopa);
#ifdef MEM_ALLOC_CHK
                                mexPrintf("save_op=(int*)mxRealloc(save_op,%d*sizeof(int))\n",nopa);
#endif
#ifdef PROFILER
                                nbRealloc++;
#endif
#ifdef N_MX_ALLOC
                                save_op=realloc_std(save_op, nopa);
#else
                                save_op=(int*)mxRealloc(save_op,nopa*sizeof(int));
#endif
#ifdef MEM_ALLOC_CHK
                                mexPrintf("ok\n");
#endif
                              }
                            save_op_s=(t_save_op_s*)(&(save_op[nop+j*2]));
                            //mexPrintf("save_op[%d] : operat=%d, first=%d, lag=%d omp_get_thread_num()=%d\n",nop+j*2, IFDIV, first->u_index, first->lag_index, omp_get_thread_num());
                            save_op_s->operat=IFDIV;
                            save_op_s->first=first->u_index;
                            save_op_s->lag=first->lag_index;
                          }
                        //nop+=2; ///!!
                      }
#ifdef RECORD_ALL
                    else if (record_all)
                      {
                        if (nop_all+1>=nopa_all)
                          chk_avail_mem(&save_op_all,&nop_all,&nopa_all,1,t);
                        save_op_all[nop_all]=IFDIV;
                        save_op_all[nop_all+1]=first->u_index;
                        nop_all+=2;
                      }
#endif
                    //first=first->NZE_R_N;
                  }
                nop += nb_var*2;
#ifdef PRINT_OUT
                mexPrintf("dividing at u[%d]\n",b[pivj]);
#endif
                u[b[pivj]]/=piv;
#ifdef WRITE_u
                if ((periods-t)<=y_kmax)
                  {
                    toto << i_toto << " u[" /*<< b[pivj]*/ << "]/=" << piv << "=" << u[b[pivj]] << endl;
                    i_toto++;
                  }
#endif
#ifdef PRINT_u
                mexPrintf("FDIV u[%d](%f)/=piv(%f)=(%f)   |",b[pivj],u[b[pivj]],piv,u[b[pivj]]);
                Print_u();mexPrintf("\n");
#endif
                if (symbolic)
                  {
                    if (record)
                      {
                        if (nop+1>=nopa)
                          {
                            nopa=int(1.5*nopa);
#ifdef MEM_ALLOC_CHK
                            mexPrintf("save_op=(int*)mxRealloc(save_op,%d*sizeof(int))\n",nopa);
#endif
#ifdef PROFILER
                            nbRealloc++;
#endif
#ifdef N_MX_ALLOC
                            save_op=realloc_std(save_op, nopa);
#else
                            save_op=(int*)mxRealloc(save_op,nopa*sizeof(int));
#endif
#ifdef MEM_ALLOC_CHK
                            mexPrintf("ok\n");
#endif
                          }
                        save_op_s=(t_save_op_s*)(&(save_op[nop]));
                        save_op_s->operat=IFDIV;
                        save_op_s->first=b[pivj];
                        save_op_s->lag=0;
                      }
                    nop+=2;
                  }
#ifdef RECORD_ALL
                else if (record_all)
                  {
                    if (nop_all+1>=nopa_all)
                      chk_avail_mem(&save_op_all,&nop_all,&nopa_all,1,t);
                    save_op_all[nop_all]=IFDIV;
                    save_op_all[nop_all+1]=b[pivj];
                    nop_all+=2;
                  }
#endif
                /*substract the elements on the non treated lines*/
                nb_eq=At_Col(i,&first);
                NonZeroElem *first_piva;
                int nb_var_piva=At_Row(pivj,&first_piva);
                //mexPrintf("pivj=%d\n",pivj);
#ifdef PRINT_OUT
                if(iter>0)
                  {
                mexPrintf("ok4 nb_eq=%d iter=%d\n",nb_eq,iter);
                mexEvalString("drawnow;");
                  }
#endif

                NonZeroElem* bc[nb_eq];
                for(j=0;j<nb_eq;j++)
                  {
                    bc[j]=first;
                    first=first->NZE_C_N;
                  }
								//#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) private(first, row, first_elem, nopa, save_op_s, nb_var_piv, nb_var_piva, first_piv, first_piva, first_sub, nb_var_sub, l_sub, l_piv, sub_c_index, piv_c_index, tmp_lag)
                for (j=0;j<nb_eq;j++)
                  {
                    first=bc[j];
                    row=first->r_index;
                    /*mexPrintf("-------------------\n");
                  	mexPrintf("j=%d line_done[row=%d]=%d\n",j,row, line_done[row]);*/
#ifdef PRINT_OUT
                    mexPrintf("t=%d, j=%d, line_done[%d]=%hd process=%d\n", t, j, row, line_done[row],omp_get_thread_num());
#endif
                    if (!line_done[row])
                      {
#ifdef PRINT_OUT
                        mexPrintf("Substracting from line %d lag %d\n",row,first->lag_index);
#endif
                        first_elem=u[first->u_index];
#ifdef PRINT_u
                        mexPrintf("FLD u[%d] (%f)  |",first->u_index,u[first->u_index]);
                        Print_u();mexPrintf("\n");
#endif
                        if (symbolic)
                          {
                            if (record)
                              {
                                if (nop+1>=nopa)
                                  {
                                    nopa=int(1.5*nopa);
#ifdef MEM_ALLOC_CHK
                                    mexPrintf("save_op=(int*)mxRealloc(save_op,%d*sizeof(int))\n",nopa);
#endif
#ifdef PROFILER
                                    nbRealloc++;
#endif
#ifdef N_MX_ALLOC
                                    save_op=realloc_std(save_op, nopa);
#else
                                    save_op=(int*)mxRealloc(save_op,nopa*sizeof(int));
#endif
#ifdef MEM_ALLOC_CHK
                                    mexPrintf("ok\n");
#endif
                                  }
                                save_op_s=(t_save_op_s*)(&(save_op[nop]));
                                save_op_s->operat=IFLD;
                                save_op_s->first=first->u_index;
                                save_op_s->lag=abs(first->lag_index);
                              }
                            nop+=2;
                          }
#ifdef RECORD_ALL
                        else if (record_all)
                          {
                            if (nop_all+1>=nopa_all)
                              chk_avail_mem(&save_op_all,&nop_all,&nopa_all,1,t);
                            save_op_all[nop_all]=IFLD;
                            save_op_all[nop_all+1]=first->u_index;
                            nop_all+=2;
                          }
#endif
                        /*mexPrintf("For equ=9\n");
                        int nb_var__=At_Row(9,&first_piv);
                        for(int uu=0; uu<nb_var__; uu++)
											    {
											    	mexPrintf("->   first_piv->c_index=%d\n",first_piv->c_index);
											    	first_piv=first_piv->NZE_R_N;
											    }

												first_piv = first_piva;
												mexPrintf("OK\n");
												for(int uu=0; uu<nb_var_piva; uu++)
											    {
											    	mexPrintf("->   first_piv->c_index=%d\n",first_piv->c_index);
											    	first_piv=first_piv->NZE_R_N;
											    }*/

                        nb_var_piv=nb_var_piva;
                        first_piv=first_piva;
                        nb_var_sub=At_Row(row,&first_sub);
                        l_sub=0;
                        l_piv=0;
                        sub_c_index=first_sub->c_index;
                        piv_c_index=first_piv->c_index;
                        tmp_lag=first_sub->lag_index;
#ifdef PROFILER
                        td1 = clock();
#endif
                        while (l_sub<nb_var_sub || l_piv<nb_var_piv)
                          {
#ifdef PRINT_OUT
                            if (l_piv<nb_var_piv)
                              mexPrintf(" piv eq=%d lag=%d var=%d l1=%d",first_piv->r_index, first_piv->lag_index, first_piv->c_index,l_piv);
                            mexPrintf("l_sub(%d)<nb_var_sub(%d)\n",l_sub,nb_var_sub);
                            if (l_sub<nb_var_sub)
                              mexPrintf(" sub eq=%d lag=%d var=%d l0=%d",first_sub->r_index, first_sub->lag_index, first_sub->c_index,l_sub);
#endif
                            //mexPrintf("sub_c_index=%d piv_c_index=%d, l_sub=%d nb_var_sub=%d, l_piv=%d nb_var_piv=%d\n",sub_c_index, piv_c_index, l_sub, nb_var_sub, l_piv, nb_var_piv);
                            if (l_sub<nb_var_sub && (sub_c_index<piv_c_index || l_piv>=nb_var_piv))
                              {
                              	//There is no nonzero element at line pivot for this column=> Nothing to do for the current element got to next column
                              	//mexPrintf("Nothing\n");
                                first_sub=first_sub->NZE_R_N;
                                if (first_sub)
                                  sub_c_index=first_sub->c_index;
                                else
                                  sub_c_index=Size*periods;
                                l_sub++;
                              }
                            else if (sub_c_index>piv_c_index || l_sub>=nb_var_sub)
                              {
                              	// There is an nonzero element at row pivot but not at the current row=> insert a negative element in the current row
                              	//mexPrintf("Insert\n");
                                tmp_u_count=Get_u();
#ifdef PROFILER
                                clock_t td0=clock();
#endif
                                lag=first_piv->c_index/Size-row/Size;
                                Insert(row,first_piv->c_index,tmp_u_count,lag);
#ifdef PROFILER
                                tinsert+=clock()-td0;
#endif
                                u[tmp_u_count]=-u[first_piv->u_index]*first_elem;
#ifdef PROFILER
                                nbless++;
                                ninsert++;
#endif
                                if (symbolic)
                                  {
                                    if (record)
                                      {
                                        if (nop+2>=nopa)
                                          {
                                            nopa=int(1.5*nopa);
#ifdef MEM_ALLOC_CHK
                                            mexPrintf("save_op=(int*)mxRealloc(save_op,%d*sizeof(int))\n",nopa);
#endif
#ifdef PROFILER
                                            nbRealloc++;
#endif
#ifdef N_MX_ALLOC
                                            save_op=realloc_std(save_op, nopa);
#else
                                            save_op=(int*)mxRealloc(save_op,nopa*sizeof(int));
#endif
#ifdef MEM_ALLOC_CHK
                                            mexPrintf("ok\n");
#endif
                                          }
                                        save_op_s=(t_save_op_s*)(&(save_op[nop]));
                                        save_op_s->operat=IFLESS;
                                        save_op_s->first=tmp_u_count;
                                        save_op_s->second=first_piv->u_index;
                                        save_op_s->lag=max(first_piv->lag_index,abs(tmp_lag));
                                      }
                                    nop+=3;
                                  }
#ifdef RECORD_ALL
                                else if (record_all)
                                  {
                                    if (nop_all+2>=nopa_all)
                                      chk_avail_mem(&save_op_all,&nop_all,&nopa_all,2,t);
                                    save_op_all[nop_all]=IFLESS;
                                    save_op_all[nop_all+1]=tmp_u_count;
                                    save_op_all[nop_all+2]=first_piv->u_index;
                                    nop_all+=3;
                                  }
#endif
                                first_piv=first_piv->NZE_R_N;
                                if (first_piv)
                                  piv_c_index=first_piv->c_index;
                                else
                                  piv_c_index=Size*periods;
                                l_piv++;
                              }
                            else /*first_sub->c_index==first_piv->c_index*/
                              {
                                if (i==sub_c_index)
                                  {
                                  	 //mexPrintf("Delete\n");
#ifdef PRINT_OUT
                                    /*if(iter>0)
                                      {
                                    mexPrintf("   delete element [%d, %d] lag %d u[%d]\n",first_sub->r_index,first_sub->c_index,first_sub->lag_index,first_sub->u_index);
                                        mexEvalString("drawnow;");
                                      }*/
#endif
                                    firsta=first;
                                    first_suba=first_sub->NZE_R_N;
#ifdef PROFILER
                                    clock_t td0=clock();
#endif
                                    Delete(first_sub->r_index,first_sub->c_index, Size);
#ifdef PROFILER
                                    tdelete+=clock()-td0;
#endif
                                    first=firsta->NZE_C_N;
                                    first_sub=first_suba;
                                    if (first_sub)
                                      sub_c_index=first_sub->c_index;
                                    else
                                      sub_c_index=Size*periods;
                                    l_sub++;
                                    first_piv=first_piv->NZE_R_N;
                                    if (first_piv)
                                      piv_c_index=first_piv->c_index;
                                    else
                                      piv_c_index=Size*periods;
                                    l_piv++;
#ifdef PRINT_OUT
                                    Print(Size,b);
#endif
                                  }
                                else
                                  {
                                  	//mexPrintf("Substract\n");
#ifdef PRINT_OUT
                                    mexPrintf("  u[%d]-=u[%d]*%f\n",first_sub->u_index,first_piv->u_index,double(first_elem));
#endif
                                    u[first_sub->u_index]-=u[first_piv->u_index]*first_elem;
#ifdef PROFILER
                                    nbless++;
#endif
#ifdef WRITE_u
                                    if ((periods-t)<=y_kmax)
                                      {
                                        toto << i_toto << " u[" /*<< first_sub->u_index*/ << "]-=u[" /*<< first_piv->u_index*/ << "]*" << first_elem << "=" << u[first_sub->u_index] << endl;
                                        i_toto++;
                                      }
#endif
#ifdef PRINT_u
                                    if(iter>0)
                                      {
                                    mexPrintf("FSUB u[%d]-=u[%d](%f)*r(%f)=(%f)  |",first_sub->u_index,first_piv->u_index,u[first_piv->u_index],first_elem,u[first_sub->u_index]);
                                    /*Print_u();*/mexPrintf("\n");
                                        mexEvalString("drawnow;");
                                      }
#endif
                                    if (symbolic)
                                      {
                                        if (record)
                                          {
                                            if (nop+2>=nopa)
                                              {
                                                nopa=int(1.5*nopa);
#ifdef MEM_ALLOC_CHK
                                                mexPrintf("save_op=(int*)mxRealloc(save_op,%d*sizeof(int))\n",nopa);
#endif
#ifdef PROFILER
                                                nbRealloc++;
#endif
#ifdef N_MX_ALLOC
                                                save_op=realloc_std(save_op, nopa);
#else
                                                save_op=(int*)mxRealloc(save_op,nopa*sizeof(int));
#endif
#ifdef MEM_ALLOC_CHK
                                                mexPrintf("ok\n");
#endif
                                              }
                                            save_op_s=(t_save_op_s*)(&(save_op[nop]));
                                            save_op_s->operat=IFSUB;
                                            save_op_s->first=first_sub->u_index;
                                            save_op_s->second=first_piv->u_index;
                                            save_op_s->lag=max(abs(tmp_lag),first_piv->lag_index);
                                          }
                                        nop+=3;
                                      }
#ifdef RECORD_ALL
                                    else if (record_all)
                                      {
                                        if (nop_all+2>=nopa_all)
                                          chk_avail_mem(&save_op_all,&nop_all,&nopa_all,2,t);
                                        save_op_all[nop_all]=IFSUB;
                                        save_op_all[nop_all+1]=first_sub->u_index;
                                        save_op_all[nop_all+2]=first_piv->u_index;
                                        nop_all+=3;
                                      }
#endif
                                    first_sub=first_sub->NZE_R_N;
                                    if (first_sub)
                                      sub_c_index=first_sub->c_index;
                                    else
                                      sub_c_index=Size*periods;
                                    l_sub++;
                                    first_piv=first_piv->NZE_R_N;
                                    if (first_piv)
                                      piv_c_index=first_piv->c_index;
                                    else
                                      piv_c_index=Size*periods;
                                    l_piv++;
                                  }
                              }
                          }
#ifdef PRINT_OUT
                        mexPrintf("row=%d pivj=%d\n",row,pivj);
                        mexPrintf("  u[%d](%f)-=u[%d](%f)*%f\n",b[row],double(u[b[row]]), b[pivj],double(u[b[pivj]]),double(first_elem));
#endif

                        u[b[row]]-=u[b[pivj]]*first_elem;
#ifdef PROFILER
                        nbless++;
                        tbigloop += clock() - td1;
#endif
#ifdef WRITE_u
                        if ((periods-t)<=y_kmax)
                          {
                            toto << i_toto << " u[" /*<< b[row]*/ << "]-=u[" /*<< b[pivj]*/ << "]*" << first_elem << "=" << u[b[row]] << endl;
                            i_toto++;
                          }
#endif
#ifdef PRINT_u
                        mexPrintf("FSUB u[%d]-=u[%d](%f)*r(%f)=(%f)   |",b[row],b[pivj],u[b[pivj]],first_elem,u[b[row]]);
                        Print_u();mexPrintf("\n");
#endif

                        if (symbolic)
                          {
                            if (record)
                              {
                                if (nop+2>=nopa)
                                  {
                                    nopa=int(1.5*nopa);
#ifdef MEM_ALLOC_CHK
                                    mexPrintf("save_op=(int*)mxRealloc(save_op,%d*sizeof(int))\n",nopa);
#endif
#ifdef PROFILER
                                    nbRealloc++;
#endif
#ifdef N_MX_ALLOC
                                    save_op=realloc_std(save_op, nopa);
#else
                                    save_op=(int*)mxRealloc(save_op,nopa*sizeof(int));
#endif
#ifdef MEM_ALLOC_CHK
                                    mexPrintf("ok\n");
#endif
                                  }
                                save_op_s=(t_save_op_s*)(&(save_op[nop]));
                                save_op_s->operat=IFSUB;
                                save_op_s->first=b[row];
                                save_op_s->second=b[pivj];
                                save_op_s->lag=abs(tmp_lag);
                              }
                            nop+=3;
                          }
#ifdef RECORD_ALL
                        else if (record_all)
                          {
                            if (nop_all+2>=nopa_all)
                              chk_avail_mem(&save_op_all,&nop_all,&nopa_all,2,t);
                            save_op_all[nop_all]=IFSUB;
                            save_op_all[nop_all+1]=b[row];
                            save_op_all[nop_all+2]=b[pivj];
                            nop_all+=3;
                          }
                        if(iter>0)
                          {
                            mexPrintf("ok4g j=%d\n",j);
                            mexEvalString("drawnow;");
                          }

#endif
                      }
                    /*else
                      first=first->NZE_C_N;*/

#ifdef PRINT_OUT
                    mexPrintf(" bef first=%x\n",first);
#endif
                  }
              }
            if (symbolic)
              {
#ifdef PROFILER
                td1=clock();
                mexPrintf("at %d nop=%d ?= nop1=%d record=%d save_opa=%x save_opaa=%x \n", t, nop, nop1, record, save_opa, save_opaa);
                mexEvalString("drawnow;");
#endif
                if (record && (nop==nop1))
                  {
#ifdef PRINT_OUT
                    mexPrintf("nop=%d, nop1=%d, record=%d save_opa=%x save_opaa=%x\n",nop, nop1, record, save_opa, save_opaa);
#endif
                    if (save_opa && save_opaa)
                      {
#ifdef PROFILER
                        clock_t ta=clock();
#endif
                        if (compare( save_op, save_opa, save_opaa, t, periods, nop, Size
#ifdef PROFILER
                        , &ndiv, &nsub
#endif
                        ))
                          {
#ifdef PROFILER
                            mexPrintf("done OK\n");
                            mexPrintf("t=%d over periods=%d t0=%f t1=%f y_kmin=%d y_kmax=%d\n",t,periods,1000*(ta-t00), 1000.0*(double(clock())-double(ta))/double(CLOCKS_PER_SEC),y_kmin, y_kmax);
                            mexPrintf("compare time %f ms\n",1000.0*(double(clock())-double(ta))/double(CLOCKS_PER_SEC));
                            mexEvalString("drawnow;");
#endif
                            tbreak=t;
                            tbreak_g=tbreak;
                            break;
                          }
                      }
                    if (save_opa)
                      {
                        if (save_opaa)
                          {
#ifdef N_MX_ALLOC
                            free(save_opaa);
#else
                            mxFree(save_opaa);
#endif
                            save_opaa=NULL;
                          }
#ifdef MEM_ALLOC_CHK
                        mexPrintf("save_opaa=(int*)mxMalloc(%d*sizeof(int))\n",nop1);
#endif
#ifdef N_MX_ALLOC
                        save_opaa=malloc_std(nop1);
#else
                        save_opaa=(int*)mxMalloc(nop1*sizeof(int));
#endif
#ifdef MEM_ALLOC_CHK
                        mexPrintf("ok\n");
#endif
                        memcpy(save_opaa, save_opa, nop1*sizeof(int));
                      }
                    if (save_opa)
                      {
#ifdef N_MX_ALLOC
                        free(save_opa);
#else
                        mxFree(save_opa);
#endif
                        save_opa=NULL;
                      }
#ifdef MEM_ALLOC_CHK
                    mexPrintf("save_opa=(int*)mxMalloc(%d*sizeof(int))\n",nop);
#endif
#ifdef N_MX_ALLOC
                    save_opa=malloc_std(nop);
#else
                    save_opa=(int*)mxMalloc(nop*sizeof(int));
#endif
#ifdef MEM_ALLOC_CHK
                    mexPrintf("ok\n");
#endif
                    memcpy(save_opa, save_op, nop*sizeof(int));
                  }
                else
                  {
                    if (nop==nop1)
                      record=true;
                    else
                      {
                        record=false;
                        if (save_opa)
                          {
#ifdef N_MX_ALLOC
                            free(save_opa);
#else
                            mxFree(save_opa);
#endif
                            save_opa=NULL;
                          }
                        if (save_opaa)
                          {
#ifdef N_MX_ALLOC
                            free(save_opaa);
#else
                            mxFree(save_opaa);
#endif
                            save_opaa=NULL;
                          }
                      }
                  }
                nop2=nop1;
                nop1=nop;
#ifdef PROFILER
                tcompare+=clock()-td1;
#endif
              }
            record=true;
          }
      }
  nop_all+=nop;
#ifdef PROFILER
  t01=clock();
  mexPrintf("resolution time=%f ms ndiv=%d nsub=%d ncomp=%d \n",1000.0*(double(t01)-double(t00))/double(CLOCKS_PER_SEC),ndiv,nsub,ncomp);
  mexPrintf(" tinsert=%f tdelete=%f tpivot=%f tbigloop=%f tdelete1=%f tdelete2=%f \n",double(1000*tinsert), double(1000*tdelete), double(1000*tpivot), double(1000*tbigloop), double(1000*tdelete1), double(1000*tdelete2));
  mexPrintf(" tdelete21=%f tdelete22=%f \n",double(1000*tdelete21),double(1000*tdelete22));
  mexPrintf(" tdelete221=%f tdelete222=%f \n",double(1000*tdelete221),double(1000*tdelete222));
  mexPrintf(" tcompare=%f \n",double(1000*tcompare));
  mexPrintf(" ninsert=%d\n",ninsert);
  mexPrintf("nbpivot=%d, nbdiv=%d, nbless=%d, nop=%d nbpivot_it=%d nbRealloc=%d\n", nbpivot, nbdiv, nbless, nbpivot + nbdiv + nbless, nbpivot_it, nbRealloc);
  mexEvalString("drawnow;");
#endif
  if (symbolic)
    {
#ifdef N_MX_ALLOC
      if (save_op)
        free(save_op);
      if (save_opa)
        free(save_opa);
      if (save_opaa)
        free(save_opaa);
#else
      if (save_op)
        mxFree(save_op);
      if (save_opa)
        mxFree(save_opa);
      if (save_opaa)
        mxFree(save_opaa);
#endif
    }
  close_swp_file();

  /*The backward substitution*/
#ifdef PRINT_OUT
  Print(Size,b);
#endif
  double slowc_lbx=slowc, res1bx;
#ifdef PROFILER
  t00 = clock();
#endif
  for (i=0;i<y_size*(periods+y_kmin);i++)
    ya[i]=y[i];
  slowc_save=slowc;
  res1bx=bksub( tbreak, last_period, Size, slowc_lbx
#ifdef PROFILER
  , &nmul
#endif
  );
  //t01=pctimer();
  t01=clock();
#ifdef PROFILER
  mexPrintf("backward substitution time=%f ms nmul=%d\n",1000.0*(double(t01)-double(t00))/double(CLOCKS_PER_SEC),nmul);
  mexEvalString("drawnow;");
#endif
#ifdef RECORD_ALL
  if (!((record_all && nop_all)||g_nop_all>0))
    {
      u_count=save_u_count;

    }
#endif
  End(Size);
  if (print_it)
    {
      //pctimer_t t2 = pctimer();
      clock_t t2 = clock();
      mexPrintf("(** %f milliseconds **)\n", 1000.0*(double(t2) - double(t1))/double(CLOCKS_PER_SEC));
      mexEvalString("drawnow;");
    }

#ifdef WRITE_u
  if (!symbolic)
    {
      toto.close();
      mexEvalString("st=fclose('all');clear all;");
      filename+=" stopped";
      mexErrMsgTxt(filename.c_str());
    }
#endif
  close_swp_file();
  time00=clock();
  if(tbreak_g==0)
    tbreak_g=periods;

	/*Check the solution*/
	/*Check_the_Solution(periods, y_kmin, y_kmax, Size, ua, pivot, b);
	mxFree(ua);*/

  return(0);
}

void
SparseMatrix::fixe_u(double **u, int u_count_int, int max_lag_plus_max_lead_plus_1)
{
  //mexPrintf("u_count_int=%d\n",u_count_int);
  u_count=u_count_int * periods;
  u_count_alloc = 2*u_count;
  //mexPrintf("mxMalloc(%d*sizeof(double)=%d)\n",u_count_alloc,u_count_alloc*sizeof(double));
  (*u)=(double*)mxMalloc(u_count_alloc*sizeof(double));
  memset((*u), 0, u_count_alloc*sizeof(double));
  u_count_init=max_lag_plus_max_lead_plus_1;
}

