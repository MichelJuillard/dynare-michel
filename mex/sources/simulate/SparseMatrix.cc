/*
 * Copyright (C) 2007-2008 Dynare Team
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
}





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

void SparseMatrix::Delete(int r,int c, int Size, int *b)
{
  NonZeroElem *first=FNZE_R[r], *firsta=NULL;
#ifdef PROFILER
  clock_t td0, td1, td2;
  td0=clock();
#endif
  while (first->c_index!=c /*&& first->NZE_R_N*/)
    {
#ifdef PRINT_OUT
      mexPrintf("looking not CRS [%d, %d]\n",first->r_index,first->c_index);
      mexEvalString("drawnow;");
#endif
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
  while (first->r_index!=r /*&& first->NZE_C_N*/)
    {
#ifdef PRINT_OUT
      mexPrintf("looking not CSS [%d, %d]\n",first->r_index,first->c_index);
      mexEvalString("drawnow;");
#endif
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



void SparseMatrix::Insert(int r, int c, int u_index, int lag_index)
{
#ifdef PRINT_OUT
  mexPrintf("In Insert r=%d, c=%d, u=%d, lag=%d \n",r,c,u_index,lag_index);
#endif
  NonZeroElem *first=FNZE_R[r], *firsta=NULL, *firstn=NULL;
  if (first)
    {
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
      /*if (first->c_index!=c)
        {*/
#ifdef PRINT_OUT
      mexPrintf("retain first->c_index=%d c=%d\n",first->c_index,c);
#endif
#ifdef NEW_ALLOC
      firstn=mem_mngr.mxMalloc_NZE();
#else
      firstn=(NonZeroElem*)mxMalloc(sizeof(NonZeroElem));
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
          first->NZE_R_N=firstn;
          firstn->NZE_R_N=NULL;
        }
      NbNZRow[r]++;

      /*}
      else
      mexPrintf("Error (in Insert): in CRS element r=%, c=%d already exists\n",r,c);*/
    }
  else
    {
#ifdef NEW_ALLOC
      firstn=mem_mngr.mxMalloc_NZE();
#else
      firstn=(NonZeroElem*)mxMalloc(sizeof(NonZeroElem));
#endif
      firstn->u_index=u_index;
      firstn->r_index=r;
      firstn->c_index=c;
      firstn->lag_index=lag_index;
      FNZE_R[r]=firstn;
      firstn->NZE_R_N=first;
      NbNZRow[r]++;
    }
  first=FNZE_C[c];
  firsta=NULL;
  while (first->r_index<r && first->NZE_C_N)
    {
      firsta=first;
      first=first->NZE_C_N;
    }
  /*if (first->r_index!=r)
    {*/
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
      first->NZE_C_N=firstn;
      firstn->NZE_C_N=NULL;
    }
  NbNZCol[c]++;
  /*}
  else
  mexPrintf("Error (in Insert): in CCS element r=%, c=%d already exists\n",r,c);*/
}

void SparseMatrix::Read_SparseMatrix(string file_name, int Size, int periods, int y_kmin, int y_kmax)
{
  int i,j,eq,var,lag;
  /*mexPrintf("in Read_SparseMatrix\n");
  mexPrintf("file_name=%s SaveCode.is_open()=%d\n",file_name.c_str(), SaveCode.is_open());
  mexEvalString("drawnow;");*/
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
  /*mexPrintf("u_count_init=%d\n",u_count_init);
  mexEvalString("drawnow;");*/
  for (i=0;i<u_count_init;i++)
    {
      /*mexPrintf("i=%d\n",i);
      mexEvalString("drawnow;");*/
      SaveCode.read(reinterpret_cast<char *>(&eq), sizeof(eq));
      SaveCode.read(reinterpret_cast<char *>(&var), sizeof(var));
      SaveCode.read(reinterpret_cast<char *>(&lag), sizeof(lag));
      SaveCode.read(reinterpret_cast<char *>(&j), sizeof(j));
      /*mexPrintf("0i=%d eq=%d var=%d lag=%d j=%d\n",i,eq,var,lag,j );
      mexEvalString("drawnow;");*/
      IM_i[std::make_pair(std::make_pair(eq, var), lag)] = j;
      /*mexPrintf("1i=%d\n",i);
      mexEvalString("drawnow;");*/
    }
#ifdef MEM_ALLOC_CHK
  mexPrintf("index_vara=(int*)mxMalloc(%d*sizeof(int))\n",Size*(periods+y_kmin+y_kmax));
#endif
  /*mexPrintf("Size=%d periods=%d y_kmin=%d y_kmax=%d\n",Size,periods, y_kmin, y_kmax);
  mexEvalString("drawnow;");*/
  index_vara=(int*)mxMalloc(Size*(periods+y_kmin+y_kmax)*sizeof(int));
#ifdef MEM_ALLOC_CHK
  mexPrintf("ok\n");
#endif
  for (j=0;j<Size;j++)
    {
      SaveCode.read(reinterpret_cast<char *>(&index_vara[j]), sizeof(*index_vara));
      //mexPrintf("index_vara[%d]=%d\n",j,index_vara[j]);
    }
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
  index_equa=(int*)mxMalloc(Size*sizeof(int));
  for(j=0;j<Size;j++)
    {
      SaveCode.read(reinterpret_cast<char *>(&index_equa[j]), sizeof(*index_equa));
      //mexPrintf("index_equa[%d]=%d\n",j,index_equa[j]);
    }
}


void SparseMatrix::Init(int periods, int y_kmin, int y_kmax, int Size, std::map<std::pair<std::pair<int, int> ,int>, int> IM)
{
  int t,i, eq, var, lag;
  double tmp_b=0.0;
  std::map<std::pair<std::pair<int, int> ,int>, int>::iterator it4;
  NonZeroElem* first;
  //mexPrintf("periods=%d, y_kmin=%d, y_kmax=%d, Size=%d, IM.size()=%d\n",periods, y_kmin, y_kmax, Size, IM.size());
#ifdef MEM_ALLOC_CHK
  mexPrintf("pivot=(int*)mxMalloc(%d*sizeof(int))\n",Size*periods);
#endif
  pivot=(int*)mxMalloc(Size*periods*sizeof(int));
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
  memset(line_done, 0, periods*Size*sizeof(*line_done));
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
  memset(FNZE_R, 0, i);
  memset(FNZE_C, 0, i);
#ifdef MEM_ALLOC_CHK
  mexPrintf("temp_NZE_R=(NonZeroElem**)(%d)\n",i);
#endif
  NonZeroElem** temp_NZE_R=(NonZeroElem**)mxMalloc(i);
#ifdef MEM_ALLOC_CHK
  mexPrintf("temp_NZE_R=(NonZeroElem**)(%d)\n",i);
#endif
  NonZeroElem** temp_NZE_C=(NonZeroElem**)mxMalloc(i);
  memset(temp_NZE_R, 0, i);
  memset(temp_NZE_C, 0, i);
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
  memset(NbNZRow, 0, i);
  memset(NbNZCol, 0, i);
  i=periods*Size*sizeof(*b);
  memset(b,0,i);
#ifdef PRINT_OUT
  mexPrintf("Now looping\n");
  mexEvalString("drawnow;");
#endif
  for (t=0;t<periods;t++)
    {
#ifdef PRINT_OUT
      mexPrintf("t=%d\n",t);
#endif
      int ti_y_kmin=-min( t            , y_kmin);
      int ti_y_kmax= min( periods-(t+1), y_kmax);
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
          //mexEvalString("drawnow;");
#endif
          if (var<(periods+y_kmax)*Size)
            {
              lag=it4->first.second;
#ifdef PRINT_OUT
              mexPrintf(", lag =%d, ti_y_kmin=%d, ti_y_kmax=%d ", lag, ti_y_kmin, ti_y_kmax);
#endif
              if (lag<=ti_y_kmax && lag>=ti_y_kmin)
                {
                  //mexPrintf("u_index=%d, eq=%d, var=%d, lag=%d ",it4->second+u_count_init*t, eq, var, lag);
                  var+=Size*t;
                  //mexPrintf("u_index=%d, eq=%d, var=%d, lag=%d ",it4->second+u_count_init*t, eq, var, lag);
                  NbNZRow[eq]++;
                  NbNZCol[var]++;
#ifdef NEW_ALLOC
                  first=mem_mngr.mxMalloc_NZE();
#else
                  first=(NonZeroElem*)mxMalloc(sizeof(NonZeroElem));
#endif
                  first->NZE_C_N=NULL;
                  first->NZE_R_N=NULL;
                  first->u_index=it4->second+u_count_init*t;
                  first->r_index=eq;
                  first->c_index=var;
                  first->lag_index=lag;
                  /*if(eq==0 && var==0)
                    mexPrintf("alloc FNZE_R[0]=%x\n",first);*/
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
#ifdef PRINT_OUT
                  mexPrintf("=> ");
#endif
                }
              else
                {
#ifdef PRINT_OUT
                  mexPrintf("nn ");
                  mexPrintf("tmp_b+=u[%d]*y[index_var[%d]]\n",it4->second+u_count_init*t,var+Size*(y_kmin+t));
                  mexPrintf("tmp_b+=u[%d](%f)*y[%d(%d)](%f)",it4->second+u_count_init*t,u[it4->second+u_count_init*t], index_vara[var+Size*(y_kmin+t)],var+Size*(y_kmin+t),y[index_vara[var+Size*(y_kmin+t)]]);
                  mexEvalString("drawnow;");
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
#ifdef PRINT_OUT
              mexPrintf("=> b[%d]=%f\n", eq, u[b[eq]]);
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

void SparseMatrix::ShortInit(int periods, int y_kmin, int y_kmax, int Size, std::map<std::pair<std::pair<int, int> ,int>, int> IM)
{
  int t, eq, var, lag;
  double tmp_b=0.0;
  std::map<std::pair<std::pair<int, int> ,int>, int>::iterator it4;

  for (t=0;t<periods;t++)
    {
#ifdef PRINT_OUT
      mexPrintf("t=%d\n",t);
#endif
      int ti_y_kmin=-min( t        , y_kmin);
      int ti_y_kmax= min( periods-(t+1), y_kmax);
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
  long int i,j,/*nop=nop4/4*/nop=nop4/2, t, index_d, k;
  double r=0.0;
  bool OK=true;
  t_save_op_s *save_op_s, *save_opa_s, *save_opaa_s;
  int *diff1, *diff2;
#ifdef MEM_ALLOC_CHK
  mexPrintf("diff1=(int*)mxMalloc(%d*sizeof(int))\n",nop);
#endif
  diff1=(int*)mxMalloc(nop*sizeof(int));
#ifdef MEM_ALLOC_CHK
  mexPrintf("diff1=(int*)mxMalloc(%d*sizeof(int))\n",nop);
#endif
  diff2=(int*)mxMalloc(nop*sizeof(int));
#ifdef MEM_ALLOC_CHK
  mexPrintf("ok\n");
#endif
  j=k=i=0;
  while (i<nop4 && OK)
    {
      save_op_s=(t_save_op_s*)&(save_op[i]);
      save_opa_s=(t_save_op_s*)&(save_opa[i]);
      save_opaa_s=(t_save_op_s*)&(save_opaa[i]);
      diff1[j]=save_op_s->first-save_opa_s->first;
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
#ifdef PROFILER
  if(OK)
    mexPrintf("at %d same construction\n",beg_t);
  else
    mexPrintf("at %d different construction\n",beg_t);
#endif
  // the same pivot for all remaining periods
  if (OK)
    for (i=beg_t;i<periods;i++)
      for (j=0;j<Size;j++)
        pivot[i*Size+j]=pivot[(i-1)*Size+j]+Size;
  if (OK)
    {

#ifdef WRITE_u
      long int i_toto=0;
      fstream toto;
      toto.open("compare_s.txt", std::ios::out);
#endif

      t=1;
      for (;t<periods-beg_t-max(y_kmax,y_kmin);t++)
        {
          i=j=0;
          while (i<nop4)
            {
              save_op_s=(t_save_op_s*)(&(save_op[i]));
              index_d=save_op_s->first+t*diff1[j];
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
              switch (save_op_s->operat)
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
                    u[index_d]-=u[save_op_s->second+t*diff2[j]]*r;
#ifdef PRINT_u
                    mexPrintf("FSUB u[%d]-=u[%d](%f)*r(%f)=(%f)\n",index_d,save_op_s->second+t*diff2[j],u[save_op_s->second+t*diff2[j]],r,u[index_d] );
#endif
                    i+=3;
                    break;
                  case IFLESS:
                    u[index_d]=-u[save_op_s->second+t*diff2[j]]*r;
#ifdef PRINT_u
                    mexPrintf("FLESS u[%d]=-u[%d](%f)*r(%f)=(%f)\n",index_d,save_op_s->second+t*diff2[j],u[save_op_s->second+t*diff2[j]],r,u[index_d] );
#endif
                    i+=3;
                    break;
                }
              j++;
            }
        }
      int t1=t;
      for (t=t1;t<periods-beg_t;t++)
        {
          i=j=0;
          while (i<nop4)
            {
              save_op_s=(t_save_op_s*)(&(save_op[i]));
              if (save_op_s->lag<((periods-beg_t)-t))
                {
                  index_d=save_op_s->first+t*diff1[j];
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
                  switch (save_op_s->operat)
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
                        u[index_d]-=u[save_op_s->second+t*diff2[j]]*r;
#ifdef PRINT_u
                        mexPrintf("FSUB u[%d]-=u[%d](%f)*r(%f)=(%f)\n",index_d,save_op_s->second+t*diff2[j],u[save_op_s->second+t*diff2[j]],r,u[index_d] );
#endif
                        i+=3;
                        break;
                      case IFLESS:
                        u[index_d]=-u[save_op_s->second+t*diff2[j]]*r;
#ifdef PRINT_u
                        mexPrintf("FLESS u[%d]=-u[%d](%f)*r(%f)=(%f)\n",index_d,save_op_s->second+t*diff2[j],u[save_op_s->second+t*diff2[j]],r,u[index_d] );
#endif
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
#ifdef WRITE_u
      toto.close();
      filename+=" stopped";
      mexEvalString("st=fclose('all');clear all;");
      mexErrMsgTxt(filename.c_str());
#endif
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
  mexPrintf("begining of run_triangular nop_all=%d\n",nop_all);
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
          y[eq] += slowc_l*yy;
#ifdef PRINT_OUT_y1
          mexPrintf("=%f (%f)\n",double(yy),double(y[eq]));
#endif
        }
    }
  return res1;
}


void
SparseMatrix::close_swp_file()
{
  mem_mngr.close_swp_f();
}





void
SparseMatrix::Direct_Simulate(int blck, int y_size, int it_, int y_kmin, int y_kmax, int Size, int periods, bool print_it, int iter)
{
  int i, j, k, l, m=0, nop;
  double err;
#ifdef PRINT_OUT_y
  int m1;
#endif
  int period = it_ * y_size, s_middle_count_loop = 0 ;
  if (periods>0)
    period = y_decal * y_size;
  //pctimer_t t1 = pctimer();
  clock_t t1 = clock();
  double uu, yy;
  //char tmp_s[150];
  //mexPrintf("period=%d\n",period);
#ifdef PRINT_OUT
  for (j = 0;j < it_ -y_kmin;j++)
    {
      for (i = 0;i < u_size;i++)
        {
          mexPrintf("u[%d]=%f ", j*u_size + i, double(u[j*u_size + i]));
        }
      mexPrintf("\n");
    }
#endif
  if (/*nb_first_table_u*/print_it > 0)
    {
      first_count_loop = it_;
      //mexPrintf("nb_prologue_table_y=%d Size%d\n",nb_prologue_table_y,Size);
      /*s_middle_count_loop = it_ -y_kmin - middle_count_loop + 1;*/
      s_middle_count_loop = it_ - y_kmin - (nb_prologue_table_y  / Size);
#ifdef PRINT_OUT
      mexPrintf("nb_first_table_u=%d first_count_loop=%d s_middle_count_loop=%d\n",nb_first_table_u,first_count_loop, s_middle_count_loop);
      mexPrintf("y_kmin=%d y_kmax=%d \n",y_kmin,y_kmax);
      mexPrintf("middle_count_loop=%d\n",middle_count_loop);
      mexPrintf("it_=%d\n",it_);
#endif
//#ifdef PRINT_OUT
      mexPrintf("----------------------------------------------------------------------\n");
      mexPrintf("      Simulate     iteration° %d     \n",iter+1);
      mexPrintf("      max. error=%.10e       \n",double(max_res));
      mexPrintf("      sqr. error=%.10e       \n",double(res2));
      mexPrintf("      abs. error=%.10e       \n",double(res1));
      mexPrintf("----------------------------------------------------------------------\n");
//endif
#ifdef PRINT_OUT
      mexPrintf("first_count\n");
      mexPrintf("(first_count_loop=%d - y_kmin=%d)=%d\n",first_count_loop,y_kmin, first_count_loop-y_kmin);
#endif
    }
  nop = 0;

  for (j = 0 ;j < first_count_loop - y_kmin;j++)
    {
#ifdef PRINT_OUT_b
      mexPrintf("------------------------------------------------------\n");
      mexPrintf("j=%d\n",j);
#endif
      first_table_u = F_first_table_u->pNext;
      first_i_table_u = F_first_i_table_u->pNext;
      for (i = 0;i < nb_first_table_u;i++)
        {
          switch (first_table_u->type)
            {
              case 1:
                u[first_table_u->index + j*first_i_table_u->index] = u[first_table_u->op1 + j * first_i_table_u->op1] * u[first_table_u->op2 + j * first_i_table_u->op2];
#ifdef PRINT_OUT_b
                mexPrintf("u[%d]=u[%d]*u[%d]=%f\n", first_table_u->index + j*first_i_table_u->index , first_table_u->op1 + j*first_i_table_u->op1, first_table_u->op2 + j*first_i_table_u->op2, double(u[first_table_u->index + j*first_i_table_u->index]));
#endif
                break;
              case 2:
                u[first_table_u->index + j*first_i_table_u->index] += u[first_table_u->op1 + j * first_i_table_u->op1] * u[first_table_u->op2 + j * first_i_table_u->op2];
#ifdef PRINT_OUT_b
                mexPrintf("u[%d]+=u[%d]*u[%d]=%f\n" , first_table_u->index + j*first_i_table_u->index, first_table_u->op1 + j*first_i_table_u->op1, first_table_u->op2 + j*first_i_table_u->op2, double(u[first_table_u->index + j*first_i_table_u->index]));
#endif
                break;
              case 3:
                u[first_table_u->index + j*first_i_table_u->index] = 1 / ( -u[first_table_u->op1 + j * first_i_table_u->op1]);
#ifdef PRINT_OUT_b
                mexPrintf("u[%d]=1/(-u[%d])=%f\n", first_table_u->index + j*first_i_table_u->index, first_table_u->op1 + j*first_i_table_u->op1, double(u[first_table_u->index + j*first_i_table_u->index]));
#endif
                break;
              case 5:
                Stack.push(u[first_table_u->index + j*first_i_table_u->index]);
#ifdef PRINT_OUT_b
                mexPrintf("push(u[%d])\n", first_table_u->index + j*first_i_table_u->index);
#endif
                break;
              case 6:
                u[first_table_u->index + j*first_i_table_u->index] = 1 / (1 - u[first_table_u->op1 + j * first_i_table_u->op1] * u[first_table_u->op2 + j * first_i_table_u->op2]);
#ifdef PRINT_OUT_b
                mexPrintf("u[%d]=1/(1-u[%d]*u[%d])=%f\n", first_table_u->index + j*first_i_table_u->index, first_table_u->op1 + j*first_i_table_u->op1, first_table_u->op2 + j*first_i_table_u->op2, double(u[first_table_u->index + j*first_i_table_u->index]));
#endif
                break;
              case 7:
                u[first_table_u->index + j*first_i_table_u->index] *= u[first_table_u->op1 + j * first_i_table_u->op1];
#ifdef PRINT_OUT_b
                mexPrintf("u[%d]*=u[%d]=%f\n", first_table_u->index + j*first_i_table_u->index, first_table_u->op1 + j*first_i_table_u->op1, double(u[first_table_u->index + j*first_i_table_u->index]));
#endif
                break;
            }
          if (isnan(u[first_table_u->index+ j*first_i_table_u->index]) || isinf(u[first_table_u->index+ j*first_i_table_u->index]))
            {
              mexPrintf("Error during the computation of u[%d] at time %d (in first_table_u) (operation type %d)",first_table_u->index,j,int(first_table_u->type));
              mexEvalString("st=fclose('all');clear all;");
              filename+=" stopped";
              mexErrMsgTxt(filename.c_str());
            }
          else if (fabs(u[first_table_u->index+ j*first_i_table_u->index])>very_big)
            {
              mexPrintf("(first) big u[%d]=%f>%f in type=%d",first_table_u->index+ j*first_i_table_u->index, double(u[first_table_u->index+ j*first_i_table_u->index]),very_big,first_table_u->type);
              mexEvalString("st=fclose('all');clear all;");
              filename+=" stopped";
              mexErrMsgTxt(filename.c_str());
            }
          first_table_u = first_table_u->pNext;
          first_i_table_u = first_i_table_u->pNext;
          nop++;
        }
    }
#ifdef PRINT_OUT_b
  mexPrintf("//////////////////////////////////////////////////////////////\n");
  mexPrintf("prologue\n");
#endif
  //int nb_prologue_push=0;
  prologue_table_u = F_prologue_table_u->pNext;
  for (i = 0;i < nb_prologue_table_u ;i++)
    {
      switch (prologue_table_u->type)
        {
          case 1:
            u[prologue_table_u->index ] = u[prologue_table_u->op1 ] * u[prologue_table_u->op2 ];
#ifdef PRINT_OUT_b
            mexPrintf("u[%d]=u[%d]*u[%d]=%f\n", prologue_table_u->index , prologue_table_u->op1 , prologue_table_u->op2 , double(u[prologue_table_u->index ]));
#endif
            break;
          case 2:
            u[prologue_table_u->index ] += u[prologue_table_u->op1 ] * u[prologue_table_u->op2 ];
#ifdef PRINT_OUT_b
            mexPrintf("u[%d]+=u[%d]*u[%d]=%f\n" , prologue_table_u->index , prologue_table_u->op1 , prologue_table_u->op2 , double(u[prologue_table_u->index ]));
#endif
            break;
          case 3:
            u[prologue_table_u->index ] = 1 / ( -u[prologue_table_u->op1 ]);
#ifdef PRINT_OUT_b
            mexPrintf("u[%d]=1/(-u[%d])=%f\n", prologue_table_u->index, prologue_table_u->op1, double(u[prologue_table_u->index]));
#endif
            break;
          case 5:
            //nb_prologue_push++;
            Stack.push(u[prologue_table_u->index]);
#ifdef PRINT_OUT_b
            mexPrintf("push(u[%d])\n", prologue_table_u->index );
#endif
            break;
          case 6:
            u[prologue_table_u->index ] = 1 / (1 - u[prologue_table_u->op1] * u[prologue_table_u->op2]);
#ifdef PRINT_OUT_b
            mexPrintf("u[%d]=1/(1-u[%d]*u[%d])=%f\n", prologue_table_u->index, prologue_table_u->op1, prologue_table_u->op2, double(u[prologue_table_u->index]));
#endif
            break;
          case 7:
            u[prologue_table_u->index] *= u[prologue_table_u->op1];
#ifdef PRINT_OUT_b
            mexPrintf("u[%d]*=u[%d]=%f\n", prologue_table_u->index, prologue_table_u->op1, double(u[prologue_table_u->index]));
#endif
            break;
        }
      if (isnan(u[prologue_table_u->index]) || isinf(u[prologue_table_u->index]))
        {
          mexPrintf("Error during the computation of u[%d] type=%d (in prologue_table_u)",prologue_table_u->index, prologue_table_u->type);
          mexEvalString("st=fclose('all');clear all;");
          filename+=" stopped";
          mexErrMsgTxt(filename.c_str());
        }
      else if (fabs(u[prologue_table_u->index])>very_big)
        {
          mexPrintf("(prologue) big u[%d]=%f>%f type=%d",prologue_table_u->index, double(u[prologue_table_u->index]), very_big, prologue_table_u->type);
          mexEvalString("st=fclose('all');clear all;");
          filename+=" stopped";
          mexErrMsgTxt(filename.c_str());
        }
      prologue_table_u = prologue_table_u->pNext;
      nop++;
    }
#ifdef PRINT_OUT_b
  mexPrintf("middle_u (s_middle_count_loop=%d\n", s_middle_count_loop);
#endif
  //int nb_middle_push=0;
  for (j = 0;j < s_middle_count_loop /*- y_kmin*/;j++)
    {
      //cout << "j=" << j << "\n";
#ifdef PRINT_OUT_b
      mexPrintf("-----------------------------------------------------------------\n");
#endif
      middle_table_u = F_middle_table_u->pNext;
      middle_i_table_u = F_middle_i_table_u->pNext;
      for (i = 0;i < nb_middle_table_u;i++)
        {
          switch (middle_table_u->type)
            {
              case 1:
                u[middle_table_u->index + j*middle_i_table_u->index] = u[middle_table_u->op1 + j * middle_i_table_u->op1] * u[middle_table_u->op2 + j * middle_i_table_u->op2];
#ifdef PRINT_OUT_b
                mexPrintf("u[%d+%d*%d=%d]=u[%d]*u[%d]=%f\n", middle_table_u->index, j, middle_i_table_u->index, middle_table_u->index + j*middle_i_table_u->index, middle_table_u->op1 + j*middle_i_table_u->op1, middle_table_u->op2 + j*middle_i_table_u->op2, double(u[middle_table_u->index + j*middle_i_table_u->index]));
#endif
                break;
              case 2:
                u[middle_table_u->index + j*middle_i_table_u->index] += u[middle_table_u->op1 + j * middle_i_table_u->op1] * u[middle_table_u->op2 + j * middle_i_table_u->op2];
#ifdef PRINT_OUT_b
                mexPrintf("u[%d+%d*%d=%d]+=u[%d]*u[%d]=%f\n" , middle_table_u->index, j, middle_i_table_u->index , middle_table_u->index + j*middle_i_table_u->index , middle_table_u->op1 + j*middle_i_table_u->op1, middle_table_u->op2 + j*middle_i_table_u->op2, double(u[middle_table_u->index + j*middle_i_table_u->index]));
#endif
                break;
              case 3:
                u[middle_table_u->index + middle_i_table_u->index] = 1 / ( -u[middle_table_u->op1 + j * middle_i_table_u->op1]);
#ifdef PRINT_OUT_b
                mexPrintf("u[%d+%d*%d=%d]=1/(-u[%d])=%f\n", middle_table_u->index, j, middle_i_table_u->index, middle_table_u->index + j*middle_i_table_u->index, middle_table_u->op1 + j*middle_i_table_u->op1, double(u[middle_table_u->index + j*middle_i_table_u->index]));
#endif
                break;
              case 5:
#ifdef PRINT_OUT_b
                mexPrintf("push(u[%d+%d*%d=%d])\n", middle_table_u->index, j, middle_i_table_u->index, middle_table_u->index + j*middle_i_table_u->index);
#endif
                //nb_middle_push++;
                Stack.push(u[middle_table_u->index + j*middle_i_table_u->index]);
                break;
              case 6:
                u[middle_table_u->index + j*middle_i_table_u->index] = 1 / (1 - u[middle_table_u->op1 + j * middle_i_table_u->op1] * u[middle_table_u->op2 + j * middle_i_table_u->op2]);
#ifdef PRINT_OUT_b
                mexPrintf("u[%d+%d*%d=%d]=1/(1-u[%d]*u[%d])=%f\n", middle_table_u->index, j, middle_i_table_u->index, middle_table_u->index + j*middle_i_table_u->index, middle_table_u->op1 + j*middle_i_table_u->op1, middle_table_u->op2 + j*middle_i_table_u->op2, double(u[middle_table_u->index + j*middle_i_table_u->index]));
#endif
                break;
              case 7:
                u[middle_table_u->index + j*middle_i_table_u->index] *= u[middle_table_u->op1 + j * middle_i_table_u->op1];
#ifdef PRINT_OUT_b
                mexPrintf("u[%d+%d*%d=%d]*=u[%d]=%f\n", middle_table_u->index, j, middle_i_table_u->index, middle_table_u->index + j*middle_i_table_u->index, middle_table_u->op1 + j*middle_i_table_u->op1, double(u[middle_table_u->index + j*middle_i_table_u->index]));
#endif
                break;
            }
          if (isnan(u[middle_table_u->index+ j*middle_i_table_u->index]) || isinf(u[middle_table_u->index+ j*middle_i_table_u->index]))
            {
              mexPrintf("Error during the computation of u[%d] at time %d type=%d (in middle_table_u)",middle_table_u->index,j,middle_table_u->type);
              mexEvalString("st=fclose('all');clear all;");
              filename+=" stopped";
              mexErrMsgTxt(filename.c_str());
            }
          else if (fabs(u[middle_table_u->index+ j*middle_i_table_u->index])>very_big)
            {
              mexPrintf("(middle) big u[%d]=%f>%f type=%d",middle_table_u->index+ j*middle_i_table_u->index, double(u[middle_table_u->index+ j*middle_i_table_u->index]), very_big,middle_table_u->type);
              mexEvalString("st=fclose('all');clear all;");
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
  for (i = 0;i < nb_last_table_u ;i++)
    {
      switch (last_table_u->type)
        {
          case 1:
            u[last_table_u->index] = u[last_table_u->op1] * u[last_table_u->op2];
#ifdef PRINT_OUT_b
            mexPrintf("u[%d]=u[%d]*u[%d]=%f\n", last_table_u->index, last_table_u->op1, last_table_u->op2, double(u[last_table_u->index]));
#endif
            break;
          case 2:
            u[last_table_u->index] += u[last_table_u->op1] * u[last_table_u->op2];
#ifdef PRINT_OUT_b
            mexPrintf("u[%d]+=u[%d]*u[%d]=%f\n", last_table_u->index, last_table_u->op1, last_table_u->op2, double(u[last_table_u->index]));
#endif
            break;
          case 3:
            u[last_table_u->index] = 1 / ( -u[last_table_u->op1]);
#ifdef PRINT_OUT_b
            mexPrintf("u[%d]=1/(-u[%d])=%f\n", last_table_u->index, last_table_u->op1, double(u[last_table_u->index]));
#endif
            break;
          case 5:
            Stack.push(u[last_table_u->index]);
#ifdef PRINT_OUT_b
            mexPrintf("push(u[%d])\n", last_table_u->index);
#endif
            break;
          case 6:
            u[last_table_u->index] = 1 / (1 - u[last_table_u->op1] * u[last_table_u->op2]);
#ifdef PRINT_OUT_b
            mexPrintf("u[%d]=1/(1-u[%d]*u[%d])=%f\n", last_table_u->index, last_table_u->op1, last_table_u->op2, double(u[last_table_u->index]));
#endif
            break;
          case 7:
            u[last_table_u->index] *= u[last_table_u->op1];
#ifdef PRINT_OUT_b
            mexPrintf("u[%d]*=u[%d]=%f\n", last_table_u->index, last_table_u->op1, double(u[last_table_u->index]));
#endif
            break;
        }
      if (isnan(u[last_table_u->index]) || isinf(u[last_table_u->index]))
        {
          mexPrintf("Error during the computation of u[%d] type=%d (in last_table_u)",last_table_u->index,last_table_u->type);
          mexEvalString("st=fclose('all');clear all;");
          filename+=" stopped";
          mexErrMsgTxt(filename.c_str());
        }
      else if (fabs(u[last_table_u->index])>very_big)
        {
          mexPrintf("(last) big u[%d]=%f>%f type=%d",last_table_u->index, double(u[last_table_u->index]),very_big,last_table_u->type);
          mexEvalString("st=fclose('all');clear all;");
          filename+=" stopped";
          mexErrMsgTxt(filename.c_str());
        }
      last_table_u = last_table_u->pNext;
      nop++;
    }
  res1 = res2 =max_res = 0;
  for (i = nb_last_table_y - 1;i >= 0;i--)
    {
      k = last_table_y[i].index;
      yy = 0;
      /*y[period + k] = 0;*/
#ifdef PRINT_OUT_y
      //mexPrintf("it_=%d\n", it_);
      mexPrintf("->y[it_*y_size+%d]=y[%d]=", k, period + k);
#endif
      for (j = last_table_y[i].nb - 1;j >= 0;j--)
        {
          uu = Stack.top();
          Stack.pop();
          m = last_table_y[i].y_index[j];
#ifdef PRINT_OUT_y
          if (j > 0)
            {
              if (m >= 0)
                mexPrintf("u[%d](%f)*y[%d](%f)+", last_table_y[i].u_index[j], double(uu), m + period, double(y[period + m]));
              else
                mexPrintf("u[%d](%f)+", last_table_y[i].u_index[j], double(uu));
            }
          else
            {
              if (m >= 0)
                mexPrintf("u[%d](%f)*y[%d](%f)", last_table_y[i].u_index[j], double(uu), m + period, double(y[period + m]));
              else
                mexPrintf("u[%d](%f)", last_table_y[i].u_index[j], double(uu));
            }
#endif
          if (m >= 0)
            yy/*y[period + k]*/ += uu * y[period + m];
          else
            yy/*y[period + k]*/ += uu;
        }
      if (isnan(yy) || isinf(yy))
        {
          mexPrintf("Error during the computation of y[%d] at time %d (in last_table_u)",k,period);
          mexEvalString("st=fclose('all');clear all;");
          filename+=" stopped";
          mexErrMsgTxt(filename.c_str());
        }
      /*if(((k-73) % y_size)==0)
        mexPrintf("y[it_*y_size +73]=%f \n",yy);*/
      err = fabs(yy - y[period + k]);
      res1 += err;
      if (max_res<err)
        max_res=err;
      res2 += err*err;
      y[period + k] += slowc*(yy - y[period + k]);
#ifdef PRINT_OUT_y
      mexPrintf("=%f\n", double(y[period + k]));
#endif
      nop++;
    }
#ifdef PRINT_OUT_y
  int deca=y_kmin - (nb_prologue_table_y  / Size);
#endif
  for (j = s_middle_count_loop+y_decal ;j >y_decal;j--)
    {
#ifdef PRINT_OUT_y
      mexPrintf("(per)j=%d y_decal=%d deca=%d in y compute\n",j,y_decal,deca);
#endif
      for (i = nb_middle_table_y - 1;i >= 0;i--)
        {
          k = middle_table_y[i].index + (j-1) * middle_i_table_y[i].index;
          yy = 0;
#ifdef PRINT_OUT_y
          mexPrintf("(0)y[%d]=", k );
#endif
          for (l = middle_table_y[i].nb - 1;l >= 0;l--)
            {
              uu = Stack.top();
              Stack.pop();
              m = middle_table_y[i].y_index[l] + (j/*-deca*/-1) * middle_i_table_y[i].y_index[l];
#ifdef PRINT_OUT_y
              if (/*m*/middle_table_y[i].y_index[l] >= 0)
                {
                  m1 = middle_table_y[i].u_index[l] + (j-deca-1) * middle_i_table_y[i].u_index[l];
                  if (l > 0)
                    mexPrintf("u[%d](%f)*y[%d](%f)+", m1, double(uu), m, double(y[m]));
                  else
                    mexPrintf("u[%d](%f)*y[%d](%f)", m1, double(uu), m, double(y[m]));
                }
              else
                {
                  m1 = middle_table_y[i].u_index[l] + (j-deca-1) * middle_i_table_y[i].u_index[l];
                  if (l > 0)
                    mexPrintf("u[%d](%f)*y[%d](%f)+", m1, double(uu), m, 1.0);
                  else
                    mexPrintf("u[%d](%f)*y[%d](%f)", m1, double(uu), m, 1.0);
                }
#endif
              if (/*m*/middle_table_y[i].y_index[l] >= 0)
                yy += uu * y[m];
              else
                yy += uu;
            }
          //mexPrintf("y[%d]=%f\n",k,yy);
          if (isnan(yy) || isinf(yy))
            {
              mexPrintf("Error during the computation of y[%d] at time %d (in middle_table_u)",middle_table_y[i].index % y_size,j+middle_table_y[i].index / y_size);
              mexEvalString("st=fclose('all');clear all;");
              filename+=" stopped";
              mexErrMsgTxt(filename.c_str());
            }
          /*if(((k-73) % y_size)==0)
            mexPrintf("y[it_*y_size +73]=%f \n",yy);*/
          err = fabs(yy - y[k]);
          res1 += err;
          if (max_res<err)
            max_res=err;
          res2 += err*err;
          y[k] += slowc*(yy - y[k]);
#ifdef PRINT_OUT_y
          mexPrintf("=%f\n", double(y[k]));
#endif
          nop++;
        }
    }
  for (i = nb_prologue_table_y - 1;i >= 0;i--)
    {
      k = prologue_table_y[i].index;
      yy = 0;
#ifdef PRINT_OUT_y
      mexPrintf("(1)y[%d]=", k+y_decal*y_size);
#endif
      for (j = prologue_table_y[i].nb - 1;j >= 0;j--)
        {
          //nb_prologue_pop++;
          uu = Stack.top();
          Stack.pop();
#ifdef PRINT_OUT_y
          if (prologue_table_y[i].y_index[j] >= 0)
            {
              if (j > 0)
                mexPrintf("u[%d](%f)*y[%d](%f)+", prologue_table_y[i].u_index[j], double(uu), prologue_table_y[i].y_index[j], double(y[prologue_table_y[i].y_index[j]]));
              else
                mexPrintf("u[%d](%f)*y[%d](%f)", prologue_table_y[i].u_index[j], double(uu), prologue_table_y[i].y_index[j], double(y[prologue_table_y[i].y_index[j]]));
            }
          else
            {
              if (j > 0)
                mexPrintf("u[%d](%f)*y[%d](%f)+", prologue_table_y[i].u_index[j], double(uu), prologue_table_y[i].y_index[j], 1.0);
              else
                mexPrintf("u[%d](%f)*y[%d](%f)", prologue_table_y[i].u_index[j], double(uu), prologue_table_y[i].y_index[j], 1.0);
            }
#endif
          if (prologue_table_y[i].y_index[j] >= 0)
            yy += uu * y[prologue_table_y[i].y_index[j]+y_decal*y_size];
          else
            yy += uu;
        }
      if (isnan(yy) || isinf(yy))
        {
          mexPrintf("Error during the computation of y[%d] at time %d (in prologue_table_u)",k % y_size, k / y_size);
          mexEvalString("st=fclose('all');clear all;");
          filename+=" stopped";
          mexErrMsgTxt(filename.c_str());
        }
      err = fabs(yy - y[k+y_decal*y_size]);
      res1 += err;
      if (max_res<err)
        max_res=err;
      res2 += err*err;
      y[k+y_decal*y_size] += slowc*(yy - y[k+y_decal*y_size]);
#ifdef PRINT_OUT_y
      mexPrintf("=%f\n", double(y[k+y_decal*y_size]));
#endif
      nop++;
    }
  //mexPrintf("nb_prologue_push=%d et nb_prologue_pop=%d\n",nb_prologue_push,nb_prologue_pop);
  for (i = nb_first_table_y - 1;i >= 0;i--)
    {
      k = first_table_y[i].index;
      yy = 0;
#ifdef PRINT_OUT_y
      mexPrintf("(2)y[%d]=", k);
#endif
      for (j = first_table_y[i].nb - 1;j >= 0;j--)
        {
          uu = Stack.top();
          Stack.pop();
#ifdef PRINT_OUT_y
          if (j > 0)
            mexPrintf("u[%d](%f)*y[%d](%f)+", first_table_y[i].u_index[j], double(uu), first_table_y[i].y_index[j], double(y[first_table_y[i].y_index[j]]));
          else
            mexPrintf("u[%d](%f)*y[%d](%f)", first_table_y[i].u_index[j], double(uu), first_table_y[i].y_index[j], double(y[first_table_y[i].y_index[j]]));
#endif
          if (m >= 0)
            yy += uu * y[first_table_y[i].y_index[j]];
          else
            yy += uu;
        }
      if (isnan(yy) || isinf(yy))
        {
          mexPrintf("Error during the computation of y[%d] (in first_table_u)",k);
          mexEvalString("st=fclose('all');clear all;");
          filename+=" stopped";
          mexErrMsgTxt(filename.c_str());
        }
      err = fabs(yy - y[k]);
      res1 += err;
      if (max_res<err)
        max_res=err;
      res2 += err*err;
      y[k] += slowc*(yy - y[k]);
#ifdef PRINT_OUT_y
      mexPrintf("=%f\n", double(y[k]));
#endif
      nop++;
    }
  //pctimer_t t2 = pctimer();
  clock_t t2 = clock();
  if (!Stack.empty())
    {
      mexPrintf("Error Stack not empty (%d)",Stack.empty());
      mexEvalString("st=fclose('all');clear all;");
      filename+=" stopped";
      mexErrMsgTxt(filename.c_str());
    }
  if (/*nb_first_table_u > 0*/print_it)
    {
      mexPrintf("res1    = %.10e\n", double(res1));
      mexPrintf("res2    = %.10e\n", double(res2));
      mexPrintf("max_res = %.10e\n", double(max_res));
      mexPrintf("(**%f milliseconds u_count : %d  nop : %d **)\n", 1000.0*(double(t2) - double(t1))/double(CLOCKS_PER_SEC), u_count, nop);
      mexEvalString("drawnow;");
    }
}











int
SparseMatrix::simulate_NG1(int blck, int y_size, int it_, int y_kmin, int y_kmax, int Size, int periods, bool print_it, bool cvg, int &iter)
{
  /*Triangularisation at each period of a block using a simple gaussian Elimination*/
  t_save_op_s *save_op_s;
  int start_compare=y_kmin;
  bool record=false;
  int *save_op=NULL, *save_opa=NULL, *save_opaa=NULL;
  long int nop=0, nopa=0;
  int tbreak=0, last_period=periods;
  int i,j,k;
  int pivj=0, pivk=0;
  NonZeroElem *first, *firsta, *first_sub, *first_piv, *first_suba;
  double piv_abs, first_elem;
  //SparseMatrix sparse_matrix;
  //mexPrintf("->u_count=%d &u_count=%x\n",u_count,&u_count);
#ifdef RECORD_ALL
  int save_u_count=u_count;
#endif
  u_count_alloc_save=u_count_alloc;
#ifdef PROFILER
  long int ndiv=0, nsub=0, ncomp=0, nmul=0;
  double tinsert=0, tdelete=0, tpivot=0, tbigloop=0;
  clock_t td1;
  int nbpivot=0, nbdiv=0, nbless=0, nbpivot_it=0, nbRealloc=0;
#endif
  //pctimer_t t01;
  clock_t t01;
  //pctimer_t t1 = pctimer();
  clock_t t1 = clock();
#ifdef PROFILER
  tdelete1=0; tdelete2=0; tdelete21=0; tdelete22=0; tdelete221=0; tdelete222=0;
#endif
  if (iter>0)
    mexPrintf("Sim : %f ms\n",(1000.0*(double(clock())-double(time00)))/double(CLOCKS_PER_SEC));
#ifdef MEMORY_LEAKS
  mexEvalString("feature('memstats');");
#endif
  if (isnan(res1) || isinf(res1))
    {
      if (slowc_save<1e-8)
        {
          mexPrintf("Dynare cannot improve the simulation\n");
          mexEvalString("drawnow;");
          mexEvalString("st=fclose('all');clear all;");
          filename+=" stopped";
          mexErrMsgTxt(filename.c_str());
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
  if (((res1/res1a-1)>-0.3) && symbolic)
    {
      if (start_compare==y_kmin)
        {
          mexPrintf("Divergence or slowdown occured during simulation.\nIn the next iteration, pivoting method will be applied for a longer period.\n");
          start_compare=min(tbreak_g,periods);
        }
      else
        {
          mexPrintf("Divergence or slowdown occured during simulation.\nIn the next iteration, pivoting method will be applied to all periods.\n");
          symbolic=false;
          alt_symbolic=true;
          markowitz_c_s=markowitz_c;
          markowitz_c=0;
        }
    }
  res1a=res1;
  mexPrintf("-----------------------------------\n");
  mexPrintf("      Simulate     iteration° %d     \n",iter+1);
  mexPrintf("      max. error=%.10e       \n",double(max_res));
  mexPrintf("      sqr. error=%.10e       \n",double(res2));
  mexPrintf("      abs. error=%.10e       \n",double(res1));
  mexPrintf("-----------------------------------\n");
  if (cvg)
    return(0);
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
#ifdef PRINT_OUT
      mexPrintf("ShortInit\n");
      mexEvalString("drawnow;");
#endif
      ShortInit(periods, y_kmin, y_kmax, Size, IM_i);
#ifdef PRINT_OUT
      mexPrintf("run_triangular\n");
      mexEvalString("drawnow;");
#endif
      run_triangular(nop_all,save_op_all);
#ifdef PRINT_OUT
      mexPrintf("OK\n");
      mexEvalString("drawnow;");
#endif
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
#ifdef PRINT_OUT
        mexPrintf("done\n");
        mexEvalString("drawnow;");
#endif
        for (int t=0;t<periods;t++)
          {
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
                if /*(t<=y_kmin) */((symbolic && t<=start_compare) || !symbolic)
                  {
#ifdef MARKOVITZ
                    double piv_v[Size];
                    int pivj_v[Size], pivk_v[Size], NR[Size], l=0, N_max=0;
                    bool one=false;
                    piv_abs=0;
#endif
                    for (j=0;j<nb_eq;j++)
                      {
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

                    pivot[i]=pivj;
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
                if (piv_abs<eps)
                  {
                    mexPrintf("Error: singular system\n");
                    mexEvalString("drawnow;");
                    mexEvalString("st=fclose('all');clear all;");
                    filename+=" stopped";
                    mexErrMsgTxt(filename.c_str());
                  }
                /*divide all the non zeros elements of the line pivj by the max_pivot*/
#ifdef PRINT_OUT
                mexPrintf("ok2\n");
                mexEvalString("drawnow;");
#endif
                int nb_var=At_Row(pivj,&first);
#ifdef PRINT_OUT
                mexPrintf("nb_var=%d\n",nb_var);
#endif
                for (j=0;j<nb_var;j++)
                  {
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
                            save_op_s->first=first->u_index;
                            save_op_s->lag=first->lag_index;
                          }
                        nop+=2;
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
                    first=first->NZE_R_N;
                  }
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
#ifdef PRINT_OUT
                mexPrintf("ok3\n");
                mexEvalString("drawnow;");
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
#ifdef PRINT_OUT
                if(iter>0)
                  {
                mexPrintf("ok4 nb_eq=%d iter=%d\n",nb_eq,iter);
                mexEvalString("drawnow;");
                  }
#endif
                for (j=0;j<nb_eq;j++)
                  {
                    int row=first->r_index;
#ifdef PRINT_OUT
                    mexPrintf("line_done[%d]=%d\n", row, int(line_done[row]));
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
                        int nb_var_piv=nb_var_piva;
                        first_piv=first_piva;
                        int nb_var_sub=At_Row(row,&first_sub);
                        int l_sub=0, l_piv=0;
                        int sub_c_index=first_sub->c_index, piv_c_index=first_piv->c_index;
                        int tmp_lag=first_sub->lag_index;
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
                            if (l_sub<nb_var_sub && (sub_c_index<piv_c_index || l_piv>=nb_var_piv))
                              {
#ifdef PRINT_OUT
                                if(iter>0)
                                  {
                                    mexPrintf("  u[%d]=u[%d]\n",first_sub->u_index,first_sub->u_index);
                                    mexEvalString("drawnow;");
                                  }
#endif
                                first_sub=first_sub->NZE_R_N;
                                if (first_sub)
                                  sub_c_index=first_sub->c_index;
                                else
                                  sub_c_index=Size*periods;
                                l_sub++;
                              }
                            else if (sub_c_index>piv_c_index || l_sub>=nb_var_sub)
                              {
                                int tmp_u_count=Get_u();
#ifdef PROFILER
                                clock_t td0=clock();
#endif
                                int lag=first_piv->c_index/Size-row/Size;
                                Insert(row,first_piv->c_index,tmp_u_count,lag);
#ifdef PROFILER
                                tinsert+=clock()-td0;
#endif
                                u[tmp_u_count]=-u[first_piv->u_index]*first_elem;
#ifdef PROFILER
                                nbless++;
#endif
#ifdef WRITE_u
                                if ((periods-t)<=y_kmax)
                                  {
                                    toto << i_toto << " u[" /*<< tmp_u_count*/ << "]=-u[" /*<< first_piv->u_index*/ << "]*" << first_elem << "=" << u[tmp_u_count] << endl;
                                    i_toto++;
                                  }
#endif
#ifdef PRINT_u
                                if(iter>0)
                                  {
                                mexPrintf("FLESS u[%d]=-u[%d](%f)*r(%f)=(%f)   |",tmp_u_count,first_piv->u_index,u[first_piv->u_index],first_elem,u[tmp_u_count]);
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
#ifdef PRINT_OUT
                                if(iter>0)
                                  {
                                mexPrintf("  u[%d at (%d, %d) lag %d]=-u[%d]*%f\n",tmp_u_count,row , first_piv->c_index, first_piv->c_index/Size-row/Size, first_piv->u_index ,double(first_elem));
                                    mexEvalString("drawnow;");
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
#ifdef PRINT_OUT
                                    /*if(iter>0)
                                      {
                                    mexPrintf("   delete element [%d, %d] lag %d u[%d]\n",first_sub->r_index,first_sub->c_index,first_sub->lag_index,first_sub->u_index);
                                        mexEvalString("drawnow;");
                                      }*/
#endif
                                    firsta=first;
                                    first_suba=first_sub->NZE_R_N;
                                    /*if(iter>0)
                                      {
                                        mexPrintf("just after\n");
                                        mexEvalString("drawnow;");
                                      }*/
#ifdef PROFILER
                                    clock_t td0=clock();
#endif
                                    Delete(first_sub->r_index,first_sub->c_index, Size, b);
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
                        /*if(iter>0)
                          {
                            mexPrintf("ok4f j=%d\n",j);
                            mexEvalString("drawnow;");
                          }*/

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
                    else
                      first=first->NZE_C_N;

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
                            mexPrintf("t=%d over periods=%d t0=%f t1=%f y_kmin=%d y_kmax=%d\n",t,periods,1000*(ta-t00), 1000.0*(double(clock())-double(ta))/double(CLOCKS_PER_SEC),y_kmin, y_kmax);
                            mexPrintf("compare time %f ms\n",1000.0*(double(clock())-double(ta))/double(CLOCKS_PER_SEC));
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
  mexPrintf("nbpivot=%d, nbdiv=%d, nbless=%d, nop=%d nbpivot_it=%d nbRealloc=%d\n", nbpivot, nbdiv, nbless, nbpivot + nbdiv + nbless, nbpivot_it, nbRealloc);
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
#endif
#ifdef RECORD_ALL
  if (!((record_all && nop_all)||g_nop_all>0))
    {
      u_count=save_u_count;
      End(Size);
    }
#endif
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
  return(0);
}

void
SparseMatrix::fixe_u(double **u, int u_count_int, int max_lag_plus_max_lead_plus_1)
{
  u_count=u_count_int * periods;
  u_count_alloc = 2*u_count;
  (*u)=(double*)mxMalloc(u_count_alloc*sizeof(double));
  memset((*u), 0, u_count_alloc*sizeof(double));
  u_count_init=max_lag_plus_max_lead_plus_1;
}


void
SparseMatrix::read_file_table_u(t_table_u **table_u, t_table_u **F_table_u, t_table_u **i_table_u, t_table_u **F_i_table_u, int *nb_table_u, bool i_to_do, bool shifting, int *nb_add_u_count, int y_kmin, int y_kmax, int u_size)
{
  char type;
  int i;
  i=SaveCode.tellp();
#ifdef PRINT_OUT
  mexPrintf("SaveCode.tellp()=%d\n",i);
#endif
  SaveCode.read(reinterpret_cast<char *>(nb_table_u), sizeof(*nb_table_u));
#ifdef PRINT_OUT
  mexPrintf("->*nb_table_u=%d\n", *nb_table_u);
#endif
  *table_u = (t_table_u*)mxMalloc(sizeof(t_table_u) - 2 * sizeof(int));
  *F_table_u = *table_u;
  if (i_to_do)
    {
#ifdef PRINT_OUT
      mexPrintf("=>i_table\n");
#endif
      (*i_table_u) = (t_table_u*)mxMalloc(sizeof(t_table_u) - 2 * sizeof(int));
      (*F_i_table_u) = (*i_table_u);
    }
  for (i = 0;i < *nb_table_u;i++)
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
            if ((*table_u)->index>max_u)
              max_u=(*table_u)->index;
            if ((*table_u)->index<min_u)
              min_u=(*table_u)->index;
            SaveCode.read(reinterpret_cast<char *>(&(*table_u)->op1), sizeof((*table_u)->op1));
            if (shifting)
              {
                (*table_u)->index -= y_kmin * u_size;
                if ((*table_u)->index>max_u)
                  max_u=(*table_u)->index;
                if ((*table_u)->index<min_u)
                  min_u=(*table_u)->index;
                (*table_u)->op1 -= y_kmin * u_size;
              }
#ifdef PRINT_OUT

            if ((*table_u)->type == 3)
              mexPrintf("u[%d]=-1/u[%d]\n", (*table_u)->index, (*table_u)->op1);
            else
              mexPrintf("u[%d]*=u[%d]\n", (*table_u)->index, (*table_u)->op1);
#endif
            if (i_to_do)
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
                if ((*i_table_u)->type == 3)
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
            if ((*table_u)->index>max_u)
              max_u=(*table_u)->index;
            if ((*table_u)->index<min_u)
              min_u=(*table_u)->index;
            SaveCode.read(reinterpret_cast<char *>(&(*table_u)->op1), sizeof((*table_u)->op1));
            SaveCode.read(reinterpret_cast<char *>(&(*table_u)->op2), sizeof((*table_u)->op2));
            if (shifting)
              {
                (*table_u)->index -= y_kmin * u_size;
                if ((*table_u)->index>max_u)
                  max_u=(*table_u)->index;
                if ((*table_u)->index<min_u)
                  min_u=(*table_u)->index;
                (*table_u)->op1 -= y_kmin * u_size;
                (*table_u)->op2 -= y_kmin * u_size;
              }
            if ((*table_u)->type == 1)
              {
#ifdef PRINT_OUT
                mexPrintf("u[%d]=u[%d]*u[%d]\n", (*table_u)->index, (*table_u)->op1, (*table_u)->op2);
#endif
                if (i_to_do)
                  (*nb_add_u_count)++;
              }
#ifdef PRINT_OUT
            else if ((*table_u)->type == 2)
              mexPrintf("u[%d]+=u[%d]*u[%d]\n", (*table_u)->index, (*table_u)->op1, (*table_u)->op2);
            else
              mexPrintf("u[%d]=1/(1-u[%d]*u[%d])\n", (*table_u)->index, (*table_u)->op1, (*table_u)->op2);
#endif
            if (i_to_do)
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
                if ((*i_table_u)->type == 1)
                  mexPrintf("i u[%d]=u[%d]*u[%d]\n", (*i_table_u)->index, (*i_table_u)->op1, (*i_table_u)->op2);
                else if ((*i_table_u)->type == 2)
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
            if ((*table_u)->index>max_u)
              max_u=(*table_u)->index;
            if ((*table_u)->index<min_u)
              min_u=(*table_u)->index;
            if (shifting)
              {
                (*table_u)->index -= y_kmin * u_size;
                if ((*table_u)->index>max_u)
                  max_u=(*table_u)->index;
                if ((*table_u)->index<min_u)
                  min_u=(*table_u)->index;
              }
#ifdef PRINT_OUT
            mexPrintf("push(u[%d])\n", (*table_u)->index);
#endif
            if (i_to_do)
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
SparseMatrix::read_file_table_y(t_table_y **table_y, t_table_y **i_table_y, int *nb_table_y, bool i_to_do, bool shifting, int y_kmin, int y_kmax, int u_size, int y_size)
{
  int i, k;
  SaveCode.read(reinterpret_cast<char *>(nb_table_y), sizeof(*nb_table_y));
#ifdef PRINT_OUT
  mexPrintf("nb_table_y=%d\n", *nb_table_y);
  //mexPrintf("y_size=%d, u_size=%d, y_kmin=%d, y_kmax=%d\n", y_size, u_size, y_kmin, y_kmax);
#endif
  (*table_y) = (t_table_y*)mxMalloc((*nb_table_y) * sizeof(t_table_y));
  for (i = 0;i < *nb_table_y;i++)
    {
      SaveCode.read(reinterpret_cast<char *>(&((*table_y)[i].index)), sizeof((*table_y)[i].index));
      SaveCode.read(reinterpret_cast<char *>(&((*table_y)[i].nb)), sizeof((*table_y)[i].nb));
      if (shifting)
        (*table_y)[i].index -= y_kmin * y_size;
#ifdef PRINT_OUT
      mexPrintf("table_y[i].nb=%d\n", (*table_y)[i].nb);
      mexPrintf("y[%d]=", (*table_y)[i].index);
#endif
      (*table_y)[i].u_index = (int*)mxMalloc((*table_y)[i].nb * sizeof(int));
      (*table_y)[i].y_index = (int*)mxMalloc((*table_y)[i].nb * sizeof(int));
      for (k = 0;k < (*table_y)[i].nb;k++)
        {
          SaveCode.read(reinterpret_cast<char *>(&((*table_y)[i].u_index[k])), sizeof((*table_y)[i].u_index[k]));
          SaveCode.read(reinterpret_cast<char *>(&((*table_y)[i].y_index[k])), sizeof((*table_y)[i].y_index[k]));
          if (shifting)
            {
              (*table_y)[i].u_index[k] -= y_kmin * u_size;
              if (((*table_y)[i].y_index[k] > y_size*y_kmin) && ((*table_y)[i].y_index[k] < y_size*(2*y_kmin + y_kmax + 2)))
                {
                  (*table_y)[i].y_index[k] -= y_kmin * y_size;
                }
            }
#ifdef PRINT_OUT
          if (k < (*table_y)[i].nb - 1)
            mexPrintf("u[%d]*y[%d]+", (*table_y)[i].u_index[k], (*table_y)[i].y_index[k]);
          else
            mexPrintf("u[%d]*y[%d]\n", (*table_y)[i].u_index[k], (*table_y)[i].y_index[k]);
#endif
        }
    }
#ifdef PRINT_OUT
  mexPrintf("*nb_table_y=%d\n", *nb_table_y);
#endif
  if (i_to_do)
    {
      *i_table_y = (t_table_y*)mxMalloc((*nb_table_y) * sizeof(t_table_y));
      for (i = 0;i < *nb_table_y;i++)
        {
          SaveCode.read(reinterpret_cast<char *>(&((*i_table_y)[i].index)), sizeof((*i_table_y)[i].index));
          SaveCode.read(reinterpret_cast<char *>(&((*i_table_y)[i].nb)), sizeof((*i_table_y)[i].nb));
#ifdef PRINT_OUT
          mexPrintf("(*i_table_y)[i].nb=%d\n", (*i_table_y)[i].nb);
          mexPrintf("y_i[%d]=", (*i_table_y)[i].index);
#endif
          (*i_table_y)[i].u_index = (int*)mxMalloc((*i_table_y)[i].nb * sizeof(int));
          (*i_table_y)[i].y_index = (int*)mxMalloc((*i_table_y)[i].nb * sizeof(int));
          for (k = 0;k < (*i_table_y)[i].nb;k++)
            {
              SaveCode.read(reinterpret_cast<char *>(&((*i_table_y)[i].u_index[k])), sizeof((*i_table_y)[i].u_index[k]));
              SaveCode.read(reinterpret_cast<char *>(&((*i_table_y)[i].y_index[k])), sizeof((*i_table_y)[i].y_index[k]));
#ifdef PRINT_OUT
              if (k < (*i_table_y)[i].nb - 1)
                mexPrintf("u[%d]*y[%d]+", (*i_table_y)[i].u_index[k], (*i_table_y)[i].y_index[k]);
              else
                mexPrintf("u[%d]*y[%d]\n", (*i_table_y)[i].u_index[k], (*i_table_y)[i].y_index[k]);
#endif
            }
        }
    }
}





void
SparseMatrix::Read_file(std::string file_name, int periods, int u_size, int y_size, int y_kmin, int y_kmax, int &nb_endo, int &u_count, int &u_count_init, double *u)
{
  int nb_add_u_count = 0;
  filename=file_name;
  //mexPrintf("-------------------------------------------------\n");
#ifdef PRINT_OUT
  mexPrintf("min_u(initial)=%d\n",min_u);
  mexPrintf("%s\n", file_name.c_str());
#endif
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
  SaveCode.read(reinterpret_cast<char *>(&nb_endo), sizeof(nb_endo));
  SaveCode.read(reinterpret_cast<char *>(&u_count), sizeof(u_count));
  SaveCode.read(reinterpret_cast<char *>(&u_count_init), sizeof(u_count_init));
#ifdef PRINT_OUT
  mexPrintf("nb_endo=%d\n", nb_endo);
  mexPrintf("u_count=%d\n", u_count);
  mexPrintf("u_count_init=%d\n", u_count_init);
  //mexPrintf("first table_u\n");
#endif
  read_file_table_u(&first_table_u, &F_first_table_u, &first_i_table_u, &F_first_i_table_u, &nb_first_table_u, true, false, &nb_add_u_count, y_kmin, y_kmax, u_size);
#ifdef PRINT_OUT
  mexPrintf("nb_first_table_u=%d\n", nb_first_table_u);
#endif
//mexErrMsgTxt("Exit from Dynare");
#ifdef PRINT_OUT
  mexPrintf("prologue table_u\n");
#endif
  read_file_table_u(&prologue_table_u, &F_prologue_table_u, NULL, NULL, &nb_prologue_table_u, false, false, &nb_add_u_count, y_kmin, y_kmax, u_size);
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
  read_file_table_u(&middle_table_u, &F_middle_table_u, &middle_i_table_u, &F_middle_i_table_u, &nb_middle_table_u, true,  /*true*/false, &nb_add_u_count, y_kmin, y_kmax, u_size);
#ifdef PRINT_OUT
  mexPrintf("nb_middle_table_u=%d\n",nb_middle_table_u);
  //mexPrintf("last table_u\n");
#endif
  read_file_table_u(&last_table_u, &F_last_table_u, NULL, NULL, &nb_last_table_u, false, false, &nb_add_u_count, y_kmin, y_kmax, u_size);
#ifdef PRINT_OUT
  mexPrintf("->nb_last_table_u=%d\n", nb_last_table_u);
  mexPrintf("i=%d\n", i);
  mexPrintf("going to read prologue_table_y\n");
#endif
  read_file_table_y(&prologue_table_y, NULL, &nb_prologue_table_y, false, false, y_kmin, y_kmax, u_size, y_size);
#ifdef PRINT_OUT
  mexPrintf("nb_prologue_table_y=%d\n", nb_prologue_table_y);
  mexPrintf("going to read first_table_y\n");
#endif
  read_file_table_y(&first_table_y, NULL, &nb_first_table_y, false, false, y_kmin, y_kmax, u_size, y_size);
#ifdef PRINT_OUT
  mexPrintf("nb_first_table_y=%d\n", nb_first_table_y);
  mexPrintf("going to read middle_table_y\n");
#endif
  read_file_table_y(&middle_table_y, &middle_i_table_y, &nb_middle_table_y, true,  /*true*/false, y_kmin, y_kmax, u_size, y_size);
#ifdef PRINT_OUT
  mexPrintf("nb_middle_table_y=%d\n", nb_middle_table_y);
  mexPrintf("going to read last_table_y\n");
#endif
  read_file_table_y(&last_table_y, NULL, &nb_last_table_y, false, false, y_kmin, y_kmax, u_size, y_size);
#ifdef PRINT_OUT
  mexPrintf("nb_last_table_y=%d\n", nb_last_table_y);
  mexPrintf("->nb_last_table_y=%d\n", nb_last_table_y);
  mexPrintf("max_u=%d\n",max_u);
  mexPrintf("min_u=%d\n",min_u);
#endif
  if (nb_last_table_u > 0)
    {
#ifdef PRINT_OUT
      mexPrintf("y_size=%d, periods=%d, y_kmin=%d, y_kmax=%d\n", y_size, periods, y_kmin, y_kmax);
      mexPrintf("u=mxMalloc(%d)\n", max(u_count + 1,max_u+1));
#endif
      u = (double*)mxMalloc(max(u_count + 1,max_u+1) * sizeof(double));
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
  if (u == NULL)
    {
      mexPrintf("memory exhausted\n");
      mexEvalString("st=fclose('all');clear all;");
      mexErrMsgTxt("Exit from Dynare");
    }
  // mexErrMsgTxt("Exit from Dynare");
}


void
SparseMatrix::close_SaveCode()
{
  if (SaveCode.is_open())
    SaveCode.close();
}


/*void
SparseMatrix::initialize(int periods_arg, int nb_endo_arg, int y_kmin_arg, int y_kmax_arg, int y_size_arg, int u_count_arg, int u_count_init_arg, double *u_arg, double *y_arg, double *ya_arg, double slowc_arg, int y_decal_arg, double markowitz_c_arg, double res1_arg, double res2_arg, double max_res_arg)
{
  periods=periods_arg;
  nb_endo=nb_endo_arg;
  y_kmin=y_kmin_arg;
  y_kmax=y_kmax_arg;
  y_size=y_size_arg;
  u_count=u_count_arg;
  u_count_init=u_count_init_arg;
  u=u_arg;
  y=y_arg;
  ya=ya_arg;
  slowc=slowc_arg;
  slowc_save=slowc;
  y_decal=y_decal_arg;
  markowitz_c=markowitz_c_arg;
  res1=res1_arg;
  res2=res2_arg;
  max_res=max_res_arg;
}
*/
