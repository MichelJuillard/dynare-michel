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

#include "Mem_Mngr.hh"

Mem_Mngr::Mem_Mngr()
{
  swp_f=false;
  swp_f_b=0;
  //verbose=false;
}
void
Mem_Mngr::Print_heap()
{
  int i;
  mexPrintf("i   :");
  for (i=0;i<CHUNK_SIZE;i++)
    mexPrintf("%3d ",i);
  mexPrintf("\n");
}

void
Mem_Mngr::init_Mem()
{
  Chunk_Stack.clear();
  CHUNK_SIZE=0;
  Nb_CHUNK=0;
  NZE_Mem=NULL;
  NZE_Mem_add=NULL;
  CHUNK_heap_pos=0;
}

void Mem_Mngr::fixe_file_name(string filename_arg)
{
  filename=filename_arg;
}

void
Mem_Mngr::init_CHUNK_BLCK_SIZE(int u_count)
{
  CHUNK_BLCK_SIZE=u_count;
}

NonZeroElem*
Mem_Mngr::mxMalloc_NZE()
{
  int i;
  if (!Chunk_Stack.empty())           /*An unused block of memory available inside the heap*/
    {
      NonZeroElem* p1 = Chunk_Stack.back();
      Chunk_Stack.pop_back();
      return(p1);
    }
  else if (CHUNK_heap_pos<CHUNK_SIZE) /*there is enough allocated memory space available we keep it at the top of the heap*/
    {
      int i=CHUNK_heap_pos++;
      return(NZE_Mem_add[i]);
    }
  else                                /*We have to allocate extra memory space*/
    {
      CHUNK_SIZE+=CHUNK_BLCK_SIZE;
      /*mexPrintf("Allocate %f Ko\n",double(CHUNK_BLCK_SIZE)*double(sizeof(NonZeroElem))/double(1024));
      mexEvalString("drawnow;");*/
      Nb_CHUNK++;
#ifdef MEM_ALLOC_CHK
      mexPrintf("CHUNK_BLCK_SIZE=%d\n",CHUNK_BLCK_SIZE);
#endif
      NZE_Mem=(NonZeroElem*)mxMalloc(CHUNK_BLCK_SIZE*sizeof(NonZeroElem));
      if(!NZE_Mem)
        {
          mexPrintf("Not enough memory available\n");
          mexEvalString("drawnow;");
        }
#ifdef MEM_ALLOC_CHK
      mexPrintf("CHUNK_SIZE=%d\n",CHUNK_SIZE);
#endif
      NZE_Mem_add=(NonZeroElem**)mxRealloc(NZE_Mem_add, CHUNK_SIZE*sizeof(NonZeroElem*));
#ifdef MEM_ALLOC_CHK
      mexPrintf("ok\n");
#endif
      for (i=CHUNK_heap_pos;i<CHUNK_SIZE;i++)
        {
          NZE_Mem_add[i]=(NonZeroElem*)(NZE_Mem+(i-CHUNK_heap_pos));
        }
      i=CHUNK_heap_pos++;
      return(NZE_Mem_add[i]);
    }
}


void
Mem_Mngr::mxFree_NZE(void* pos)
{
  int i, gap;
  /*if(verbose)
    {
      mexPrintf("pos=%x Nb_CHUNK=%d CHUNK_BLCK_SIZE=%d\n",pos,Nb_CHUNK, CHUNK_BLCK_SIZE);
      mexEvalString("drawnow;");
    }
  */
  for (i=0;i<Nb_CHUNK;i++)
    {
      /*if(verbose)
        {
          mexPrintf("i=%d\n",i);
          mexEvalString("drawnow;");
          mexPrintf("NZE_Mem_add[i*CHUNK_BLCK_SIZE]=%d\n",NZE_Mem_add[i*CHUNK_BLCK_SIZE]);
          mexEvalString("drawnow;");
        }*/
      gap=((long int)(pos)-(long int)(NZE_Mem_add[i*CHUNK_BLCK_SIZE]))/sizeof(NonZeroElem);
      if ((gap<CHUNK_BLCK_SIZE) && (gap>=0))
        break;
    }
  /*if(verbose)
    {
      mexPrintf("push_back()\n");
      mexEvalString("drawnow;");
    }*/
  Chunk_Stack.push_back((NonZeroElem*)pos);
  /*if(verbose)
    {
      mexPrintf("End\n");
      mexEvalString("drawnow;");
    }*/
}


void
Mem_Mngr::write_swp_f(int *save_op_all,long int *nop_all)
{
  swp_f=true;
  swp_f_b++;
  mexPrintf("writing the block %d with size=%d\n",swp_f_b,*nop_all);
  if (!SaveCode_swp.is_open())
    {
      mexPrintf("open the swp file for writing\n");
#ifdef PRINT_OUT
      mexPrintf("file opened\n");
#endif
      SaveCode_swp.open((filename + ".swp").c_str(), std::ios::out | std::ios::binary);
      if (!SaveCode_swp.is_open())
        {
          mexPrintf("Error : Can't open file \"%s\" for writing\n", (filename + ".swp").c_str());
          mexEvalString("st=fclose('all');clear all;");
          mexErrMsgTxt("Exit from Dynare");
        }
#ifdef PRINT_OUT
      mexPrintf("done\n");
#endif
    }
  SaveCode_swp.write(reinterpret_cast<char *>(nop_all), sizeof(*nop_all));
  SaveCode_swp.write(reinterpret_cast<char *>(save_op_all), (*nop_all)*sizeof(int));
  (*nop_all)=0;
}

bool
Mem_Mngr::read_swp_f(int **save_op_all,long int *nop_all)
{
  int j;
  swp_f=true;
  if (!SaveCode_swp.is_open())
    {
#ifdef PRINT_OUT
      mexPrintf("file opened\n");
#endif
      mexPrintf("open the file %s\n",(filename + ".swp").c_str());
      SaveCode_swp.open((filename + ".swp").c_str(), std::ios::in | std::ios::binary);
      j=SaveCode_swp.is_open();
      mexPrintf("is_open()=%d\n",j);

      if (!SaveCode_swp.is_open())
        {
          mexPrintf("Error : Can't open file \"%s\" for reading\n", (filename + ".swp").c_str());
          mexEvalString("st=fclose('all');clear all;");
          mexErrMsgTxt("Exit from Dynare");
        }
#ifdef PRINT_OUT
      mexPrintf("done\n");
#endif
      SaveCode_swp.seekg(0);
    }

  j=SaveCode_swp.tellg();
  SaveCode_swp.read(reinterpret_cast<char *>(nop_all), sizeof(*nop_all));
  (*save_op_all)=(int*)mxMalloc((*nop_all)*sizeof(int));
  SaveCode_swp.read(reinterpret_cast<char *>(*save_op_all), (*nop_all)*sizeof(int));
  return(SaveCode_swp.good());
}


void
Mem_Mngr::close_swp_f()
{
  if (SaveCode_swp.is_open())
    {
      SaveCode_swp.close();
      mexPrintf("close the swp file\n");
    }
}

int*
Mem_Mngr::malloc_std(long int nop)
{
  return((int*)malloc(nop*sizeof(int)));
}

int*
Mem_Mngr::realloc_std(int* save_op_o, long int &nopa)
{
  int *save_op=(int*)realloc(save_op_o,nopa*sizeof(int));
  if (!save_op)
    {
      int nopag=int(nopa/3);
      nopa=nopa-nopag;
      while (!save_op && nopag>0)
        {
          nopag=int(nopag*0.66);
          save_op=(int*)realloc(save_op_o,nopa*sizeof(int));
        }
      if (!save_op)
        {
          mexPrintf("Memory exhausted\n");
          mexEvalString("drawnow;");
          mexEvalString("st=fclose('all');clear all;");
          filename+=" stopped";
          mexErrMsgTxt(filename.c_str());
        }
    }
  return(save_op);
}

void
Mem_Mngr::chk_avail_mem(int **save_op_all,long int *nop_all,long int *nopa_all,int add, int t)
  {
    mexPrintf("Error: out of save_op_all[%d] nopa_all=%d t=%d\n",(*nop_all)+add,(*nopa_all),t);
    int tmp_nopa_all=int(1.5*(*nopa_all));
    int *tmp_i;
    if (tmp_nopa_all*sizeof(int)<1024*1024)
      {
        mexPrintf("allocate %d bites save_op_all=%x\n",tmp_nopa_all*sizeof(int),*save_op_all);
        tmp_i=(int*)mxRealloc(*save_op_all,tmp_nopa_all*sizeof(int));
        mexPrintf("tmp_i=");
        mexPrintf("%x\n",tmp_i);
      }
    else
      tmp_i=NULL;
    if (!tmp_i)
      {
        write_swp_f((*save_op_all),nop_all);
      }
    else
      {
        mexPrintf("allocated\n");
        (*save_op_all)=tmp_i;
        (*nopa_all)=tmp_nopa_all;
      }
    mexPrintf("end of chk\n");
  }

void
Mem_Mngr::Free_All()
{
  int i;
  for (i=0;i<Nb_CHUNK;i++)
    {
      mxFree(NZE_Mem_add[i*CHUNK_BLCK_SIZE]);
    }
  mxFree(NZE_Mem_add);
  init_Mem();
}
