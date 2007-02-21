  ////////////////////////////////////////////////////////////////////////
  //                           simulate.c                               //
  //              simulate file designed for Matlab LCC C compiler      //
  //         use LCC_COMPILER option in MODEL command [defailt option]  //
  ////////////////////////////////////////////////////////////////////////


#define PRINT_OUT
//#define FIXE
#define SIZE_OF_INT sizeof(int)

typedef struct t_table_y
{
  int index, nb;
  int *u_index, *y_index;
} t_table_y;

typedef struct t_table_u
{
  struct t_table_u* pNext;
  unsigned char type;
  int index;
  int op1, op2;
} t_table_u;

FILE *SaveCode = NULL, *fopen();


void
read_file_table_u(t_table_u **table_u, t_table_u **F_table_u,
                  t_table_u **i_table_u, t_table_u **F_i_table_u,
                  int *nb_table_u, bool i_to_do, bool shifting,
                  int *nb_add_u_count)
{
  char type;
  int i;
  fread(nb_table_u, SIZE_OF_INT, 1, SaveCode);
#ifdef PRINT_OUT
  mexPrintf("->*nb_table_u=%d\n", *nb_table_u);
#endif
  *table_u = (t_table_u*)mxMalloc(sizeof(t_table_u) - 2 * sizeof(int));
  *F_table_u = *table_u;
  for(i = 0;i < *nb_table_u;i++)
    {
      fread(&type, sizeof(type), 1, SaveCode);
      switch (type)
        {
        case 3:
        case 7:
          (*table_u)->pNext = (t_table_u*)mxMalloc(sizeof(t_table_u) - sizeof(int));
          (*table_u) = (*table_u)->pNext;
          (*table_u)->type = type;
          fread(&(*table_u)->index, SIZE_OF_INT, 1, SaveCode);
          fread(&(*table_u)->op1, SIZE_OF_INT, 1, SaveCode);
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
          break;
        case 1:
        case 2:
        case 6:
          (*table_u)->pNext = (t_table_u*)mxMalloc(sizeof(t_table_u));
          (*table_u) = (*table_u)->pNext;
          (*table_u)->type = type;
          fread(&(*table_u)->index, SIZE_OF_INT, 1, SaveCode);
          fread(&(*table_u)->op1, SIZE_OF_INT, 1, SaveCode);
          fread(&(*table_u)->op2, SIZE_OF_INT, 1, SaveCode);
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
          break;
        case 5:
          (*table_u)->pNext = (t_table_u*)mxMalloc(sizeof(t_table_u) - 2 * sizeof(int));
          (*table_u) = (*table_u)->pNext;
          (*table_u)->type = type;
          fread(&(*table_u)->index, SIZE_OF_INT, 1, SaveCode);
          if(shifting)
            (*table_u)->index -= y_kmin * u_size;
#ifdef PRINT_OUT
          mexPrintf("push(u[%d])\n", (*table_u)->index);
#endif
          break;
        }
    }
  if(i_to_do)
    {
#ifdef PRINT_OUT
      mexPrintf("=>i_table\n");
#endif
      (*i_table_u) = (t_table_u*)mxMalloc(sizeof(t_table_u) - 2 * sizeof(int));
      (*F_i_table_u) = (*i_table_u);
      for(i = 0;i < *nb_table_u;i++)
        {
          fread(&type, sizeof(type), 1, SaveCode);
          switch (type)
            {
            case 3:
            case 7:
              (*i_table_u)->pNext = (t_table_u*)mxMalloc(sizeof(t_table_u) - sizeof(int));
              (*i_table_u) = (*i_table_u)->pNext;
              (*i_table_u)->type = type;
              fread(&(*i_table_u)->index, SIZE_OF_INT, 1, SaveCode);
              fread(&(*i_table_u)->op1, SIZE_OF_INT, 1, SaveCode);
#ifdef FIXE
              (*i_table_u)->index = u_size;
              (*i_table_u)->op1 = u_size;
#endif
#ifdef PRINT_OUT
              if((*i_table_u)->type == 3)
                mexPrintf("u[%d]=1/(1-u[%d])\n", (*i_table_u)->index, (*i_table_u)->op1);
              else
                mexPrintf("u[%d]*=u[%d]\n", (*i_table_u)->index, (*i_table_u)->op1);
#endif
              break;
            case 1:
            case 2:
            case 6:
              (*i_table_u)->pNext = (t_table_u*)mxMalloc(sizeof(t_table_u));
              (*i_table_u) = (*i_table_u)->pNext;
              (*i_table_u)->type = type;
              fread(&(*i_table_u)->index, SIZE_OF_INT, 1, SaveCode);
              fread(&(*i_table_u)->op1, SIZE_OF_INT, 1, SaveCode);
              fread(&(*i_table_u)->op2, SIZE_OF_INT, 1, SaveCode);
#ifdef FIXE
              (*i_table_u)->index = u_size;
              (*i_table_u)->op1 = u_size;
              (*i_table_u)->op2 = u_size;
#endif
#ifdef PRINT_OUT
              if((*i_table_u)->type == 1)
                mexPrintf("u[%d]=u[%d]*u[%d]\n", (*i_table_u)->index, (*i_table_u)->op1, (*i_table_u)->op2);
              else if((*i_table_u)->type == 2)
                mexPrintf("u[%d]+=u[%d]*u[%d]\n", (*i_table_u)->index, (*i_table_u)->op1, (*i_table_u)->op2);
              else
                mexPrintf("u[%d]=1/(1-u[%d]*u[%d])\n", (*i_table_u)->index, (*i_table_u)->op1, (*i_table_u)->op2);
#endif
              break;
            case 5:
              (*i_table_u)->pNext = (t_table_u*)mxMalloc(sizeof(t_table_u) - 2 * sizeof(int));
              (*i_table_u) = (*i_table_u)->pNext;
              (*i_table_u)->type = type;
              fread(&(*i_table_u)->index, SIZE_OF_INT, 1, SaveCode);
#ifdef FIXE
              (*i_table_u)->index = u_size;
#endif
#ifdef PRINT_OUT
              mexPrintf("push(u[%d])\n", (*i_table_u)->index);
#endif
              break;
            }
        }
    }
}

void
read_file_table_y(t_table_y **table_y, t_table_y **i_table_y, int *nb_table_y, bool i_to_do, bool shifting)
{
  int i, k;
  fread(nb_table_y, SIZE_OF_INT, 1, SaveCode);

  mexPrintf("nb_table_y=%d\n", *nb_table_y);
#ifdef PRINT_OUT
  mexPrintf("nb_table_y=%d\n", *nb_table_y);
  mexPrintf("y_size=%d, u_size=%d, y_kmin=%d, y_kmax=%d\n", y_size, u_size, y_kmin, y_kmax);
#endif
  (*table_y) = (t_table_y*)mxMalloc((*nb_table_y) * sizeof(t_table_y));
  for(i = 0;i < *nb_table_y;i++)
    {
      fread(&((*table_y)[i].index), SIZE_OF_INT, 1, SaveCode);
      fread(&((*table_y)[i].nb), SIZE_OF_INT, 1, SaveCode);
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
          fread(&((*table_y)[i].u_index[k]), SIZE_OF_INT, 1, SaveCode);
          fread(&((*table_y)[i].y_index[k]), SIZE_OF_INT, 1, SaveCode);
          /*if(shifting)
            {
              (*table_y)[i].u_index[k] -= y_kmin * u_size;
              if(((*table_y)[i].y_index[k] > y_size*y_kmin) && ((*table_y)[i].y_index[k] < y_size*(2*y_kmin + y_kmax + 2)))
                {
                  (*table_y)[i].y_index[k] -= y_kmin * y_size;
                }
            }*/
#ifdef PRINT_OUT
          //mexPrintf("sizeof((*i_table_y)[i].y_index[k])=%d\n",sizeof((*i_table_y)[i].y_index[k]));
          if(k < (*table_y)[i].nb - 1)
            mexPrintf("u[%d]*y[%d]+", (*table_y)[i].u_index[k], (*table_y)[i].y_index[k]);
          else
            mexPrintf("u[%d]*y[%d]\n", (*table_y)[i].u_index[k], (*table_y)[i].y_index[k]);
#endif
        }
    }
#ifdef PRINT_OUT
  mexPrintf("*nb_table_y=%d\n", *nb_table_y);
  mexPrintf("i_to_do=%d\n", i_to_do);
#endif
  if(i_to_do)
    {
      *i_table_y = (t_table_y*)mxMalloc((*nb_table_y) * sizeof(t_table_y));
      for(i = 0;i < *nb_table_y;i++)
        {
          fread(&((*i_table_y)[i].index), SIZE_OF_INT, 1, SaveCode);
          fread(&((*i_table_y)[i].nb), SIZE_OF_INT, 1, SaveCode);
#ifdef PRINT_OUT
          mexPrintf("(*i_table_y)[i].nb=%d\n", (*i_table_y)[i].nb);
          mexPrintf("y_i[%d]=", (*i_table_y)[i].index);
#endif
          (*i_table_y)[i].u_index = (int*)mxMalloc((*i_table_y)[i].nb * sizeof(int));
          (*i_table_y)[i].y_index = (int*)mxMalloc((*i_table_y)[i].nb * sizeof(int));
          for(k = 0;k < (*i_table_y)[i].nb;k++)
            {
              fread(&((*i_table_y)[i].u_index[k]), SIZE_OF_INT, 1, SaveCode);
              fread(&((*i_table_y)[i].y_index[k]), SIZE_OF_INT, 1, SaveCode);
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


int i, j, k, nb_endo, u_count, u_count_init;
int nb_prologue_table_u, nb_first_table_u, nb_middle_table_u, nb_last_table_u;
int nb_prologue_table_y, nb_first_table_y, nb_middle_table_y, nb_last_table_y;
int first_count_loop, middle_count_loop;
char type;
t_table_u *prologue_table_u, *first_table_u, *first_i_table_u, *middle_table_u, *middle_i_table_u, *last_table_u;
t_table_y *prologue_table_y, *first_table_y, *middle_table_y, *middle_i_table_y, *last_table_y;
t_table_u *F_prologue_table_u, *F_first_table_u, *F_first_i_table_u, *F_middle_table_u, *F_middle_i_table_u, *F_last_table_u;


void
Read_file(char* a_file_name, int periods, int u_size1, int y_size, int y_kmin, int y_kmax)
{
  int nb_add_u_count = 0;
  char *file_name=(char*)malloc(strlen(a_file_name)+5);
  file_name=strcat(a_file_name,".bin");
  u_size = u_size1;
#ifdef PRINT_OUT
  mexPrintf("file_name=%s\n", file_name);
#endif
  if(!SaveCode)
    {
#ifdef PRINT_OUT
      mexPrintf("file opened\n");
#endif
      if(!(SaveCode=fopen(file_name,"r")))
        {
          mexPrintf("Error : Can't open file \"%s\" for reading\n", file_name);
          mexErrMsgTxt("Exit from Dynare");
        }
#ifdef PRINT_OUT
      mexPrintf("done\n");
#endif
    }
  fread(&nb_endo, SIZE_OF_INT, 1, SaveCode);
  fread(&u_count, SIZE_OF_INT, 1, SaveCode);
  fread(&u_count_init, SIZE_OF_INT, 1, SaveCode);

#ifdef PRINT_OUT
  mexPrintf("nb_endo=%d, u_count=%d, u_count_init=%d\n",nb_endo,u_count, u_count_init);
  mexPrintf("prologue table_u\n");
#endif
  read_file_table_u(&prologue_table_u, &F_prologue_table_u, NULL, NULL, &nb_prologue_table_u, false, false, &nb_add_u_count);
#ifdef PRINT_OUT
  mexPrintf("nb_prologue_table_u=%d\n", nb_prologue_table_u);
  mexPrintf("first table_u\n");
#endif
  read_file_table_u(&first_table_u, &F_first_table_u, &first_i_table_u, &F_first_i_table_u, &nb_first_table_u, true, false, &nb_add_u_count);
#ifdef PRINT_OUT
  mexPrintf("nb_first_table_u=%d\n", nb_first_table_u);
  mexPrintf("nb_endo=%d\n", nb_endo);
  mexPrintf("nb_first_table_u=%d\n", nb_first_table_u);
  mexPrintf("u_count=%d\n", u_count);
#endif
  fread(&middle_count_loop, SIZE_OF_INT, 1, SaveCode);
#ifdef PRINT_OUT
  mexPrintf("middle table_u\n");
#endif
  read_file_table_u(&middle_table_u, &F_middle_table_u, &middle_i_table_u, &F_middle_i_table_u, &nb_middle_table_u, true,  /*true*/false, &nb_add_u_count);
#ifdef PRINT_OUT
  mexPrintf("last table_u\n");
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
}


typedef struct t_Stack
{
  double uu;
  struct t_Stack *pPrev;
} t_Stack;

t_Stack *Stack=NULL, *tmp_Stack;

void push(double uu)
{
  if(!Stack)
    {
      Stack=(t_Stack*)malloc(sizeof(t_Stack));
      Stack->pPrev=NULL;
    }
  else
    {
      tmp_Stack=Stack;
      Stack=(t_Stack*)malloc(sizeof(t_Stack));
      Stack->pPrev=tmp_Stack;
    }
   Stack->uu=uu;
}

double pop()
{
  double uu;
  uu=Stack->uu;
  Stack=Stack->pPrev;
  return(uu);
}

void
simulate(int blck, int y_size, int it_, int y_kmin, int y_kmax)
{
  int i, j, k, l, m, m1, nop;
  int period = it_ * y_size, s_middle_count_loop = 0 ;
  pctimer_t t1 = pctimer(), t2;
  double uu;
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
      mexPrintf("//                          Simulate                                  \\\n");
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
              push(u[first_table_u->index + j*first_i_table_u->index]);
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
          first_table_u = first_table_u->pNext;
          first_i_table_u = first_i_table_u->pNext;
          nop++;
        }
    }
#ifdef PRINT_OUT
  mexPrintf("prologue\n");
#endif
  prologue_table_u = F_prologue_table_u->pNext;
  for(i = 0;i < nb_prologue_table_u;i++)
    {
      switch (prologue_table_u->type)
        {
        case 1:
          u[prologue_table_u->index ] = u[prologue_table_u->op1 ] * u[prologue_table_u->op2 ];
#ifdef PRINT_OUT
          mexPrintf("u[%d]=u[%d]*u[%d]=%f\n", prologue_table_u->index , prologue_table_u->op1 , prologue_table_u->op2 , u[prologue_table_u->index ]);
#endif
          break;
        case 2:
          u[prologue_table_u->index ] += u[prologue_table_u->op1 ] * u[prologue_table_u->op2 ];
#ifdef PRINT_OUT
          mexPrintf("u[%d]+=u[%d]*u[%d]=%f\n" , prologue_table_u->index , prologue_table_u->op1 , prologue_table_u->op2 , u[prologue_table_u->index ]);
#endif
          break;
        case 3:
          u[prologue_table_u->index ] = 1 / ( -u[prologue_table_u->op1 ]);
#ifdef PRINT_OUT
          mexPrintf("u[%d]=1/(-u[%d])=%f\n", prologue_table_u->index, prologue_table_u->op1, u[prologue_table_u->index]);
#endif
          break;
        case 5:
          push(u[prologue_table_u->index]);
#ifdef PRINT_OUT
          mexPrintf("push(u[%d])\n", prologue_table_u->index );
#endif
          break;
        case 6:
          u[prologue_table_u->index ] = 1 / (1 - u[prologue_table_u->op1] * u[prologue_table_u->op2]);
#ifdef PRINT_OUT
          mexPrintf("u[%d]=1/(1-u[%d]*u[%d])=%f\n", prologue_table_u->index, prologue_table_u->op1, prologue_table_u->op2, u[prologue_table_u->index]);
#endif
          break;
        case 7:
          u[prologue_table_u->index] *= u[prologue_table_u->op1];
#ifdef PRINT_OUT
          mexPrintf("u[%d]*=u[%d]=%f\n", prologue_table_u->index, prologue_table_u->op1, u[prologue_table_u->index]);
#endif
          break;
        }
      prologue_table_u = prologue_table_u->pNext;
      nop++;
    }
#ifdef PRINT_OUT
  mexPrintf("middle_u (s_middle_count_loop=%d\n", s_middle_count_loop);
#endif
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
              push(u[middle_table_u->index + j*middle_i_table_u->index]);
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
          push(u[last_table_u->index]);
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
      last_table_u = last_table_u->pNext;
      nop++;
    }
#ifdef PRINT_OUT
  mexPrintf("nb_last_table_y=%d\n",nb_last_table_y);
#endif
  for(i = nb_last_table_y - 1;i >= 0;i--)
    {
      k = last_table_y[i].index;
      y[period + k] = 0;
#ifdef PRINT_OUT
      mexPrintf("it_=%d\n", it_);
      mexPrintf("->y[it_*y_size+%d]=y[%d]=", k, it_*y_size + k);
#endif
      for(j = last_table_y[i].nb - 1;j >= 0;j--)
        {
          uu = pop();
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
            y[period + k] += uu * y[period + m];
          else
            y[period + k] += uu;
        }
#ifdef PRINT_OUT
      mexPrintf("=%f\n", y[period + k]);
#endif
      nop++;
    }
#ifdef PRINT_OUT
  mexPrintf("s_middle_count_loop=%d\n",s_middle_count_loop);
  mexPrintf("nb_middle_table_y=%d\n",nb_middle_table_y);
#endif
  for(j = s_middle_count_loop - y_kmin - 1;j >= 0;j--)
    {
      for(i = nb_middle_table_y - 1;i >= 0;i--)
        {
          k = middle_table_y[i].index + j * middle_i_table_y[i].index;
          y[k] = 0;
#ifdef PRINT_OUT
          mexPrintf("y[%d]=", k );
#endif
          for(l = middle_table_y[i].nb - 1;l >= 0;l--)
            {
              uu = pop();
              m = middle_table_y[i].y_index[l] + j * middle_i_table_y[i].y_index[l];
#ifdef PRINT_OUT
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
#endif
              if(m >= 0)
                y[k] += uu * y[m];
              else
                y[k] += uu;
            }
#ifdef PRINT_OUT
          mexPrintf("=%f\n", y[k]);
#endif
          nop++;
        }
    }
  for(i = nb_prologue_table_y - 1;i >= 0;i--)
    {
      k = prologue_table_y[i].index;
      y[k] = 0;
#ifdef PRINT_OUT
      mexPrintf("y[%d]=", k);
#endif
      for(j = prologue_table_y[i].nb - 1;j >= 0;j--)
        {
          uu = pop();
#ifdef PRINT_OUT
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
#endif
          if(prologue_table_y[i].y_index[j] >= 0)
            y[k] += uu * y[prologue_table_y[i].y_index[j]];
          else
            y[k] += uu;
        }
#ifdef PRINT_OUT
      mexPrintf("=%f\n", y[k]);
#endif
      nop++;
    }
  for(i = nb_first_table_y - 1;i >= 0;i--)
    {
      k = first_table_y[i].index;
      y[k] = 0;
#ifdef PRINT_OUT
      mexPrintf("y[%d]=", k);
#endif
      for(j = first_table_y[i].nb - 1;j >= 0;j--)
        {
          uu = pop();
#ifdef PRINT_OUT
          if(j > 0)
            mexPrintf("u[%d](%f)*y[%d](%f)+", first_table_y[i].u_index[j], uu, first_table_y[i].y_index[j], y[first_table_y[i].y_index[j]]);
          else
            mexPrintf("u[%d](%f)*y[%d](%f)", first_table_y[i].u_index[j], uu, first_table_y[i].y_index[j], y[first_table_y[i].y_index[j]]);
#endif
          if(m >= 0)
            y[k] += uu * y[first_table_y[i].y_index[j]];
          else
            y[k] += uu;
        }
#ifdef PRINT_OUT
      mexPrintf("=%f\n", y[k]);
#endif
      nop++;
    }
  t2 = pctimer();
  if(nb_first_table_u > 0)
    mexPrintf("(**%f milliseconds u_count : %d  nop : %d **)\n", 1000*(t2 - t1), u_count, nop);
}

