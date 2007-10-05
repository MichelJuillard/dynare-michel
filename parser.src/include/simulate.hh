#ifndef SIMULATE_HH_INCLUDED
#define SIMULATE_HH_INCLUDED

#include <math>
#include <stack>
#include <set>
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>
//#include "pctimer_h.hh"
#include <time.h>
#include <string>
#include <map>
#include <algorithm>
#include "CodeInterpreter.hh"
#include "SymbolTableTypes.hh"
#include "mex.h"
#include "ExprNode.hh"
#define pow_ pow
//#define pow pow1

// typedef struct IM_compact
// {
//   int size, u_init, u_finish, nb_endo;
//   int *u, *Var, *Equ, *Var_Index, *Equ_Index, *Var_dyn_Index;
// };
typedef struct Variable_l
{
  int* Index;
};
typedef struct tBlock
{
    int Size, Sized, Type, Max_Lead, Max_Lag, Simulation_Type, /*icc1_size,*/ Nb_Lead_Lag_Endo;
    int *Variable, *dVariable, *Equation/*, *icc1, *ics*/;
    int *variable_dyn_index, *variable_dyn_leadlag;
    IM_compact *IM_lead_lag;
};

typedef struct tModel_Block
{
    int Size;
    tBlock * List;
};
#define INDIRECT_SIMULATE
#define MARKOVITZ
#define PRINT_OUT_p
//#define RECORD_ALL
//#define DEBUGC
//#define PRINT_OUT
//#define PRINT_OUT_y1
//#define PRINT_u
//#define PRINT_OUT_b
//#define PRINT_OUT_y
//#define DEBUG
//#define EXTENDED
//#define FLOAT
//#define WRITE_u
//#define MEMORY_LEAKS
//#define N_MX_ALLOC
//#define MEM_ALLOC_CHK
#define NEW_ALLOC
//#define PROFILER
#ifdef EXTENDED
typedef long double longd;
#else
typedef double longd;
#endif


using namespace std;

const int IFLD  =0;
const int IFDIV =1;
const int IFLESS=2;
const int IFSUB =3;
const int IFLDZ =4;
const int IFMUL =5;
const int IFSTP =6;
const int IFADD =7;
const longd eps=1e-7;

typedef struct NonZeroElem
  {
    int u_index;
    int r_index, c_index, lag_index;
    NonZeroElem *NZE_R_N, *NZE_C_N;
  };

  typedef struct t_save_op_s
  {
    short int lag, operat;
    int first, second;
  };




std::map<std::pair<std::pair<int, int> ,int>, int> IM_i;
std::multimap<std::pair<int,int>, int> var_in_equ_and_lag_i, equ_in_var_and_lag_i;
int  *index_vara, *pivota=NULL, *save_op_all=NULL, *b, *g_save_op=NULL;
int *pivot, *pivotk;
longd *pivotv, *pivotva=NULL;
bool *line_done;
bool symbolic=true, alt_symbolic=false;
#ifdef RECORD_ALL
bool record_all=false;
#endif
long int nop_all, nopa_all, nop1, nop2;
const longd very_big=1e24;
int Per_y_, Per_u_, it_, nb_row_x, nb_row_xd, u_size, y_size, x_size, y_kmin, y_kmax, y_decal;
int periods, maxit_, max_u=0, min_u=0x7FFFFFFF, g_nop_all=0;
longd *params, markowitz_c, markowitz_c_s;
longd  *u, *y, *x, *r, *g1, *g2, *ya, *direction;
longd slowc, slowc_save, solve_tolf, max_res, res1, res2, res1a=9.0e60;
bool cvg, print_err, swp_f;
int swp_f_b=0;
//pctimer_t t0, t1;
clock_t t0, t1;
int u_count_alloc, u_count_alloc_save, size_of_direction;
int i, j, k, nb_endo, u_count, u_count_init, iter;
const int alt_symbolic_count_max=1;
int alt_symbolic_count=0;
longd err;
std::string filename;
NonZeroElem **FNZE_R, **FNZE_C;
int *NbNZRow, *NbNZCol;

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

class SparseMatrix
  {
  public:
    SparseMatrix();
    void Init(int periods, int y_kmin, int y_kmax, int Size, std::map<std::pair<std::pair<int, int> ,int>, int> IM);
    void ShortInit(int periods, int y_kmin, int y_kmax, int Size, std::map<std::pair<std::pair<int, int> ,int>, int> IM);
    void End(int Size);
    void Insert(int r, int c, int u_index, int lag_index);
    void Delete(int r,int c, int Size, int *b);
    int At_Row(int r, NonZeroElem **first);
    int At_Pos(int r, int c, NonZeroElem **first);
    int At_Col(int c, NonZeroElem **first);
    int At_Col(int c, int lag, NonZeroElem **first);
    int NRow(int r);
    int NCol(int c);
    int Union_Row(int row1, int row2);
    void Print(int Size,int *b);
    int Get_u();
    void Delete_u(int pos);
    void Clear_u();
    void Print_u();
    void mxFree_NZE(void* pos);
    NonZeroElem* mxMalloc_NZE();
    int complete(int beg_t, int Size, int periods, int *b);
    longd bksub( int tbreak, int last_period, int Size, longd slowc_l
#ifdef PROFILER
    , long int *nmul
#endif
    );
    void Print_heap();
    void init_Mem();
    void run_triangular(int nop_all,int *op_all);
    void run_it(int nop_all,int *op_all);
    void run_u_period1(int periods);


    vector<int> u_liste;
    //set<NonZeroElem*> mem_NZE;
    //vector<void*> used_in_chunk;
    vector<NonZeroElem*> Chunk_Stack;
    //map<vector,pair<void*,size>> chunk_vector_list;
    int CHUNK_SIZE, CHUNK_BLCK_SIZE, Nb_CHUNK;
    int CHUNK_heap_pos/*, CHUNK_heap_max_size*/;
    //int *NZE_pos, *NZE_rpos;
    NonZeroElem** NZE_Mem_add;
    //bool *NZE_Mem_available;
    NonZeroElem* NZE_Mem;
  };

#define get_code_int          *((int*)Code); Code+=sizeof(int);
#define get_code_double       *((double*)Code); Code+=sizeof(double);
#define get_code_pdouble      (double*)Code; Code+=sizeof(double);
#define get_code_bool         *((bool*)(Code++))
#define get_code_char         *((char*)(Code++))
#define get_code_pos          (int)Code
#define get_code_pointer      Code
#define set_code_pointer(pos) Code=pos
typedef struct Block_type
{
  int begin, end, size, type;
};
typedef struct Block_contain_type
{
  int Equation, Variable, Own_Derivative;
};

class Interpreter
{
  protected :
    void compute_block_time();
    void simulate_a_block(int size,int type, string file_name, string bin_basename);
    longd *T;
    vector<Block_contain_type> Block_Contain;
    vector<Block_type> Block;
    char *Code;
    stack<longd> Stack;
    int Block_Count;
    int _it;
    longd *g1, *r;
    bool GaussSeidel;
  public :
    Interpreter();
    void compute_blocks(string file_name, string bin_basename);
};


std::fstream SaveCode, SaveCode_swp;

#endif // SIMULATE_HH_INCLUDED
