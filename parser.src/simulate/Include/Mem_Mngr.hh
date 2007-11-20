#ifndef MEM_MNGR_HH_INCLUDED
#define MEM_MNGR_HH_INCLUDED

#include <vector>
#include <fstream>
#include "mex.h"
using namespace std;

typedef struct NonZeroElem
  {
    int u_index;
    int r_index, c_index, lag_index;
    NonZeroElem *NZE_R_N, *NZE_C_N;
  };

typedef vector<NonZeroElem*> v_NonZeroElem;

class Mem_Mngr
{
public:
    void write_swp_f(int *save_op_all,long int *nop_all);
    bool read_swp_f(int **save_op_all,long int *nop_all);
    void close_swp_f();
    void Print_heap();
    void init_Mem();
    void mxFree_NZE(void* pos);
    NonZeroElem* mxMalloc_NZE();
    void init_CHUNK_BLCK_SIZE(int u_count);
    void Free_All();
    Mem_Mngr();
    void fixe_file_name(string filename_arg);
    int* malloc_std(long int nop);
    int* realloc_std(int* save_op_o, long int &nopa);
    void chk_avail_mem(int **save_op_all,long int *nop_all,long int *nopa_all,int add, int t);
    bool swp_f;
    //bool verbose;
private:
    v_NonZeroElem Chunk_Stack;
    int CHUNK_SIZE, CHUNK_BLCK_SIZE, Nb_CHUNK;
    int CHUNK_heap_pos/*, CHUNK_heap_max_size*/;
    NonZeroElem** NZE_Mem_add;
    NonZeroElem* NZE_Mem;
    int swp_f_b;
    fstream  SaveCode_swp;
    string filename;
};

#endif
