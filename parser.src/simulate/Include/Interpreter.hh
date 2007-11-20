#ifndef INTERPRETER_HH_INCLUDED
#define INTERPRETER_HH_INCLUDED

#include <stack>
#include <vector>
#include <string>
#include <cmath>
#include "CodeInterpreter.hh"
#include "SymbolTableTypes.hh"
#include "ExprNode.hh"
//#include "Mem_Mngr.hh"
#include "SparseMatrix.hh"
#include "LinBCG.hh"
#include "mex.h"

//#define DEBUGC

using namespace std;

typedef struct Block_contain_type
{
  int Equation, Variable, Own_Derivative;
};

typedef struct Block_type
{
  int begin, end, size, type;
};

#define pow_ pow
#define get_code_int          *((int*)Code); Code+=sizeof(int);
#define get_code_double       *((double*)Code); Code+=sizeof(double);
#define get_code_pdouble      (double*)Code; Code+=sizeof(double);
#define get_code_bool         *((bool*)(Code++))
#define get_code_char         *((char*)(Code++))
#define get_code_pos          (int)Code
#define get_code_pointer      Code
#define set_code_pointer(pos) Code=pos


class Interpreter : SparseMatrix
{
  protected :
    double pow1(double a, double b);
    void compute_block_time();
    void simulate_a_block(int size,int type, string file_name, string bin_basename, bool Gaussian_Elimination);
    double *T;
    vector<Block_contain_type> Block_Contain;
    vector<Block_type> Block;
    char *Code;
    stack<double> Stack;
    int Block_Count, Per_u_, Per_y_;
    int it_, nb_row_x, nb_row_xd, maxit_, size_of_direction;
    double *g1, *r;
    //double max_res, res1, res2, slowc, markowitz_c;
    double solve_tolf;
    bool GaussSeidel;
    double *x, *params;
    //double *y, *ya, *x, *direction;
    map<pair<pair<int, int> ,int>, int> IM_i;
    string filename;
  public :
    //ReadBinFile read_bin_file;

    Interpreter(double *params_arg, double *y_arg, double *ya_arg, double *x_arg, double *direction_arg, int y_size_arg, int nb_row_x_arg,
                int nb_row_xd_arg, int periods_arg, int y_kmin_arg, int y_kmax_arg, int maxit_arg_, double solve_tolf_arg, int size_o_direction_arg,
                double slowc_arg, int y_decal_arg, double markowitz_c_arg, string &filename_arg);
    void compute_blocks(string file_name, string bin_basename);
};


#endif
