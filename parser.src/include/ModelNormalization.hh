#ifndef MODELNORMALIZATION
#define MODELNORMALIZATION
#include "SymbolTable.hh"
#include "CodeInterpreter.hh"


typedef struct Edge
{
  Edge *next;
  int Vertex_Index;
};

typedef struct Equation_vertex
{
  Edge *First_Edge;
  Edge *Next_Edge;
  int  matched,index;
};

typedef struct Equation_set
{
  Equation_vertex *Number;
  int size;
  int edges;
};

typedef struct simple
{
  int index, block;
  bool available;
};

typedef std::string (*t_getNameByID)(Type type, int id);

class Normalization
{
private:
  typedef struct Variable_vertex
  {
    int  matched;
  };
  typedef struct Variable_set
  {
    Variable_vertex *Number;
    int size;
  };
  typedef struct t_Heap
  {
    int u;        /* vertex */
    int i_parent; /* index in t_Heap of parent vertex in tree of u */
    int v;        /* current matched of u */
  };
public:
  Normalization(const SymbolTable &symbol_table_arg);
  ~Normalization();
  bool Normalize(int n, int prologue, int epilogue, bool* IM, simple* Index_Var_IM, Equation_set* Equation,bool mixing, bool* IM_s);
  void Gr_to_IM_basic(int n0, int prologue, int epilogue, bool* IM, Equation_set *Equation,bool transpose);
  t_getNameByID getnamebyID;
  const SymbolTable &symbol_table;
  void Set_fp_verbose(bool ok);
private:
  void IM_to_Gr(int n0, int prologue, int epilogue, bool* IM, Equation_set *Equation, Variable_set *Variable );
  void Inits(Equation_set *Equation);
  void UpdatePath(Equation_set *Equation, Variable_set *Variable, int i1, int i2);
  void FindAugmentingPaths(Equation_set *Equation, Variable_set *Variable);
  void CheapMatching(Equation_set *Equation, Variable_set *Variable);
  void MaximumMatching(Equation_set *Equation, Variable_set *Variable);
  int MeasureMatching(Equation_set *Equation);
  void OutputMatching(Equation_set* Equation);
  void Gr_to_IM(int n0, int prologue, int epilogue, bool* IM, simple* Index_Var_IM, Equation_set *Equation,bool mixing, bool* IM_s);
  void Free_Equation(int n, Equation_set* Equation);
  void Free_Other(Variable_set* Variable);
  void Free_All(int n, Equation_set* Equation, Variable_set* Variable);
  void ErrorHandling(int n, bool* IM, simple* Index_Equ_IM);
  int eq, eex;
  int IndexUnmatched;
  bool fp_verbose;
  bool* visited;
  t_Heap* Local_Heap;
};
#endif
