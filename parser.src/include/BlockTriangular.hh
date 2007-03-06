#ifndef BLOCKTRIANGULAR_H
#define BLOCKTRIANGULAR_H
//------------------------------------------------------------------------------
/*! \file
  \version 1.0
  \date 16/07/2006
  \par This file defines the BlockTriangular class.
*/
//------------------------------------------------------------------------------
#include <string>
#include "ExprNode.hh"
#include "SymbolTable.hh"
#include "ModelNormalization.hh"
#include "ModelBlocks.hh"
//------------------------------------------------------------------------------
/*!
  \class  BlockTriangular
  \brief  Creat the incidence matrice and reorder the model's equations.
*/

#include "ExprNode.hh"

typedef struct List_IM
{
  List_IM* pNext;
  int lead_lag;
  bool* IM;
};


typedef struct vari
{
  int Size;
  int* arc;
  int used_arc;
  int available;
};


class BlockTriangular
{
public:
  Normalization normalization;
  BlockTriangular(const  SymbolTable &symbol_table_arg);
  ~BlockTriangular();
  /*! The incidence matrix for each lead and lags */
  Blocks blocks;
  SymbolTable symbol_table;
  List_IM* Build_IM(int lead_lag);
  List_IM* Get_IM(int lead_lag);
  bool* bGet_IM(int lead_lag);
  void fill_IM(int equation, int variable_endo, int lead_lag);
  void unfill_IM(int equation, int variable_endo, int lead_lag);
  void incidence_matrix() const;
  void init_incidence_matrix(int nb_endo);
  void Print_IM(int n);
  void Free_IM(List_IM* First_IM);
  void Print_SIM(bool* IM, int n);
  void Normalize_and_BlockDecompose_Static_Model();
  void Normalize_and_BlockDecompose_Static_0_Model();
  bool Normalize_and_BlockDecompose(bool* IM, Model_Block* ModelBlock, int n, int* prologue, int* epilogue, simple* Index_Var_IM, simple* Index_Equ_IM, bool Do_Normalization, bool mixing, bool* IM_0 );
  void Normalize_and_BlockDecompose_0();
  void Normalize_and_BlockDecompose_Inside_Earth();
  void Prologue_Epilogue(bool* IM, int* prologue, int* epilogue, int n, simple* Index_Var_IM, simple* Index_Equ_IM);
  void Sort_By_Cols(bool* IM, int d, int f);
  void getMax_Lead_Lag(int var, int equ, int *lead, int *lag);
  void getMax_Lead_Lag_B(int size, int* Equation, int *Variable, int *lead, int *lag);
  void swap_IM_c(bool *SIM, int pos1, int pos2, int pos3, simple* Index_Var_IM, simple* Index_Equ_IM, int n);
  void Allocate_Block(int size, int *count_Equ, int *count_Block, int type, Model_Block * ModelBlock, int* Table, int TableSize);
  void Free_Block(Model_Block* ModelBlock);
  void SetVariableTable(int *Table,int Size,int HSize);
  List_IM *First_IM ;
  List_IM *Last_IM ;
  simple *Index_Equ_IM;
  simple *Index_Var_IM;
  int prologue, epilogue;
  int Model_Max_Lead, Model_Max_Lag, periods;
  bool bt_verbose;
  int endo_nbr, TableSize;
  int* Table;
  Model_Block* ModelBlock;
  inline static std::string BlockType0(int type)
  {
    switch (type)
      {
      case 0:
        return ("SIMULTANEOUS TIME SEPARABLE  ");
        break;
      case 1:
        return ("PROLOGUE                     ");
        break;
      case 2:
        return ("EPILOGUE                     ");
        break;
      case 3:
        return ("SIMULTANEOUS TIME UNSEPARABLE");
        break;
      default:
        return ("UNKNOWN                      ");
        break;
      }
  };
  inline static std::string BlockSim(int type)
  {
    switch (type)
      {
      case 0:
        return ("EVALUATE FOREWARD            ");
        break;
      case 1:
        return ("EVALUATE BACKWARD            ");
        break;
      case 2:
        return ("SOLVE FOREWARD SIMPLE        ");
        break;
      case 3:
        return ("SOLVE BACKWARD SIMPLE        ");
        break;
      case 4:
        return ("SOLVE TWO BOUNDARIES SIMPLE  ");
        break;
      case 5:
        return ("SOLVE FOREWARD COMPLETE      ");
        break;
      case 6:
        return ("SOLVE BACKWARD COMPLETE      ");
        break;
      case 7:
        return ("SOLVE TWO BOUNDARIES COMPLETE");
        break;
      default:
        return ("UNKNOWN                      ");
        break;
      }
  };
  inline static std::string BlockSim_d(int type)
  {
    switch (type)
      {
      case 0:
        return ("EVALUATE_FOREWARD            ");
        break;
      case 1:
        return ("EVALUATE_BACKWARD            ");
        break;
      case 2:
        return ("SOLVE_FOREWARD_SIMPLE        ");
        break;
      case 3:
        return ("SOLVE_BACKWARD_SIMPLE        ");
        break;
      case 4:
        return ("SOLVE_TWO_BOUNDARIES_SIMPLE  ");
        break;
      case 5:
        return ("SOLVE_FOREWARD_COMPLETE      ");
        break;
      case 6:
        return ("SOLVE_BACKWARD_COMPLETE      ");
        break;
      case 7:
        return ("SOLVE_TWO_BOUNDARIES_COMPLETE");
        break;
      default:
        return ("UNKNOWN                      ");
        break;
      }
  };

};
//------------------------------------------------------------------------------
#endif
