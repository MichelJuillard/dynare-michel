/*! \file 
 \version 1.0
 \date 04/09/2004
 \par This file implements the OutputFile class methodes.
*/
//------------------------------------------------------------------------------
#include <iostream>
#include <sstream>
using namespace std;
//------------------------------------------------------------------------------
#include "OutputFile.h"
#include "SymbolTable.h"
#include "ModelTree.h"
//------------------------------------------------------------------------------
OutputFile::OutputFile()
{
   	clear_all = true;
}
//------------------------------------------------------------------------------
OutputFile::~OutputFile() 
{
	// Empty
}
//------------------------------------------------------------------------------
void OutputFile::Open(string iFileName)
{
	if (iFileName.size())
	{
		mOutputFile.open(iFileName.c_str(),ios::out|ios::binary);
		if (!mOutputFile.is_open())
		{
			cout << "OutputFile::Open : Error : Can't open file " << iFileName
				<< " for writing\n";
			exit(-1);
		}
	}
	else
	{
			cout << "OutputFile::Open : Error : Missing file name\n";
			exit(-1);
	}
	mOutputFile << "%\n% Status : main Dynare file \n%\n";
	mOutputFile << "% Warning : this file is generated automatically by Dynare\n";
	mOutputFile << "%			from model file (.mod)\n\n";
	if (clear_all)
		mOutputFile << "clear all\n";
	mOutputFile << "tic;\n";
	mOutputFile << "global M_ oo_ exedet_ exdet_ recur_ recurs_ \n";
	mOutputFile << "global options_ endval_\n";
	mOutputFile << "global ys0_ recurs0_ ex0_ ct_\n";
	mOutputFile << "options_ = [];\n";
	iFileName.erase(iFileName.size()-2);
	mOutputFile << "M_.fname = '" << iFileName << "';\n";
	mOutputFile << "%\n% Some global variables initialisation\n%\n";
	mOutputFile << "options_.smpl=0;\noptions_.dynatol = 0.00001;\noptions_.maxit_=10;\noptions_.slowc=1;\noptions_.timing=0;\nct_=0;\noptions_.gstep=1e-2;\n";
	mOutputFile << "options_.debug=0\n";
	mOutputFile << "endval_=0;rplottype_=0;\noptions_.initval_file=0;\n";
	mOutputFile << "oo_.exo_simul = [];\n";
	mOutputFile << "oo_.y_simul = [];\n";
	mOutputFile << "oo_.dr = struct([]);\n";
	mOutputFile << "diary off;\nwarning off;\nwarning on;\nwarning backtrace;\n";
	mOutputFile << "\ndelete " << iFileName << ".log;\n";
	mOutputFile << "logname_ = '" << iFileName << ".log';\n";
	mOutputFile << "diary '" << iFileName << ".log';\n";
	if (ModelTree::offset == 0)
	{	
		mOutputFile << "if exist('" << iFileName << "_static.c', 'file') == 2,\n";
		mOutputFile << "   clear " << iFileName << "_static\n";
		mOutputFile << "   mex -O " << iFileName << "_static.c\n";
		mOutputFile << "end\n";
		mOutputFile << "if exist('" << iFileName << "_dynamic.c', 'file') == 2,\n";
		mOutputFile << "   clear " << iFileName << "_dynamic\n";
		mOutputFile << "   mex -O " << iFileName << "_dynamic.c\n";
		mOutputFile << "end\n";
	}
	else
	{
		mOutputFile << "if exist('" << iFileName << "_static.dll', 'file') == 3,\n";
		mOutputFile << "  clear " << iFileName << "_static\n";
		mOutputFile << "  !del " << iFileName << "_static.dll\n";
		mOutputFile << "end\n";
		mOutputFile << "if exist('" << iFileName << "_dynamic.dll', 'file') == 3,\n";
		mOutputFile << "  clear " << iFileName << "_dynamic\n";
		mOutputFile << "  !del " << iFileName << "_dynamic.dll\n";
		mOutputFile << "end\n";
	}
}
//------------------------------------------------------------------------------
void OutputFile::Save(ostringstream& iOutput)
{
	mOutputFile << SymbolTable::get();
	mOutputFile << ModelTree::get();
	mOutputFile << iOutput.str();
	mOutputFile << "\ndisp(['Total computing time : ' sec2hms(round(toc)) ]);\n";
	mOutputFile.close();
}
//------------------------------------------------------------------------------
#ifdef TEST_OUTPUTFILE
#include "NumericalInitialization.h"
#include "ComputingTasks.h"
#include "Expression.h"
#include "Shocks.h"
#include "SigmaeInitialization.h" 
#include "TmpSymbolTable.h"
int main(void)
{
	OutputFile outputfile;
	SymbolTable st;
	Expression	exp;
	NumericalConstants numconst;
	NumericalInitialization numinit;	
	
	outputfile.Open("Test.m");
	
	SymbolTable::AddSymbolDeclar("a",eExogenous);//0
	SymbolTable::AddSymbolDeclar("b",eParameter);//1
	SymbolTable::AddSymbolDeclar("c",eExogenous);//2
	SymbolTable::AddSymbolDeclar("d",eExogenousDet);//3
	SymbolTable::AddSymbolDeclar("x",eParameter);//3      
	SymbolTable::AddSymbolDeclar("y",eExogenous);//3  
		                                                        
	
	numconst.AddConstant("alpha");
	exp.AddToken(0,eExogenous,EXP);				//0
	exp.AddToken(0,eParameter,0,eExogenousDet,PLUS);			//1
	exp.AddToken(0,eTempResult,UMINUS);	//2
	exp.AddToken(1,eExogenous,1,eTempResult,TIMES);			//3
	exp.AddToken(3,eTempResult,0,eNumericalConstant,TIMES);		//4	
	exp.AddToken(4,eTempResult,0,eTempResult,COMMA);		//5	
	exp.AddToken(5,eTempResult,0,eExogenous,COMMA);		//6	  
	exp.AddToken(6,eTempResult,"function1");		//6	  
	exp.set();	
	//cout << exp.get();

	
	numinit.SetConstant("x","1");
	numinit.InitInitval();
	numinit.SetInit("y",exp.get());
	numinit.EndInitval();
	numinit.InitEndval();
	numinit.EndEndval();
	numinit.InitHistval();		
	numinit.SetHist("y",3, exp.get());			
	//cout << numinit.get();

	
	SigmaeInitialization siginit;
	
	siginit.AddExpression("00");
	siginit.EndOfRow();
	siginit.AddExpression("10");
	siginit.AddExpression("11");
	siginit.EndOfRow();
	siginit.AddExpression("20");
	siginit.AddExpression("21");
	siginit.AddExpression("22");
	siginit.EndOfRow();
	siginit.AddExpression("30");
	siginit.AddExpression("31");
	siginit.AddExpression("32");
	siginit.AddExpression("33");
	siginit.EndOfRow();
	siginit.set();
	
	TmpSymbolTable tmp_symbol_table1, tmp_symbol_table2;
	

	
	ComputingTasks computing_tasks;
	
	computing_tasks.set();	
	computing_tasks.SetSteady();
	computing_tasks.SetCheck();
	computing_tasks.SetSimul();
	
	tmp_symbol_table1.AddTempSymbol("tmp1");
	tmp_symbol_table1.AddTempSymbol("tmp2");
	tmp_symbol_table1.AddTempSymbol("tmp3");
	tmp_symbol_table1.set("var_list_");

	computing_tasks.SetStochSimul(tmp_symbol_table1.get());
	
	computing_tasks.SetOption("DROP", "500");				
	computing_tasks.RunEstimation();					
	computing_tasks.SetEstimationInit();
	computing_tasks.SetEstimation("d", "init_val", "lo_bound", "up_bound", "prior", "p1", "p2", "p3","p4");
	computing_tasks.SetEstimation("a", "c", "init_val", "lo_bound", "up_bound", "prior", "p1", "p2", "p3", "p4");
	computing_tasks.SetCalibInit();
	computing_tasks.SetCalibVariance();
	computing_tasks.SetCalibCovariance();
	computing_tasks.SetCalibAutoCorrelation();
	computing_tasks.SetCalib();
	
	tmp_symbol_table1.AddTempSymbol("tmp11");
	tmp_symbol_table1.AddTempSymbol("tmp22"); 
	tmp_symbol_table1.AddTempSymbol("tmp33");
	tmp_symbol_table1.set("varl_list_");
	computing_tasks.SetOsr(tmp_symbol_table1.get());				
	
	tmp_symbol_table1.AddTempSymbol("tmp11");
	tmp_symbol_table1.AddTempSymbol("tmp22");
	tmp_symbol_table1.AddTempSymbol("tmp33");
	tmp_symbol_table1.set("var_list_");
	tmp_symbol_table2.AddTempSymbol("tmp4");
	tmp_symbol_table2.AddTempSymbol("tmp5");
	tmp_symbol_table2.AddTempSymbol("tmp6");
	tmp_symbol_table2.set("olr_inst_");
	computing_tasks.SetOlr(tmp_symbol_table1.get(), tmp_symbol_table2.get());
	
	computing_tasks.SetOptimWeightsInit();
	computing_tasks.SetOptimWeights1();
	computing_tasks.SetOptimWeights2();
	outputfile.Save();	
}
#endif
//------------------------------------------------------------------------------
