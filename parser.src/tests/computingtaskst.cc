#include <iostream>
#include <list>
#include <sstream>

using namespace std;

#include "ComputingTasks.h"


int main(void)
{
	TmpSymbolTable tst1, tst2;
	SymbolTable st;
	
	ComputingTasks ct;
	
	tst1.AddTempSymbol("a");
	tst1.AddTempSymbol("b");
	tst1.AddTempSymbol("c");
	tst1.AddTempSymbol("d");
	tst1.AddTempSymbol("e");
	
	tst2.AddTempSymbol("aa");
	tst2.AddTempSymbol("bb");
	tst2.AddTempSymbol("cc");
	tst2.AddTempSymbol("dd");
	tst2.AddTempSymbol("ee");

	ct.Set();
	ct.SetSteady();					  
	ct.SetCheck();						  
	ct.SetSimul();						  
	ct.SetStochSimul(tst1);	 
	ct.SetOption("o1", "o2");			  			  
	ct.SetEstimationInit();                                
	//ct.SetEstimation();                                
	//ct.SetEstimation();                                
	//ct.SetEstimation();                                
	ct.SetCalibInit();                                
	ct.SetCalibVariance();                                
	ct.SetCalibCovariance();                                
	ct.SetCalibAutoCorrelation();                                
	ct.SetCalib();                                
	ct.SetOsr(tst1);				
	ct.SetOlr(tst1,tst2);                       
										    
	ct.SetOptimWeightsInit();                                
	ct.SetOptimWeights1();                                
	ct.SetOptimWeights2();                
	ct.RunEstimation();	                
	cout << ct.get();                                
}                                
