#include <iostream>
#include <map>
#include <string>
#include <algorithm>
#include <vector>
#include <sstream>
using namespace std;
#include "d_tab.h"
#include "Expression.h"
#include "NumericalInitialization.h"



int main(void)
{


	SymbolTable st;
	Expression	exp;
	
	
	SymbolTable::AddSymbolDeclar("a",eExogenous);//0
	SymbolTable::AddSymbolDeclar("b",eParameter);//1
	SymbolTable::AddSymbolDeclar("c",eExogenous);//2
	SymbolTable::AddSymbolDeclar("d",eExogenousDet);//3
	SymbolTable::AddSymbolDeclar("x",eParameter);//3      
	SymbolTable::AddSymbolDeclar("y",eExogenous);//3  
		                                                        
	
	exp.AddConstant("alpha");
	exp.AddToken(0,eExogenous,EXP);				//0
	exp.AddToken(0,eParameter,0,eExogenousDet,PLUS);			//1
	exp.AddToken(0,eTempResult,UMINUS);	//2
	exp.AddToken(1,eExogenous,1,eTempResult,TIMES);			//3
	exp.AddToken(3,eTempResult,0,eNumericalConstant,TIMES);		//4	
	exp.AddToken(4,eTempResult,0,eTempResult,COMMA);		//5	
	exp.AddToken(5,eTempResult,0,eExogenous,COMMA);		//6	  
	exp.AddToken(6,eTempResult,"function1");		//6	  	
	cout << exp.get();
	//Testing if symbol exists 
	if (!SymbolTable::Exist("x"))
    {
      		cout << "Error : Unknown parameter: " << "x" << endl;
      		exit(-1);
    }
	
	NumericalInitialization numinit;	
	
	numinit.SetConstant("x","ok");
	numinit.InitInitval();
	numinit.SetInit("y",exp.get());
	numinit.EndInitval();
	numinit.InitEndval();
	numinit.EndEndval();
	numinit.InitHistval();		
	numinit.SetHist("y",3, exp.get());			
	cout << numinit.get();
	
	
}
