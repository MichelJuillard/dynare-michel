#include <iostream>
#include <map>
#include <string>
#include <algorithm>
#include <vector>
#include <sstream>
using namespace std;
#include "d_tab.h"
#include "Expression.h"


int main(void)
{
	Expression exp; 
	exp.AddConstant("alpha");
	exp.AddToken(0,eExogenous,EXP);				//0
	exp.AddToken(0,eParameter,0,eExogenousDet,PLUS);			//1
	exp.AddToken(0,eTempResult,UMINUS);	//2
	exp.AddToken(1,eExogenous,1,eTempResult,TIMES);			//3
	exp.AddToken(3,eTempResult,0,eNumericalConstant,TIMES);		//4	
	exp.AddToken(4,eTempResult,0,eTempResult,COMMA);		//5	
	exp.AddToken(5,eTempResult,0,eExogenous,COMMA);		//6	  
	exp.AddToken(6,eTempResult,"function1");		//7	  	
	cout << exp.get() << endl;
	
}
