/*! \file 
 \version 1.0
 \date 04/09/2004
 \par This file implements the OperatorTable class methodes.
*/
//------------------------------------------------------------------------------
#include <map>
#include <string>
using namespace std;
//------------------------------------------------------------------------------
#include "OperatorTable.h"
//------------------------------------------------------------------------------
OperatorTable::OperatorTable()
{			
	operator_table[COMMA].str	= ",";
	operator_table[EQUAL].str 	= "="; 
	operator_table[PLUS].str 	= "+";        
	operator_table[MINUS].str 	= "-";        
	operator_table[TIMES].str 	= "*";        
	operator_table[DIVIDE].str 	= "/";        
	operator_table[UMINUS].str 	= "-";        
	operator_table[POWER].str 	= "^";	
	operator_table[EXP].str 	= "exp";      
	operator_table[LOG].str 	= "log";      
	operator_table[LOG10].str 	= "log10";    
	operator_table[COS].str 	= "cos";      
	operator_table[SIN].str 	= "sin";      
	operator_table[TAN].str 	= "tan";      
	operator_table[ACOS].str 	= "acos";     
	operator_table[ASIN].str 	= "asin";     
	operator_table[ATAN].str 	= "atan";     
	operator_table[COSH].str 	= "cosh";     
	operator_table[SINH].str 	= "sinh";     
	operator_table[TANH].str 	= "tanh";     
	operator_table[ACOSH].str 	= "acosh";    
	operator_table[ASINH].str 	= "asinh";    
	operator_table[ATANH].str 	= "atanh";    
	operator_table[SQRT].str 	= "sqrt";
	operator_table[NAME].str 	= ""; 

 
	operator_table[COMMA].precedence 	= -1;
	operator_table[EQUAL].precedence 	= 0;
	operator_table[PLUS].precedence 	= 0;	
    operator_table[MINUS].precedence 	= 1;
    operator_table[TIMES].precedence 	= 2;
    operator_table[DIVIDE].precedence 	= 3;  
    operator_table[UMINUS].precedence 	= 4;  
    operator_table[POWER].precedence 	= 5;    
    operator_table[EXP].precedence 		=     
    operator_table[LOG].precedence 		=     
    operator_table[LOG10].precedence 	= 
    operator_table[COS].precedence 		=     
    operator_table[SIN].precedence 		=     
    operator_table[TAN].precedence 		=     
    operator_table[ACOS].precedence 	= 
    operator_table[ASIN].precedence 	= 
    operator_table[ATAN].precedence 	= 
    operator_table[COSH].precedence 	= 
    operator_table[SINH].precedence 	= 
    operator_table[TANH].precedence 	= 
    operator_table[ACOSH].precedence 	= 
    operator_table[ASINH].precedence 	= 
    operator_table[ATANH].precedence 	= 
    operator_table[SQRT].precedence 	= 
    operator_table[NAME].precedence 	= 6;


	// Operator costs for M files
	operator_table[COMMA].cost[1] 	= 0;
	operator_table[EQUAL].cost[1] 	= 0;
	operator_table[PLUS].cost[1] 	= 90;	
    operator_table[MINUS].cost[1] 	= 90;
    operator_table[TIMES].cost[1] 	= 90;
    operator_table[DIVIDE].cost[1] 	= 990;  
    operator_table[UMINUS].cost[1] 	= 70;  
    operator_table[POWER].cost[1] 	= 1160;    
    operator_table[EXP].cost[1] 	= 160;    
    operator_table[LOG].cost[1] 	= 300;    
    operator_table[LOG10].cost[1] 	= 16000;
    operator_table[COS].cost[1] 	= 210;    
    operator_table[SIN].cost[1] 	= 210;    
    operator_table[TAN].cost[1] 	= 230;    
    operator_table[ACOS].cost[1] 	= 300;
    operator_table[ASIN].cost[1] 	= 310;
    operator_table[ATAN].cost[1] 	= 140;
    operator_table[COSH].cost[1] 	= 210;
    operator_table[SINH].cost[1] 	= 240;
    operator_table[TANH].cost[1] 	= 190;
    operator_table[ACOSH].cost[1] 	= 770;
    operator_table[ASINH].cost[1] 	= 460;
    operator_table[ATANH].cost[1] 	= 350;
    operator_table[SQRT].cost[1] 	= 570;
    operator_table[NAME].cost[1] 	= 0;

	// Operator costs for C files
	operator_table[COMMA].cost[0] 	= 0;
	operator_table[EQUAL].cost[0] 	= 0;
	operator_table[PLUS].cost[0] 	= 4;	
    operator_table[MINUS].cost[0] 	= 4;
    operator_table[TIMES].cost[0] 	= 4;
    operator_table[DIVIDE].cost[0] 	= 15;  
    operator_table[UMINUS].cost[0] 	= 3;  
    operator_table[POWER].cost[0] 	= 520;    
    operator_table[EXP].cost[0] 	= 210;    
    operator_table[LOG].cost[0] 	= 137;    
    operator_table[LOG10].cost[0] 	= 139;
    operator_table[COS].cost[0] 	= 160;    
    operator_table[SIN].cost[0] 	= 160;    
    operator_table[TAN].cost[0] 	= 170;    
    operator_table[ACOS].cost[0] 	= 190;
    operator_table[ASIN].cost[0] 	= 180;
    operator_table[ATAN].cost[0] 	= 190;
    operator_table[COSH].cost[0] 	= 240;
    operator_table[SINH].cost[0] 	= 240;
    operator_table[TANH].cost[0] 	= 240;
    operator_table[ACOSH].cost[0] 	= 210;
    operator_table[ASINH].cost[0] 	= 220;
    operator_table[ATANH].cost[0] 	= 150;
    operator_table[SQRT].cost[0] 	= 90;
    operator_table[NAME].cost[0] 	= 0;
	
	operator_table[COMMA].isfunction 	= false;
    operator_table[EQUAL].isfunction 	= false;
	operator_table[PLUS].isfunction 	= false; 
    operator_table[MINUS].isfunction 	= false;
    operator_table[TIMES].isfunction 	= false;
    operator_table[DIVIDE].isfunction 	= false;
    operator_table[UMINUS].isfunction 	= false;
    operator_table[POWER].isfunction 	= false;
    operator_table[EXP].isfunction 		= true;
    operator_table[LOG].isfunction 		= true;
    operator_table[LOG10].isfunction 	= true;
    operator_table[COS].isfunction 		= true;   
    operator_table[SIN].isfunction 		= true;   
    operator_table[TAN].isfunction 		= true;   
    operator_table[ACOS].isfunction 	= true;
    operator_table[ASIN].isfunction 	= true;
    operator_table[ATAN].isfunction 	= true;
    operator_table[COSH].isfunction 	= true;
    operator_table[SINH].isfunction 	= true;
    operator_table[TANH].isfunction 	= true;
    operator_table[ACOSH].isfunction 	= true;
    operator_table[ASINH].isfunction 	= true;
    operator_table[ATANH].isfunction 	= true;
    operator_table[SQRT].isfunction 	= true;
    operator_table[NAME].isfunction 	= true;
}
//------------------------------------------------------------------------------
OperatorTable::~OperatorTable()
{
	// Empty
}
//------------------------------------------------------------------------------
