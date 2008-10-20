/*1:*/

#ifndef TWOD_MATRIX_H
#define TWOD_MATRIX_H

#include "GeneralMatrix.h"

#include <stdio.h> 

class TwoDMatrix;
/*2:*/

class ConstTwoDMatrix:public ConstGeneralMatrix{
public:
	ConstTwoDMatrix(int m,int n,const double*d)
		:ConstGeneralMatrix(d,m,n){}
	ConstTwoDMatrix(const TwoDMatrix&m);
	ConstTwoDMatrix(const TwoDMatrix&m,int first_col,int num);
	ConstTwoDMatrix(const ConstTwoDMatrix&m,int first_col,int num);
	ConstTwoDMatrix(int first_row,int num,const TwoDMatrix&m);
	ConstTwoDMatrix(int first_row,int num,const ConstTwoDMatrix&m);
	ConstTwoDMatrix(const ConstTwoDMatrix&m,int first_row,int first_col,int rows,int cols)
		:ConstGeneralMatrix(m,first_row,first_col,rows,cols){}
	virtual~ConstTwoDMatrix(){}
	
	int nrows()const
	{return numRows();}
	int ncols()const
	{return numCols();}
	void writeMat4(FILE*fd,const char*vname)const;
};

/*:2*/
;
/*3:*/

class TwoDMatrix:public GeneralMatrix{
public:
	TwoDMatrix(int r,int c)
		:GeneralMatrix(r,c){}
	TwoDMatrix(int r,int c,double*d)
		:GeneralMatrix(d,r,c){}
	TwoDMatrix(int r,int c,const double*d)
		:GeneralMatrix(d,r,c){}
	TwoDMatrix(const GeneralMatrix&m)
		:GeneralMatrix(m){}
	TwoDMatrix(const GeneralMatrix&m,char*dummy)
		:GeneralMatrix(m,dummy){}
	TwoDMatrix(const TwoDMatrix&m,int first_col,int num)
		:GeneralMatrix(m,0,first_col,m.numRows(),num){}
	TwoDMatrix(TwoDMatrix&m,int first_col,int num)
		:GeneralMatrix(m,0,first_col,m.numRows(),num){}
	TwoDMatrix(int first_row,int num,const TwoDMatrix&m)
		:GeneralMatrix(m,first_row,0,num,m.ncols()){}
	TwoDMatrix(int first_row,int num,TwoDMatrix&m)
		:GeneralMatrix(m,first_row,0,num,m.ncols()){}
	TwoDMatrix(TwoDMatrix&m,int first_row,int first_col,int rows,int cols)
		:GeneralMatrix(m,first_row,first_col,rows,cols){}
	TwoDMatrix(const TwoDMatrix&m,int first_row,int first_col,int rows,int cols)
		:GeneralMatrix(m,first_row,first_col,rows,cols){}
	TwoDMatrix(const ConstTwoDMatrix&a,const ConstTwoDMatrix&b)
		:GeneralMatrix(a,b){}
	virtual~TwoDMatrix(){}
	
	int nrows()const
	{return numRows();}
	int ncols()const
	{return numCols();}
	
	/*4:*/
	
	void copyRow(int from,int to);
	void copyRow(const ConstTwoDMatrix&m,int from,int to);
	void copyRow(const TwoDMatrix&m,int from,int to)
	{copyRow(ConstTwoDMatrix(m),from,to);}
	void addRow(const ConstTwoDMatrix&m,int from,int to)
	{addRow(1.0,m,from,to);}
	void addRow(const TwoDMatrix&m,int from,int to)
	{addRow(1.0,ConstTwoDMatrix(m),from,to);}
	void addRow(double d,const ConstTwoDMatrix&m,int from,int to);
	void addRow(double d,const TwoDMatrix&m,int from,int to)
	{addRow(d,ConstTwoDMatrix(m),from,to);}
	
	
	/*:4*/
	;
	/*5:*/
	
	void copyColumn(int from,int to);
	void copyColumn(const ConstTwoDMatrix&m,int from,int to);
	void copyColumn(const TwoDMatrix&m,int from,int to)
	{copyColumn(ConstTwoDMatrix(m),from,to);}
	void addColumn(const ConstTwoDMatrix&m,int from,int to)
	{addColumn(1.0,ConstTwoDMatrix(m),from,to);}
	void addColumn(const TwoDMatrix&m,int from,int to)
	{addColumn(1.0,ConstTwoDMatrix(m),from,to);}
	void addColumn(double d,const ConstTwoDMatrix&m,int from,int to);
	void addColumn(double d,const TwoDMatrix&m,int from,int to)
	{addColumn(d,ConstTwoDMatrix(m),from,to);}
	
	/*:5*/
	;
	void save(const char*fname)const;
	void writeMat4(FILE*fd,const char*vname)const
	{ConstTwoDMatrix(*this).writeMat4(fd,vname);}
};

/*:3*/
;
/*6:*/

class Mat4Header{
	int type;
	int rows;
	int cols;
	int imagf;
	int namelen;
	const char*vname;
public:
	Mat4Header(const ConstTwoDMatrix&m,const char*vname);
	Mat4Header(const ConstTwoDMatrix&m,const char*vname,const char*dummy);
	void write(FILE*fd)const;
};



/*:6*/
;

#endif


/*:1*/
