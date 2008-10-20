/*1:*/

#include "twod_matrix.h"
#include "tl_exception.h"


/*2:*/

ConstTwoDMatrix::ConstTwoDMatrix(const TwoDMatrix&m)
:ConstGeneralMatrix(m){}

ConstTwoDMatrix::ConstTwoDMatrix(const TwoDMatrix&m,int first_col,int num)
:ConstGeneralMatrix(m,0,first_col,m.nrows(),num){}

ConstTwoDMatrix::ConstTwoDMatrix(const ConstTwoDMatrix&m,int first_col,int num)
:ConstGeneralMatrix(m,0,first_col,m.nrows(),num){}

ConstTwoDMatrix::ConstTwoDMatrix(int first_row,int num,const TwoDMatrix&m)
:ConstGeneralMatrix(m,first_row,0,num,m.ncols()){}

ConstTwoDMatrix::ConstTwoDMatrix(int first_row,int num,const ConstTwoDMatrix&m)
:ConstGeneralMatrix(m,first_row,0,num,m.ncols()){}

/*:2*/
;
/*3:*/

void ConstTwoDMatrix::writeMat4(FILE*fd,const char*vname)const
{
	Mat4Header header(*this,vname);
	header.write(fd);
	for(int j= 0;j<ncols();j++)
		for(int i= 0;i<nrows();i++)
			fwrite(&(get(i,j)),sizeof(double),1,fd);
}

/*:3*/
;
/*4:*/

void TwoDMatrix::copyRow(int from,int to)
{
	if(from!=to)
		copyRow(ConstTwoDMatrix(*this),from,to);
}

void TwoDMatrix::copyRow(const ConstTwoDMatrix&m,int from,int to)
{
	ConstVector fr_row(from,m);
	Vector to_row(to,*this);
	to_row= fr_row;
}

void TwoDMatrix::addRow(double d,const ConstTwoDMatrix&m,int from,int to)
{
	ConstVector fr_row(from,m);
	Vector to_row(to,*this);
	to_row.add(d,fr_row);
}


/*:4*/
;
/*5:*/

void TwoDMatrix::copyColumn(int from,int to)
{
	if(from!=to)
		copyColumn(ConstTwoDMatrix(*this),from,to);
}

void TwoDMatrix::copyColumn(const ConstTwoDMatrix&m,int from,int to)
{
	ConstVector fr_col(m,from);
	Vector to_col(*this,to);
	to_col= fr_col;
}

void TwoDMatrix::addColumn(double d,const ConstTwoDMatrix&m,int from,int to)
{
	ConstVector fr_col(m,from);
	Vector to_col(*this,to);
	to_col.add(d,fr_col);
}

/*:5*/
;
/*6:*/

void TwoDMatrix::save(const char*fname)const
{
	FILE*fd;
	if(NULL==(fd= fopen(fname,"w"))){
		TL_RAISE("Cannot open file for writing in TwoDMatrix::save");
	}
	for(int row= 0;row<nrows();row++){
		for(int col= 0;col<ncols();col++)
			fprintf(fd," %20.10g",get(row,col));
		fprintf(fd,"\n");
	}
	fclose(fd);
}

/*:6*/
;
/*7:*/

Mat4Header::Mat4Header(const ConstTwoDMatrix&m,const char*vn)
:type(0),rows(m.nrows()),cols(m.ncols()),imagf(0),namelen(strlen(vn)+1),
vname(vn)
{}


/*:7*/
;
/*8:*/

Mat4Header::Mat4Header(const ConstTwoDMatrix&m,const char*vn,const char*dummy)
:type(1),rows(m.nrows()),cols(m.ncols()),imagf(0),namelen(strlen(vn)+1),
vname(vn)
{}


/*:8*/
;
/*9:*/

void Mat4Header::write(FILE*fd)const
{
	fwrite(&type,sizeof(int),1,fd);
	fwrite(&rows,sizeof(int),1,fd);
	fwrite(&cols,sizeof(int),1,fd);
	fwrite(&imagf,sizeof(int),1,fd);
	fwrite(&namelen,sizeof(int),1,fd);
	fwrite(vname,1,namelen,fd);
}


/*:9*/
;

/*:1*/
