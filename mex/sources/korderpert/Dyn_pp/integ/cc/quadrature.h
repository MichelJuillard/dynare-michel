/*1:*/

#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "vector_function.h"
#include "int_sequence.h"
#include "sthread.h"

/*2:*/

class OneDQuadrature{
public:
virtual~OneDQuadrature(){}
virtual int numLevels()const= 0;
virtual int numPoints(int level)const= 0;
virtual double point(int level,int i)const= 0;
virtual double weight(int lelel,int i)const= 0;
};

/*:2*/
;
/*3:*/

class Quadrature{
protected:
int dim;
public:
Quadrature(int d):dim(d){}
virtual~Quadrature(){}
int dimen()const
{return dim;}
virtual void integrate(const VectorFunction&func,int level,
int tn,Vector&out)const= 0;
virtual void integrate(VectorFunctionSet&fs,int level,Vector&out)const= 0;
virtual int numEvals(int level)const= 0;
};

/*:3*/
;
/*4:*/

template<typename _Tpit> 
class QuadratureImpl;

template<typename _Tpit> 
class IntegrationWorker:public THREAD{
const QuadratureImpl<_Tpit> &quad;
VectorFunction&func;
int level;
int ti;
int tn;
Vector&outvec;
public:
IntegrationWorker(const QuadratureImpl<_Tpit> &q,VectorFunction&f,int l,
int tii,int tnn,Vector&out)
:quad(q),func(f),level(l),ti(tii),tn(tnn),outvec(out){}
/*5:*/

void operator()(){
_Tpit beg= quad.begin(ti,tn,level);
_Tpit end= quad.begin(ti+1,tn,level);
Vector tmpall(outvec.length());
tmpall.zeros();
Vector tmp(outvec.length());



for(_Tpit run= beg;run!=end;++run){
func.eval(run.point(),run.signal(),tmp);
tmpall.add(run.weight(),tmp);
}

{
SYNCHRO syn(&outvec,"IntegrationWorker");
outvec.add(1.0,tmpall);
}
}


/*:5*/
;
};


/*:4*/
;
/*6:*/

template<typename _Tpit> 
class QuadratureImpl:public Quadrature{
friend class IntegrationWorker<_Tpit> ;
public:
QuadratureImpl(int d):Quadrature(d){}
/*7:*/

void integrate(VectorFunctionSet&fs,int level,Vector&out)const{


out.zeros();
THREAD_GROUP gr;
for(int ti= 0;ti<fs.getNum();ti++){
gr.insert(new IntegrationWorker<_Tpit> (*this,fs.getFunc(ti),
level,ti,fs.getNum(),out));
}
gr.run();
}


/*:7*/
;
void integrate(const VectorFunction&func,
int level,int tn,Vector&out)const{
VectorFunctionSet fs(func,tn);
integrate(fs,level,out);
}
/*8:*/

void savePoints(const char*fname,int level)const
{
FILE*fd;
if(NULL==(fd= fopen(fname,"w"))){

fprintf(stderr,"Cannot open file %s for writing.\n",fname);
exit(1);
}
_Tpit beg= begin(0,1,level);
_Tpit end= begin(1,1,level);
for(_Tpit run= beg;run!=end;++run){
fprintf(fd,"%16.12g",run.weight());
for(int i= 0;i<dimen();i++)
fprintf(fd,"\t%16.12g",run.point()[i]);
fprintf(fd,"\n");
}
fclose(fd);
}


/*:8*/
;
_Tpit start(int level)const
{return begin(0,1,level);}
_Tpit end(int level)const
{return begin(1,1,level);}
protected:
virtual _Tpit begin(int ti,int tn,int level)const= 0;
};

/*:6*/
;
/*9:*/

class OneDPrecalcQuadrature:public OneDQuadrature{
int num_levels;
const int*num_points;
const double*weights;
const double*points;
IntSequence offsets;
public:
OneDPrecalcQuadrature(int nlevels,const int*npoints,
const double*wts,const double*pts)
:num_levels(nlevels),num_points(npoints),
weights(wts),points(pts),offsets(num_levels)
{calcOffsets();}
virtual~OneDPrecalcQuadrature(){}
int numLevels()const
{return num_levels;}
int numPoints(int level)const
{return num_points[level-1];}
double point(int level,int i)const
{return points[offsets[level-1]+i];}
double weight(int level,int i)const
{return weights[offsets[level-1]+i];}
protected:
void calcOffsets();
};

/*:9*/
;
/*10:*/

class GaussHermite:public OneDPrecalcQuadrature{
public:
GaussHermite();
};

/*:10*/
;
/*11:*/

class GaussLegendre:public OneDPrecalcQuadrature{
public:
GaussLegendre();
};

/*:11*/
;
/*12:*/

class NormalICDF{
public:
static double get(double x);
};

/*:12*/
;

#endif

/*:1*/
