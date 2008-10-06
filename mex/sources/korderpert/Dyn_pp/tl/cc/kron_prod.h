/*1:*/
#line 30 "./kron_prod.hweb"


#ifndef KRON_PROD_H
#define KRON_PROD_H

#include "twod_matrix.h"
#include "permutation.h"
#include "int_sequence.h"

class KronProdAll;
class KronProdAllOptim;
class KronProdIA;
class KronProdIAI;
class KronProdAI;

/*2:*/
#line 58 "./kron_prod.hweb"

class KronProdDimens{
friend class KronProdAll;
friend class KronProdAllOptim;
friend class KronProdIA;
friend class KronProdIAI;
friend class KronProdAI;
private:
IntSequence rows;
IntSequence cols;
public:
/*3:*/
#line 86 "./kron_prod.hweb"

KronProdDimens(int dim)
:rows(dim,0),cols(dim,0){}
KronProdDimens(const KronProdDimens&kd)
:rows(kd.rows),cols(kd.cols){}
KronProdDimens(const KronProdDimens&kd,int i);

/*:3*/
#line 69 "./kron_prod.hweb"
;
/*4:*/
#line 94 "./kron_prod.hweb"

const KronProdDimens&operator= (const KronProdDimens&kd)
{rows= kd.rows;cols= kd.cols;return*this;}
bool operator==(const KronProdDimens&kd)const
{return rows==kd.rows&&cols==kd.cols;}

/*:4*/
#line 70 "./kron_prod.hweb"
;
/*5:*/
#line 101 "./kron_prod.hweb"

int dimen()const
{return rows.size();}
void setRC(int i,int r,int c)
{rows[i]= r;cols[i]= c;}
void getRC(int i,int&r,int&c)const
{r= rows[i];c= cols[i];}
void getRC(int&r,int&c)const
{r= rows.mult();c= cols.mult();}
int nrows()const
{return rows.mult();}
int ncols()const
{return cols.mult();}
int nrows(int i)const
{return rows[i];}
int ncols(int i)const
{return cols[i];}

/*:5*/
#line 71 "./kron_prod.hweb"
;
};

/*:2*/
#line 45 "./kron_prod.hweb"
;
/*6:*/
#line 130 "./kron_prod.hweb"

class KronProd{
protected:
KronProdDimens kpd;
public:
KronProd(int dim)
:kpd(dim){}
KronProd(const KronProdDimens&kd)
:kpd(kd){}
KronProd(const KronProd&kp)
:kpd(kp.kpd){}
virtual~KronProd(){}

int dimen()const
{return kpd.dimen();}

virtual void mult(const ConstTwoDMatrix&in,TwoDMatrix&out)const= 0;
void mult(const TwoDMatrix&in,TwoDMatrix&out)const
{mult(ConstTwoDMatrix(in),out);}

void checkDimForMult(const ConstTwoDMatrix&in,const TwoDMatrix&out)const;
void checkDimForMult(const TwoDMatrix&in,const TwoDMatrix&out)const
{checkDimForMult(ConstTwoDMatrix(in),out);}

static void kronMult(const ConstVector&v1,const ConstVector&v2,
Vector&res);

int nrows()const
{return kpd.nrows();}
int ncols()const
{return kpd.ncols();}
int nrows(int i)const
{return kpd.nrows(i);}
int ncols(int i)const
{return kpd.ncols(i);}
};

/*:6*/
#line 46 "./kron_prod.hweb"
;
/*7:*/
#line 183 "./kron_prod.hweb"

class KronProdAll:public KronProd{
friend class KronProdIA;
friend class KronProdIAI;
friend class KronProdAI;
protected:
const TwoDMatrix**const matlist;
public:
KronProdAll(int dim)
:KronProd(dim),matlist(new const TwoDMatrix*[dim]){}
virtual~KronProdAll()
{delete[]matlist;}
void setMat(int i,const TwoDMatrix&m);
void setUnit(int i,int n);
const TwoDMatrix&getMat(int i)const
{return*(matlist[i]);}

void mult(const ConstTwoDMatrix&in,TwoDMatrix&out)const;
Vector*multRows(const IntSequence&irows)const;
private:
bool isUnit()const;
};

/*:7*/
#line 47 "./kron_prod.hweb"
;
/*8:*/
#line 235 "./kron_prod.hweb"

class KronProdAllOptim:public KronProdAll{
protected:
Permutation oper;
public:
KronProdAllOptim(int dim)
:KronProdAll(dim),oper(dim){}
void optimizeOrder();
const Permutation&getPer()const
{return oper;}
};

/*:8*/
#line 48 "./kron_prod.hweb"
;
/*9:*/
#line 250 "./kron_prod.hweb"

class KronProdIA:public KronProd{
friend class KronProdAll;
const TwoDMatrix&mat;
public:
KronProdIA(const KronProdAll&kpa)
:KronProd(KronProdDimens(kpa.kpd,kpa.dimen()-1)),
mat(kpa.getMat(kpa.dimen()-1))
{}
void mult(const ConstTwoDMatrix&in,TwoDMatrix&out)const;
};

/*:9*/
#line 49 "./kron_prod.hweb"
;
/*10:*/
#line 265 "./kron_prod.hweb"

class KronProdAI:public KronProd{
friend class KronProdIAI;
friend class KronProdAll;
const TwoDMatrix&mat;
public:
KronProdAI(const KronProdAll&kpa)
:KronProd(KronProdDimens(kpa.kpd,0)),
mat(kpa.getMat(0))
{}
KronProdAI(const KronProdIAI&kpiai);

void mult(const ConstTwoDMatrix&in,TwoDMatrix&out)const;
};

/*:10*/
#line 50 "./kron_prod.hweb"
;
/*11:*/
#line 282 "./kron_prod.hweb"

class KronProdIAI:public KronProd{
friend class KronProdAI;
friend class KronProdAll;
const TwoDMatrix&mat;
public:
KronProdIAI(const KronProdAll&kpa,int i)
:KronProd(KronProdDimens(kpa.kpd,i)),
mat(kpa.getMat(i))
{}
void mult(const ConstTwoDMatrix&in,TwoDMatrix&out)const;
};


/*:11*/
#line 51 "./kron_prod.hweb"
;

#endif

/*:1*/
