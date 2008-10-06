/*1:*/
#line 6 "./korder.cweb"


#include "kord_exception.h"
#include "korder.h"

#include "cpplapack.h"

/*2:*/
#line 25 "./korder.cweb"

PLUMatrix::PLUMatrix(const PLUMatrix&plu)
:TwoDMatrix(plu),inv(plu.inv),ipiv(new int[nrows()])
{
memcpy(ipiv,plu.ipiv,nrows()*sizeof(int));
}


/*:2*/
#line 13 "./korder.cweb"
;
/*3:*/
#line 37 "./korder.cweb"

void PLUMatrix::calcPLU()
{
int info;
int rows= nrows();
inv= (const Vector&)getData();
LAPACK_dgetrf(&rows,&rows,inv.base(),&rows,ipiv,&info);
}

/*:3*/
#line 14 "./korder.cweb"
;
/*4:*/
#line 48 "./korder.cweb"

void PLUMatrix::multInv(TwoDMatrix&m)const
{
KORD_RAISE_IF(m.nrows()!=ncols(),
"The matrix is not square in PLUMatrix::multInv");
int info;
int mcols= m.ncols();
int mrows= m.nrows();
double*mbase= m.getData().base();
LAPACK_dgetrs("N",&mrows,&mcols,inv.base(),&mrows,ipiv,
mbase,&mrows,&info);
KORD_RAISE_IF(info!=0,
"Info!=0 in PLUMatrix::multInv");
}

/*:4*/
#line 15 "./korder.cweb"
;
/*5:*/
#line 69 "./korder.cweb"

MatrixA::MatrixA(const FSSparseTensor&f,const IntSequence&ss,
const TwoDMatrix&gy,const PartitionY&ypart)
:PLUMatrix(ypart.ny())
{
zeros();

IntSequence c(1);c[0]= 1;
FGSTensor f_y(f,ss,c,TensorDimens(ss,c));
add(1.0,f_y);

ConstTwoDMatrix gss_ys(ypart.nstat+ypart.npred,ypart.nyss(),gy);
c[0]= 0;
FGSTensor f_yss(f,ss,c,TensorDimens(ss,c));
TwoDMatrix sub(*this,ypart.nstat,ypart.nys());
sub.multAndAdd(ConstTwoDMatrix(f_yss),gss_ys);

calcPLU();
}

/*:5*/
#line 16 "./korder.cweb"
;
/*6:*/
#line 97 "./korder.cweb"

MatrixS::MatrixS(const FSSparseTensor&f,const IntSequence&ss,
const TwoDMatrix&gy,const PartitionY&ypart)
:PLUMatrix(ypart.ny())
{
zeros();

IntSequence c(1);c[0]= 1;
FGSTensor f_y(f,ss,c,TensorDimens(ss,c));
add(1.0,f_y);

ConstTwoDMatrix gss_ys(ypart.nstat+ypart.npred,ypart.nyss(),gy);
c[0]= 0;
FGSTensor f_yss(f,ss,c,TensorDimens(ss,c));
TwoDMatrix sub(*this,ypart.nstat,ypart.nys());
sub.multAndAdd(ConstTwoDMatrix(f_yss),gss_ys);

TwoDMatrix sub2(*this,ypart.nstat+ypart.npred,ypart.nyss());
sub2.add(1.0,f_yss);

calcPLU();
}


/*:6*/
#line 17 "./korder.cweb"
;
/*13:*/
#line 281 "./korder.cweb"

template<> ctraits<KOrder::unfold> ::Tg&KOrder::g<KOrder::unfold> ()
{return _ug;}
template<> const ctraits<KOrder::unfold> ::Tg&KOrder::g<KOrder::unfold> ()const
{return _ug;}
template<> ctraits<KOrder::fold> ::Tg&KOrder::g<KOrder::fold> ()
{return _fg;}
template<> const ctraits<KOrder::fold> ::Tg&KOrder::g<KOrder::fold> ()const
{return _fg;}
template<> ctraits<KOrder::unfold> ::Tgs&KOrder::gs<KOrder::unfold> ()
{return _ugs;}
template<> const ctraits<KOrder::unfold> ::Tgs&KOrder::gs<KOrder::unfold> ()const
{return _ugs;}
template<> ctraits<KOrder::fold> ::Tgs&KOrder::gs<KOrder::fold> ()
{return _fgs;}
template<> const ctraits<KOrder::fold> ::Tgs&KOrder::gs<KOrder::fold> ()const
{return _fgs;}
template<> ctraits<KOrder::unfold> ::Tgss&KOrder::gss<KOrder::unfold> ()
{return _ugss;}
template<> const ctraits<KOrder::unfold> ::Tgss&KOrder::gss<KOrder::unfold> ()const
{return _ugss;}
template<> ctraits<KOrder::fold> ::Tgss&KOrder::gss<KOrder::fold> ()
{return _fgss;}
template<> const ctraits<KOrder::fold> ::Tgss&KOrder::gss<KOrder::fold> ()const
{return _fgss;}
template<> ctraits<KOrder::unfold> ::TG&KOrder::G<KOrder::unfold> ()
{return _uG;}
template<> const ctraits<KOrder::unfold> ::TG&KOrder::G<KOrder::unfold> ()const
{return _uG;}
template<> ctraits<KOrder::fold> ::TG&KOrder::G<KOrder::fold> ()
{return _fG;}
template<> const ctraits<KOrder::fold> ::TG&KOrder::G<KOrder::fold> ()const
{return _fG;}
template<> ctraits<KOrder::unfold> ::TZstack&KOrder::Zstack<KOrder::unfold> ()
{return _uZstack;}
template<> const ctraits<KOrder::unfold> ::TZstack&KOrder::Zstack<KOrder::unfold> ()const
{return _uZstack;}
template<> ctraits<KOrder::fold> ::TZstack&KOrder::Zstack<KOrder::fold> ()
{return _fZstack;}
template<> const ctraits<KOrder::fold> ::TZstack&KOrder::Zstack<KOrder::fold> ()const
{return _fZstack;}
template<> ctraits<KOrder::unfold> ::TGstack&KOrder::Gstack<KOrder::unfold> ()
{return _uGstack;}
template<> const ctraits<KOrder::unfold> ::TGstack&KOrder::Gstack<KOrder::unfold> ()const
{return _uGstack;}
template<> ctraits<KOrder::fold> ::TGstack&KOrder::Gstack<KOrder::fold> ()
{return _fGstack;}
template<> const ctraits<KOrder::fold> ::TGstack&KOrder::Gstack<KOrder::fold> ()const
{return _fGstack;}
template<> ctraits<KOrder::unfold> ::Tm&KOrder::m<KOrder::unfold> ()
{return _um;}
template<> const ctraits<KOrder::unfold> ::Tm&KOrder::m<KOrder::unfold> ()const
{return _um;}
template<> ctraits<KOrder::fold> ::Tm&KOrder::m<KOrder::fold> ()
{return _fm;}
template<> const ctraits<KOrder::fold> ::Tm&KOrder::m<KOrder::fold> ()const
{return _fm;}


/*:13*/
#line 18 "./korder.cweb"
;
/*10:*/
#line 211 "./korder.cweb"

template<> 
void KOrder::sylvesterSolve<KOrder::unfold> (ctraits<unfold> ::Ttensor&der)const
{
JournalRecordPair pa(journal);
pa<<"Sylvester equation for dimension = "<<der.getSym()[0]<<endrec;
if(ypart.nys()> 0&&ypart.nyss()> 0){
KORD_RAISE_IF(!der.isFinite(),
"RHS of Sylverster is not finite");
TwoDMatrix gs_y(*(gs<unfold> ().get(Symmetry(1,0,0,0))));
GeneralSylvester sylv(der.getSym()[0],ny,ypart.nys(),
ypart.nstat+ypart.npred,
matA.getData().base(),matB.getData().base(),
gs_y.getData().base(),der.getData().base());
sylv.solve();
}else if(ypart.nys()> 0&&ypart.nyss()==0){
matA.multInv(der);
}
}

/*:10*/
#line 19 "./korder.cweb"
;
/*11:*/
#line 235 "./korder.cweb"

template<> 
void KOrder::sylvesterSolve<KOrder::fold> (ctraits<fold> ::Ttensor&der)const
{
ctraits<unfold> ::Ttensor tmp(der);
sylvesterSolve<unfold> (tmp);
ctraits<fold> ::Ttensor ftmp(tmp);
der.getData()= (const Vector&)(ftmp.getData());
}

/*:11*/
#line 20 "./korder.cweb"
;
/*12:*/
#line 246 "./korder.cweb"

void KOrder::switchToFolded()
{
JournalRecordPair pa(journal);
pa<<"Switching from unfolded to folded"<<endrec;

int maxdim= g<unfold> ().getMaxDim();
for(int dim= 1;dim<=maxdim;dim++){
SymmetrySet ss(dim,4);
for(symiterator si(ss);!si.isEnd();++si){
if((*si)[2]==0&&g<unfold> ().check(*si)){
FGSTensor*ft= new FGSTensor(*(g<unfold> ().get(*si)));
insertDerivative<fold> (ft);
if(dim> 1){
gss<unfold> ().remove(*si);
gs<unfold> ().remove(*si);
g<unfold> ().remove(*si);
}
}
if(G<unfold> ().check(*si)){
FGSTensor*ft= new FGSTensor(*(G<unfold> ().get(*si)));
G<fold> ().insert(ft);
if(dim> 1){
G<fold> ().remove(*si);
}
}
}
}
}



/*:12*/
#line 21 "./korder.cweb"
;
/*7:*/
#line 133 "./korder.cweb"

KOrder::KOrder(int num_stat,int num_pred,int num_both,int num_forw,
const TensorContainer<FSSparseTensor> &fcont,
const TwoDMatrix&gy,const TwoDMatrix&gu,const TwoDMatrix&v,
Journal&jr)
:ypart(num_stat,num_pred,num_both,num_forw),
ny(ypart.ny()),nu(gu.ncols()),maxk(fcont.getMaxDim()),
nvs(4),
_ug(4),_fg(4),_ugs(4),_fgs(4),_ugss(4),_fgss(4),
_uG(4),_fG(4),
_uZstack(&_uG,ypart.nyss(),&_ug,ny,ypart.nys(),nu),
_fZstack(&_fG,ypart.nyss(),&_fg,ny,ypart.nys(),nu),
_uGstack(&_ugs,ypart.nys(),nu),
_fGstack(&_fgs,ypart.nys(),nu),
_um(maxk,v),_fm(_um),f(fcont),
matA(*(f.get(Symmetry(1))),_uZstack.getStackSizes(),gy,ypart),
matS(*(f.get(Symmetry(1))),_uZstack.getStackSizes(),gy,ypart),
matB(*(f.get(Symmetry(1))),_uZstack.getStackSizes()),
journal(jr)
{
KORD_RAISE_IF(gy.ncols()!=ypart.nys(),
"Wrong number of columns in gy in KOrder constructor");
KORD_RAISE_IF(v.ncols()!=nu,
"Wrong number of columns of Vcov in KOrder constructor");
KORD_RAISE_IF(nu!=v.nrows(),
"Wrong number of rows of Vcov in KOrder constructor");
KORD_RAISE_IF(maxk<2,
"Order of approximation must be at least 2 in KOrder constructor");
KORD_RAISE_IF(gy.nrows()!=ypart.ny(),
"Wrong number of rows in gy in KOrder constructor");
KORD_RAISE_IF(gu.nrows()!=ypart.ny(),
"Wrong number of rows in gu in KOrder constuctor");
KORD_RAISE_IF(gu.ncols()!=nu,
"Wrong number of columns in gu in KOrder constuctor");


nvs[0]= ypart.nys();nvs[1]= nu;nvs[2]= nu;nvs[3]= 1;

/*8:*/
#line 178 "./korder.cweb"

UGSTensor*tgy= new UGSTensor(ny,TensorDimens(Symmetry(1,0,0,0),nvs));
tgy->getData()= gy.getData();
insertDerivative<unfold> (tgy);
UGSTensor*tgu= new UGSTensor(ny,TensorDimens(Symmetry(0,1,0,0),nvs));
tgu->getData()= gu.getData();
insertDerivative<unfold> (tgu);

/*:8*/
#line 171 "./korder.cweb"
;
/*9:*/
#line 187 "./korder.cweb"

UGSTensor*tGy= faaDiBrunoG<unfold> (Symmetry(1,0,0,0));
G<unfold> ().insert(tGy);
UGSTensor*tGu= faaDiBrunoG<unfold> (Symmetry(0,1,0,0));
G<unfold> ().insert(tGu);
UGSTensor*tGup= faaDiBrunoG<unfold> (Symmetry(0,0,1,0));
G<unfold> ().insert(tGup);



/*:9*/
#line 172 "./korder.cweb"
;
}

/*:7*/
#line 22 "./korder.cweb"
;

/*:1*/
