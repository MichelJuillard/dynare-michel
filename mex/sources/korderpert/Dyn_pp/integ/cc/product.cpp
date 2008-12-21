/*1:*/

#include "product.h"
#include "symmetry.h"

/*2:*/

prodpit::prodpit()
:prodq(NULL),level(0),npoints(0),jseq(NULL),
end_flag(true),sig(NULL),p(NULL)
{
}

/*:2*/
;
/*3:*/

prodpit::prodpit(const ProductQuadrature&q,int j0,int l)
:prodq(&q),level(l),npoints(q.uquad.numPoints(l)),jseq(new IntSequence(q.dimen(),0)),
end_flag(false),sig(new ParameterSignal(q.dimen())),p(new Vector(q.dimen()))
{
if(j0<npoints){
(*jseq)[0]= j0;
setPointAndWeight();
}else{
end_flag= true;
}
}

/*:3*/
;
/*4:*/

prodpit::prodpit(const prodpit&ppit)
:prodq(ppit.prodq),level(ppit.level),npoints(ppit.npoints),
end_flag(ppit.end_flag),w(ppit.w)
{
if(ppit.jseq)
jseq= new IntSequence(*(ppit.jseq));
else
jseq= NULL;
if(ppit.sig)
sig= new ParameterSignal(*(ppit.sig));
else
sig= NULL;
if(ppit.p)
p= new Vector(*(ppit.p));
else
p= NULL;
}

/*:4*/
;
/*5:*/

prodpit::~prodpit()
{
if(jseq)
delete jseq;
if(sig)
delete sig;
if(p)
delete p;
}

/*:5*/
;
/*6:*/

bool prodpit::operator==(const prodpit&ppit)const
{
bool ret= true;
ret= ret&prodq==ppit.prodq;
ret= ret&end_flag==ppit.end_flag;
ret= ret&((jseq==NULL&&ppit.jseq==NULL)||
(jseq!=NULL&&ppit.jseq!=NULL&&*jseq==*(ppit.jseq)));
return ret;
}

/*:6*/
;
/*7:*/

const prodpit&prodpit::operator= (const prodpit&ppit)
{
prodq= ppit.prodq;
end_flag= ppit.end_flag;
w= ppit.w;

if(jseq)
delete jseq;
if(sig)
delete sig;
if(p)
delete p;

if(ppit.jseq)
jseq= new IntSequence(*(ppit.jseq));
else
jseq= NULL;
if(ppit.sig)
sig= new ParameterSignal(*(ppit.sig));
else
sig= NULL;
if(ppit.p)
p= new Vector(*(ppit.p));
else
p= NULL;

return*this;
}

/*:7*/
;
/*8:*/

prodpit&prodpit::operator++()
{

int i= prodq->dimen()-1;
(*jseq)[i]++;
while(i>=0&&(*jseq)[i]==npoints){
(*jseq)[i]= 0;
i--;
if(i>=0)
(*jseq)[i]++;
}
sig->signalAfter(std::max(i,0));

if(i==-1)
end_flag= true;

if(!end_flag)
setPointAndWeight();

return*this;
}


/*:8*/
;
/*9:*/

void prodpit::setPointAndWeight()
{


w= 1.0;
for(int i= 0;i<prodq->dimen();i++){
(*p)[i]= (prodq->uquad).point(level,(*jseq)[i]);
w*= (prodq->uquad).weight(level,(*jseq)[i]);
}
}

/*:9*/
;
/*10:*/

void prodpit::print()const
{
printf("j=[");
for(int i= 0;i<prodq->dimen();i++)
printf("%2d ",(*jseq)[i]);
printf("] %+4.3f*(",w);
for(int i= 0;i<prodq->dimen()-1;i++)
printf("%+4.3f ",(*p)[i]);
printf("%+4.3f)\n",(*p)[prodq->dimen()-1]);
}

/*:10*/
;
/*11:*/

ProductQuadrature::ProductQuadrature(int d,const OneDQuadrature&uq)
:QuadratureImpl<prodpit> (d),uquad(uq)
{

}

/*:11*/
;
/*12:*/

prodpit ProductQuadrature::begin(int ti,int tn,int l)const
{


int npoints= uquad.numPoints(l);
return prodpit(*this,ti*npoints/tn,l);
}

/*:12*/
;
/*13:*/

void ProductQuadrature::designLevelForEvals(int max_evals,int&lev,int&evals)const
{
int last_evals;
evals= 1;
lev= 1;
do{
lev++;
last_evals= evals;
evals= numEvals(lev);
}while(lev<uquad.numLevels()-2&&evals<max_evals);
lev--;
evals= last_evals;

}

/*:13*/
;

/*:1*/
