/*1:*/

#include "smolyak.h"
#include "symmetry.h"

/*2:*/

smolpit::smolpit()
:smolq(NULL),isummand(0),jseq(NULL),sig(NULL),p(NULL)
{
}

/*:2*/
;
/*3:*/

smolpit::smolpit(const SmolyakQuadrature&q,unsigned int isum)
:smolq(&q),isummand(isum),jseq(new IntSequence(q.dimen(),0)),
sig(new ParameterSignal(q.dimen())),p(new Vector(q.dimen()))
{
if(isummand<q.numSummands()){
setPointAndWeight();
}
}

/*:3*/
;
/*4:*/

smolpit::smolpit(const smolpit&spit)
:smolq(spit.smolq),isummand(spit.isummand),w(spit.w)
{
if(spit.jseq)
jseq= new IntSequence(*(spit.jseq));
else
jseq= NULL;
if(spit.sig)
sig= new ParameterSignal(*(spit.sig));
else
sig= NULL;
if(spit.p)
p= new Vector(*(spit.p));
else
p= NULL;
}

/*:4*/
;
/*5:*/

smolpit::~smolpit()
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

bool smolpit::operator==(const smolpit&spit)const
{
bool ret= true;
ret= ret&smolq==spit.smolq;
ret= ret&isummand==spit.isummand;
ret= ret&((jseq==NULL&&spit.jseq==NULL)||
(jseq!=NULL&&spit.jseq!=NULL&&*jseq==*(spit.jseq)));
return ret;
}

/*:6*/
;
/*7:*/

const smolpit&smolpit::operator= (const smolpit&spit)
{
smolq= spit.smolq;
isummand= spit.isummand;
w= spit.w;

if(jseq)
delete jseq;
if(sig)
delete sig;
if(p)
delete p;

if(spit.jseq)
jseq= new IntSequence(*(spit.jseq));
else
jseq= NULL;
if(spit.sig)
sig= new ParameterSignal(*(spit.sig));
else
sig= NULL;
if(spit.p)
p= new Vector(*(spit.p));
else
p= NULL;

return*this;
}

/*:7*/
;
/*8:*/

smolpit&smolpit::operator++()
{

const IntSequence&levpts= smolq->levpoints[isummand];
int i= smolq->dimen()-1;
(*jseq)[i]++;
while(i>=0&&(*jseq)[i]==levpts[i]){
(*jseq)[i]= 0;
i--;
if(i>=0)
(*jseq)[i]++;
}
sig->signalAfter(std::max(i,0));

if(i<0)
isummand++;

if(isummand<smolq->numSummands())
setPointAndWeight();

return*this;
}


/*:8*/
;
/*9:*/

void smolpit::setPointAndWeight()
{


int l= smolq->level;
int d= smolq->dimen();
int sumk= (smolq->levels[isummand]).sum();
int m1exp= l+d-sumk-1;
w= (2*(m1exp/2)==m1exp)?1.0:-1.0;
w*= smolq->psc.noverk(d-1,sumk-l);
for(int i= 0;i<d;i++){
int ki= (smolq->levels[isummand])[i];
(*p)[i]= (smolq->uquad).point(ki,(*jseq)[i]);
w*= (smolq->uquad).weight(ki,(*jseq)[i]);
}
}

/*:9*/
;
/*10:*/

void smolpit::print()const
{
printf("isum=%-3d: [",isummand);
for(int i= 0;i<smolq->dimen();i++)
printf("%2d ",(smolq->levels[isummand])[i]);
printf("] j=[");
for(int i= 0;i<smolq->dimen();i++)
printf("%2d ",(*jseq)[i]);
printf("] %+4.3f*(",w);
for(int i= 0;i<smolq->dimen()-1;i++)
printf("%+4.3f ",(*p)[i]);
printf("%+4.3f)\n",(*p)[smolq->dimen()-1]);
}

/*:10*/
;
/*11:*/

SmolyakQuadrature::SmolyakQuadrature(int d,int l,const OneDQuadrature&uq)
:QuadratureImpl<smolpit> (d),level(l),uquad(uq),psc(d-1,d-1)
{


int cum= 0;
SymmetrySet ss(l-1,d+1);
for(symiterator si(ss);!si.isEnd();++si){
if((*si)[d]<=d-1){
IntSequence lev((const IntSequence&)*si,0,d);
lev.add(1);
levels.push_back(lev);
IntSequence levpts(d);
for(int i= 0;i<d;i++)
levpts[i]= uquad.numPoints(lev[i]);
levpoints.push_back(levpts);
cum+= levpts.mult();
cumevals.push_back(cum);
}
}
}

/*:11*/
;
/*12:*/

int SmolyakQuadrature::numEvals(int l)const
{
if(l!=level)
return calcNumEvaluations(l);
else
return cumevals[numSummands()-1];
}


/*:12*/
;
/*13:*/

smolpit SmolyakQuadrature::begin(int ti,int tn,int l)const
{

if(ti==tn)
return smolpit(*this,numSummands());

int totevals= cumevals[numSummands()-1];
int evals= (totevals*ti)/tn;
unsigned int isum= 0;
while(isum+1<numSummands()&&cumevals[isum+1]<evals)
isum++;
return smolpit(*this,isum);
}

/*:13*/
;
/*14:*/

int SmolyakQuadrature::calcNumEvaluations(int lev)const
{
int cum= 0;
SymmetrySet ss(lev-1,dim+1);
for(symiterator si(ss);!si.isEnd();++si){
if((*si)[dim]<=dim-1){
IntSequence lev((const IntSequence&)*si,0,dim);
lev.add(1);
IntSequence levpts(dim);
for(int i= 0;i<dim;i++)
levpts[i]= uquad.numPoints(lev[i]);
cum+= levpts.mult();
}
}
return cum;
}

/*:14*/
;
/*15:*/

void SmolyakQuadrature::designLevelForEvals(int max_evals,int&lev,int&evals)const
{
int last_evals;
evals= 1;
lev= 1;
do{
lev++;
last_evals= evals;
evals= calcNumEvaluations(lev);
}while(lev<uquad.numLevels()&&evals<=max_evals);
lev--;
evals= last_evals;
}


/*:15*/
;

/*:1*/
