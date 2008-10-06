/*1:*/
#line 6 "./global_check.cweb"

#include "SymSchurDecomp.h"

#include "global_check.h"

#include "smolyak.h"
#include "product.h"
#include "quasi_mcarlo.h"

#include "cpplapack.h"

/*2:*/
#line 31 "./global_check.cweb"

ResidFunction::ResidFunction(const Approximation&app)
:VectorFunction(app.getModel().nexog(),app.getModel().numeq()),approx(app),
model(app.getModel().clone()),
yplus(NULL),ystar(NULL),u(NULL),hss(NULL)
{
}

/*:2*/
#line 17 "./global_check.cweb"
;
/*3:*/
#line 40 "./global_check.cweb"

ResidFunction::ResidFunction(const ResidFunction&rf)
:VectorFunction(rf),approx(rf.approx),model(rf.model->clone()),
yplus(NULL),ystar(NULL),u(NULL),hss(NULL)
{
if(rf.yplus)
yplus= new Vector(*(rf.yplus));
if(rf.ystar)
ystar= new Vector(*(rf.ystar));
if(rf.u)
u= new Vector(*(rf.u));
if(rf.hss)
hss= new FTensorPolynomial(*(rf.hss));
}


/*:3*/
#line 18 "./global_check.cweb"
;
/*4:*/
#line 57 "./global_check.cweb"

ResidFunction::~ResidFunction()
{
delete model;
/*5:*/
#line 65 "./global_check.cweb"

if(yplus)
delete yplus;
if(ystar)
delete ystar;
if(u)
delete u;
if(hss)
delete hss;


/*:5*/
#line 61 "./global_check.cweb"
;
}

/*:4*/
#line 19 "./global_check.cweb"
;
/*6:*/
#line 79 "./global_check.cweb"

void ResidFunction::setYU(const Vector&ys,const Vector&xx)
{
/*5:*/
#line 65 "./global_check.cweb"

if(yplus)
delete yplus;
if(ystar)
delete ystar;
if(u)
delete u;
if(hss)
delete hss;


/*:5*/
#line 82 "./global_check.cweb"
;

ystar= new Vector(ys);
u= new Vector(xx);
yplus= new Vector(model->numeq());
approx.getFoldDecisionRule().evaluate(DecisionRule::horner,
*yplus,*ystar,*u);

/*7:*/
#line 103 "./global_check.cweb"

union{const FoldDecisionRule*c;FoldDecisionRule*n;}dr;
dr.c= &(approx.getFoldDecisionRule());
FTensorPolynomial dr_ss(model->nstat()+model->npred(),model->nboth()+model->nforw(),
*(dr.n));

/*:7*/
#line 90 "./global_check.cweb"
;
/*8:*/
#line 110 "./global_check.cweb"

Vector ytmp_star(ConstVector(*yplus,model->nstat(),model->npred()+model->nboth()));
ConstVector ysteady_star(dr.c->getSteady(),model->nstat(),
model->npred()+model->nboth());
ytmp_star.add(-1.0,ysteady_star);

/*:8*/
#line 91 "./global_check.cweb"
;
/*9:*/
#line 117 "./global_check.cweb"

hss= new FTensorPolynomial(dr_ss,ytmp_star);
ConstVector ysteady_ss(dr.c->getSteady(),model->nstat()+model->npred(),
model->nboth()+model->nforw());
if(hss->check(Symmetry(0))){
hss->get(Symmetry(0))->getData().add(1.0,ysteady_ss);
}else{
FFSTensor*ten= new FFSTensor(hss->nrows(),hss->nvars(),0);
ten->getData()= ysteady_ss;
hss->insert(ten);
}

/*:9*/
#line 92 "./global_check.cweb"
;
}

/*:6*/
#line 20 "./global_check.cweb"
;
/*10:*/
#line 132 "./global_check.cweb"

void ResidFunction::eval(const Vector&point,const ParameterSignal&sig,Vector&out)
{
KORD_RAISE_IF(point.length()!=hss->nvars(),
"Wrong dimension of input vector in ResidFunction::eval");
KORD_RAISE_IF(out.length()!=model->numeq(),
"Wrong dimension of output vector in ResidFunction::eval");
Vector yss(hss->nrows());
hss->evalHorner(yss,point);
model->evaluateSystem(out,*ystar,*yplus,yss,*u);
}

/*:10*/
#line 21 "./global_check.cweb"
;
/*11:*/
#line 148 "./global_check.cweb"

void GlobalChecker::check(const Quadrature&quad,int level,
const ConstVector&ys,const ConstVector&x,Vector&out)
{
for(int ifunc= 0;ifunc<vfs.getNum();ifunc++)
((GResidFunction&)(vfs.getFunc(ifunc))).setYU(ys,x);
quad.integrate(vfs,level,out);
}

/*:11*/
#line 22 "./global_check.cweb"
;
/*12:*/
#line 165 "./global_check.cweb"

void GlobalChecker::check(int max_evals,const ConstTwoDMatrix&y,
const ConstTwoDMatrix&x,TwoDMatrix&out)
{
JournalRecordPair pa(journal);
pa<<"Checking approximation error for "<<y.ncols()
<<" states with at most "<<max_evals<<" evaluations"<<endrec;

/*13:*/
#line 182 "./global_check.cweb"

GaussHermite gh;

SmolyakQuadrature dummy_sq(model.nexog(),1,gh);
int smol_evals;
int smol_level;
dummy_sq.designLevelForEvals(max_evals,smol_level,smol_evals);

ProductQuadrature dummy_pq(model.nexog(),gh);
int prod_evals;
int prod_level;
dummy_pq.designLevelForEvals(max_evals,prod_level,prod_evals);

bool take_smolyak= (smol_evals<prod_evals)&&(smol_level>=prod_level-1);

/*:13*/
#line 173 "./global_check.cweb"
;
Quadrature*quad;
int lev;
/*14:*/
#line 198 "./global_check.cweb"

if(take_smolyak){
quad= new SmolyakQuadrature(model.nexog(),smol_level,gh);
lev= smol_level;
JournalRecord rec(journal);
rec<<"Selected Smolyak (level,evals)=("<<smol_level<<","
<<smol_evals<<") over product ("<<prod_level<<","
<<prod_evals<<")"<<endrec;
}else{
quad= new ProductQuadrature(model.nexog(),gh);
lev= prod_level;
JournalRecord rec(journal);
rec<<"Selected product (level,evals)=("<<prod_level<<","
<<prod_evals<<") over Smolyak ("<<smol_level<<","
<<smol_evals<<")"<<endrec;
}

/*:14*/
#line 176 "./global_check.cweb"
;
/*15:*/
#line 216 "./global_check.cweb"

int first_row= (y.nrows()==model.numeq())?model.nstat():0;
ConstTwoDMatrix ysmat(y,first_row,0,model.npred()+model.nboth(),y.ncols());
for(int j= 0;j<y.ncols();j++){
ConstVector yj(ysmat,j);
ConstVector xj(x,j);
Vector outj(out,j);
check(*quad,lev,yj,xj,outj);
}



/*:15*/
#line 177 "./global_check.cweb"
;
delete quad;
}

/*:12*/
#line 23 "./global_check.cweb"
;
/*16:*/
#line 233 "./global_check.cweb"

void GlobalChecker::checkAlongShocksAndSave(FILE*fd,const char*prefix,
int m,double mult,int max_evals)
{
JournalRecordPair pa(journal);
pa<<"Calculating errors along shocks +/- "
<<mult<<" std errors, granularity "<<m<<endrec;
/*17:*/
#line 250 "./global_check.cweb"

TwoDMatrix y_mat(model.numeq(),2*m*model.nexog()+1);
for(int j= 0;j<2*m*model.nexog()+1;j++){
Vector yj(y_mat,j);
yj= (const Vector&)model.getSteady();
}

/*:17*/
#line 240 "./global_check.cweb"
;
/*18:*/
#line 258 "./global_check.cweb"

TwoDMatrix exo_mat(model.nexog(),2*m*model.nexog()+1);
exo_mat.zeros();
for(int ishock= 0;ishock<model.nexog();ishock++){
double max_sigma= sqrt(model.getVcov().get(ishock,ishock));
for(int j= 0;j<2*m;j++){
int jmult= (j<m)?j-m:j-m+1;
exo_mat.get(ishock,1+2*m*ishock+j)= 
mult*jmult*max_sigma/m;
}
}

/*:18*/
#line 241 "./global_check.cweb"
;

TwoDMatrix errors(model.numeq(),2*m*model.nexog()+1);
check(max_evals,y_mat,exo_mat,errors);

/*19:*/
#line 271 "./global_check.cweb"

TwoDMatrix res(model.nexog(),2*m+1);
JournalRecord rec(journal);
rec<<"Shock    value         error"<<endrec;
ConstVector err0(errors,0);
char shock[9];
char erbuf[17];
for(int ishock= 0;ishock<model.nexog();ishock++){
TwoDMatrix err_out(model.numeq(),2*m+1);
sprintf(shock,"%-8s",model.getExogNames().getName(ishock));
for(int j= 0;j<2*m+1;j++){
int jj;
Vector error(err_out,j);
if(j!=m){
if(j<m)
jj= 1+2*m*ishock+j;
else
jj= 1+2*m*ishock+j-1;
ConstVector coljj(errors,jj);
error= coljj;
}else{
jj= 0;
error= err0;
}
JournalRecord rec1(journal);
sprintf(erbuf,"%12.7g    ",error.getMax());
rec1<<shock<<" "<<exo_mat.get(ishock,jj)
<<"\t"<<erbuf<<endrec;
}
char tmp[100];
sprintf(tmp,"%s_shock_%s_errors",prefix,model.getExogNames().getName(ishock));
err_out.writeMat4(fd,tmp);
}


/*:19*/
#line 246 "./global_check.cweb"
;
}

/*:16*/
#line 24 "./global_check.cweb"
;
/*20:*/
#line 320 "./global_check.cweb"

void GlobalChecker::checkOnEllipseAndSave(FILE*fd,const char*prefix,
int m,double mult,int max_evals)
{
JournalRecordPair pa(journal);
pa<<"Calculating errors at "<<m
<<" ellipse points scaled by "<<mult<<endrec;
/*21:*/
#line 338 "./global_check.cweb"

TwoDMatrix*ycov= approx.calcYCov();
TwoDMatrix ycovpred((const TwoDMatrix&)*ycov,model.nstat(),model.nstat(),
model.npred()+model.nboth(),model.npred()+model.nboth());
delete ycov;
SymSchurDecomp ssd(ycovpred);
ssd.correctDefinitness(1.e-05);
TwoDMatrix ycovfac(ycovpred.nrows(),ycovpred.ncols());
KORD_RAISE_IF(!ssd.isPositiveSemidefinite(),
"Covariance matrix of the states not positive \
				  semidefinite in GlobalChecker::checkOnEllipseAndSave");
ssd.getFactor(ycovfac);


/*:21*/
#line 327 "./global_check.cweb"
;
/*22:*/
#line 363 "./global_check.cweb"

int d= model.npred()+model.nboth()-1;
TwoDMatrix ymat(model.npred()+model.nboth(),(d==0)?2:m);
if(d==0){
ymat.get(0,0)= 1;
ymat.get(0,1)= -1;
}else{
int icol= 0;
ReversePerScheme ps;
QMCarloCubeQuadrature qmc(d,m,ps);
qmcpit beg= qmc.start(m);
qmcpit end= qmc.end(m);
for(qmcpit run= beg;run!=end;++run,icol++){
Vector ycol(ymat,icol);
Vector x(run.point());
x.mult(2*M_PI);
ycol[0]= 1;
for(int i= 0;i<d;i++){
Vector subsphere(ycol,0,i+1);
subsphere.mult(cos(x[i]));
ycol[i+1]= sin(x[i]);
}
}
}

/*:22*/
#line 328 "./global_check.cweb"
;
/*23:*/
#line 392 "./global_check.cweb"

TwoDMatrix umat(model.nexog(),ymat.ncols());
umat.zeros();
ymat.mult(mult);
ymat.multLeft(ycovfac);
ConstVector ys(model.getSteady(),model.nstat(),
model.npred()+model.nboth());
for(int icol= 0;icol<ymat.ncols();icol++){
Vector ycol(ymat,icol);
ycol.add(1.0,ys);
}

/*:23*/
#line 329 "./global_check.cweb"
;
/*24:*/
#line 405 "./global_check.cweb"

TwoDMatrix out(model.numeq(),ymat.ncols());
check(max_evals,ymat,umat,out);

char tmp[100];
sprintf(tmp,"%s_ellipse_points",prefix);
ymat.writeMat4(fd,tmp);
sprintf(tmp,"%s_ellipse_errors",prefix);
out.writeMat4(fd,tmp);

/*:24*/
#line 330 "./global_check.cweb"
;
}


/*:20*/
#line 25 "./global_check.cweb"
;
/*25:*/
#line 418 "./global_check.cweb"

void GlobalChecker::checkAlongSimulationAndSave(FILE*fd,const char*prefix,
int m,int max_evals)
{
JournalRecordPair pa(journal);
pa<<"Calculating errors at "<<m
<<" simulated points"<<endrec;
RandomShockRealization sr(model.getVcov(),system_random_generator.int_uniform());
TwoDMatrix*y= approx.getFoldDecisionRule().simulate(DecisionRule::horner,
m,model.getSteady(),sr);
TwoDMatrix x(model.nexog(),m);
x.zeros();
TwoDMatrix out(model.numeq(),m);
check(max_evals,*y,x,out);

char tmp[100];
sprintf(tmp,"%s_simul_points",prefix);
y->writeMat4(fd,tmp);
sprintf(tmp,"%s_simul_errors",prefix);
out.writeMat4(fd,tmp);

delete y;
}


/*:25*/
#line 26 "./global_check.cweb"
;

/*:1*/
