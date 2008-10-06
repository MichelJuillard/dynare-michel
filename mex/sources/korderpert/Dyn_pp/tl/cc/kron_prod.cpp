/*1:*/
#line 5 "./kron_prod.cweb"

#include "kron_prod.h"
#include "tl_exception.h"

#include <stdio.h> 

/*2:*/
#line 39 "./kron_prod.cweb"

KronProdDimens::KronProdDimens(const KronProdDimens&kd,int i)
:rows((i==0||i==kd.dimen()-1)?(2):(3)),
cols((i==0||i==kd.dimen()-1)?(2):(3))
{
TL_RAISE_IF(i<0||i>=kd.dimen(),
"Wrong index for pickup in KronProdDimens constructor");

int kdim= kd.dimen();
if(i==0){
/*3:*/
#line 61 "./kron_prod.cweb"

rows[0]= kd.rows[0];
rows[1]= kd.rows.mult(1,kdim);
cols[0]= kd.cols[0];
cols[1]= rows[1];

/*:3*/
#line 49 "./kron_prod.cweb"
;
}else if(i==kdim-1){
/*4:*/
#line 71 "./kron_prod.cweb"

rows[0]= kd.cols.mult(0,kdim-1);
rows[1]= kd.rows[kdim-1];
cols[0]= rows[0];
cols[1]= kd.cols[kdim-1];

/*:4*/
#line 51 "./kron_prod.cweb"
;
}else{
/*5:*/
#line 83 "./kron_prod.cweb"

rows[0]= kd.cols.mult(0,i);
cols[0]= rows[0];
rows[1]= kd.rows[i];
cols[1]= kd.cols[i];
cols[2]= kd.rows.mult(i+1,kdim);
rows[2]= cols[2];


/*:5*/
#line 53 "./kron_prod.cweb"
;
}
}

/*:2*/
#line 11 "./kron_prod.cweb"
;
/*6:*/
#line 95 "./kron_prod.cweb"

void KronProd::checkDimForMult(const ConstTwoDMatrix&in,const TwoDMatrix&out)const
{
int my_rows;
int my_cols;
kpd.getRC(my_rows,my_cols);
TL_RAISE_IF(in.nrows()!=out.nrows()||in.ncols()!=my_rows,
"Wrong dimensions for KronProd in KronProd::checkDimForMult");
}

/*:6*/
#line 12 "./kron_prod.cweb"
;
/*7:*/
#line 108 "./kron_prod.cweb"

void KronProd::kronMult(const ConstVector&v1,const ConstVector&v2,
Vector&res)
{
TL_RAISE_IF(res.length()!=v1.length()*v2.length(),
"Wrong vector lengths in KronProd::kronMult");
res.zeros();
for(int i= 0;i<v1.length();i++){
Vector sub(res,i*v2.length(),v2.length());
sub.add(v1[i],v2);
}
}


/*:7*/
#line 13 "./kron_prod.cweb"
;
/*8:*/
#line 123 "./kron_prod.cweb"

void KronProdAll::setMat(int i,const TwoDMatrix&m)
{
matlist[i]= &m;
kpd.setRC(i,m.nrows(),m.ncols());
}

/*:8*/
#line 14 "./kron_prod.cweb"
;
/*9:*/
#line 131 "./kron_prod.cweb"

void KronProdAll::setUnit(int i,int n)
{
matlist[i]= NULL;
kpd.setRC(i,n,n);
}

/*:9*/
#line 15 "./kron_prod.cweb"
;
/*10:*/
#line 139 "./kron_prod.cweb"

bool KronProdAll::isUnit()const
{
int i= 0;
while(i<dimen()&&matlist[i]==NULL)
i++;
return i==dimen();
}

/*:10*/
#line 16 "./kron_prod.cweb"
;
/*22:*/
#line 359 "./kron_prod.cweb"

Vector*KronProdAll::multRows(const IntSequence&irows)const
{
TL_RAISE_IF(irows.size()!=dimen(),
"Wrong length of row indices in KronProdAll::multRows");

Vector*last= NULL;
ConstVector*row;
vector<Vector*> to_delete;
for(int i= 0;i<dimen();i++){
int j= dimen()-1-i;
/*23:*/
#line 387 "./kron_prod.cweb"

if(matlist[j])
row= new ConstVector(irows[j],*(matlist[j]));
else{
Vector*aux= new Vector(ncols(j));
aux->zeros();
(*aux)[irows[j]]= 1.0;
to_delete.push_back(aux);
row= new ConstVector(*aux);
}

/*:23*/
#line 370 "./kron_prod.cweb"
;
/*24:*/
#line 402 "./kron_prod.cweb"

if(last){
Vector*newlast;
newlast= new Vector(last->length()*row->length());
kronMult(*row,ConstVector(*last),*newlast);
delete last;
last= newlast;
}else{
last= new Vector(*row);
}


/*:24*/
#line 371 "./kron_prod.cweb"
;
delete row;
}

for(unsigned int i= 0;i<to_delete.size();i++)
delete to_delete[i];

return last;
}

/*:22*/
#line 17 "./kron_prod.cweb"
;
/*11:*/
#line 156 "./kron_prod.cweb"

void KronProdIA::mult(const ConstTwoDMatrix&in,TwoDMatrix&out)const
{
checkDimForMult(in,out);

int id_cols= kpd.cols[0];
ConstTwoDMatrix a(mat);

for(int i= 0;i<id_cols;i++){
TwoDMatrix outi(out,i*a.ncols(),a.ncols());
ConstTwoDMatrix ini(in,i*a.nrows(),a.nrows());
outi.mult(ini,a);
}
}

/*:11*/
#line 18 "./kron_prod.cweb"
;
/*12:*/
#line 172 "./kron_prod.cweb"

KronProdAI::KronProdAI(const KronProdIAI&kpiai)
:KronProd(KronProdDimens(2)),mat(kpiai.mat)
{
kpd.rows[0]= mat.nrows();
kpd.cols[0]= mat.ncols();
kpd.rows[1]= kpiai.kpd.rows[2];
kpd.cols[1]= kpiai.kpd.cols[2];
}


/*:12*/
#line 19 "./kron_prod.cweb"
;
/*13:*/
#line 202 "./kron_prod.cweb"

void KronProdAI::mult(const ConstTwoDMatrix&in,TwoDMatrix&out)const
{
checkDimForMult(in,out);

int id_cols= kpd.cols[1];
ConstTwoDMatrix a(mat);

if(in.getLD()==in.nrows()){
ConstTwoDMatrix in_resh(in.nrows()*id_cols,a.nrows(),in.getData().base());
TwoDMatrix out_resh(in.nrows()*id_cols,a.ncols(),out.getData().base());
out_resh.mult(in_resh,a);
}else{
out.zeros();
for(int i= 0;i<a.ncols();i++){
TwoDMatrix outi(out,i*id_cols,id_cols);
for(int j= 0;j<a.nrows();j++){
ConstTwoDMatrix ini(in,j*id_cols,id_cols);
outi.add(a.get(j,i),ini);
}
}
}
}


/*:13*/
#line 20 "./kron_prod.cweb"
;
/*14:*/
#line 240 "./kron_prod.cweb"

void KronProdIAI::mult(const ConstTwoDMatrix&in,TwoDMatrix&out)const
{
checkDimForMult(in,out);

int id_cols= kpd.cols[0];

KronProdAI akronid(*this);
int in_bl_width;
int out_bl_width;
akronid.kpd.getRC(in_bl_width,out_bl_width);

for(int i= 0;i<id_cols;i++){
TwoDMatrix outi(out,i*out_bl_width,out_bl_width);
ConstTwoDMatrix ini(in,i*in_bl_width,in_bl_width);
akronid.mult(ini,outi);
}
}

/*:14*/
#line 21 "./kron_prod.cweb"
;
/*15:*/
#line 273 "./kron_prod.cweb"

void KronProdAll::mult(const ConstTwoDMatrix&in,TwoDMatrix&out)const
{
/*16:*/
#line 287 "./kron_prod.cweb"

if(isUnit()){
out.zeros();
out.add(1.0,in);
return;
}

/*:16*/
#line 276 "./kron_prod.cweb"
;
/*17:*/
#line 297 "./kron_prod.cweb"

bool is_zero= false;
for(int i= 0;i<dimen()&&!is_zero;i++)
is_zero= matlist[i]&&matlist[i]->isZero();
if(is_zero||in.isZero()){
out.zeros();
return;
}

/*:17*/
#line 277 "./kron_prod.cweb"
;
/*18:*/
#line 307 "./kron_prod.cweb"

if(dimen()==1){
if(matlist[0])
out.mult(in,ConstTwoDMatrix(*(matlist[0])));
return;
}

/*:18*/
#line 278 "./kron_prod.cweb"
;
int c;
TwoDMatrix*last= NULL;
/*19:*/
#line 317 "./kron_prod.cweb"

if(matlist[0]){
KronProdAI akronid(*this);
c= akronid.kpd.ncols();
last= new TwoDMatrix(in.nrows(),c);
akronid.mult(in,*last);
}else{
last= new TwoDMatrix(in.nrows(),in.ncols(),in.getData().base());
}

/*:19*/
#line 281 "./kron_prod.cweb"
;
/*20:*/
#line 331 "./kron_prod.cweb"

for(int i= 1;i<dimen()-1;i++){
if(matlist[i]){
KronProdIAI interkron(*this,i);
c= interkron.kpd.ncols();
TwoDMatrix*newlast= new TwoDMatrix(in.nrows(),c);
interkron.mult(*last,*newlast);
delete last;
last= newlast;
}
}

/*:20*/
#line 282 "./kron_prod.cweb"
;
/*21:*/
#line 346 "./kron_prod.cweb"

if(matlist[dimen()-1]){
KronProdIA idkrona(*this);
idkrona.mult(*last,out);
}else{
out= *last;
}
delete last;

/*:21*/
#line 283 "./kron_prod.cweb"
;
}

/*:15*/
#line 22 "./kron_prod.cweb"
;
/*25:*/
#line 420 "./kron_prod.cweb"

void KronProdAllOptim::optimizeOrder()
{
for(int i= 0;i<dimen();i++){
int swaps= 0;
for(int j= 0;j<dimen()-1;j++){
if(((double)kpd.rows[j])/kpd.cols[j]<((double)kpd.rows[j+1])/kpd.cols[j+1]){
/*26:*/
#line 438 "./kron_prod.cweb"

int s= kpd.rows[j+1];
kpd.rows[j+1]= kpd.rows[j];
kpd.rows[j]= s;
s= kpd.cols[j+1];
kpd.cols[j+1]= kpd.cols[j];
kpd.cols[j]= s;
const TwoDMatrix*m= matlist[j+1];
matlist[j+1]= matlist[j];
matlist[j]= m;

/*:26*/
#line 427 "./kron_prod.cweb"
;
/*27:*/
#line 450 "./kron_prod.cweb"

s= oper.getMap()[j+1];
oper.getMap()[j+1]= oper.getMap()[j];
oper.getMap()[j]= s;
swaps++;


/*:27*/
#line 428 "./kron_prod.cweb"
;
}
}
if(swaps==0){
return;
}
}
}

/*:25*/
#line 23 "./kron_prod.cweb"
;

/*:1*/
