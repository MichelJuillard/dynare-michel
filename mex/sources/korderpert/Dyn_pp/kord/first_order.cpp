/*1:*/


#include "kord_exception.h"
#include "first_order.h"
#include "cpplapack.h"

/*2:*/

int order_eigs(const double*alphar,const double*alphai,const double*beta)
{
	return(*alphar**alphar+*alphai**alphai<*beta**beta);
}


/*:2*/
;
/*3:*/

void FirstOrder::solve(const TwoDMatrix&fd)
{
	JournalRecordPair pa(journal);
	pa<<"Recovering first order derivatives "<<endrec;
	
	/*4:*/
	
	/*5:*/
	
	int off= 0;
	ConstTwoDMatrix fyplus(fd,off,ypart.nyss());
	off+= ypart.nyss();
	ConstTwoDMatrix fyszero(fd,off,ypart.nstat);
	off+= ypart.nstat;
	ConstTwoDMatrix fypzero(fd,off,ypart.npred);
	off+= ypart.npred;
	ConstTwoDMatrix fybzero(fd,off,ypart.nboth);
	off+= ypart.nboth;
	ConstTwoDMatrix fyfzero(fd,off,ypart.nforw);
	off+= ypart.nforw;
	ConstTwoDMatrix fymins(fd,off,ypart.nys());
	off+= ypart.nys();
	ConstTwoDMatrix fuzero(fd,off,nu);
	off+= nu;
	
	/*:5*/
	;
	/*6:*/
	
	int n= ypart.ny()+ypart.nboth;
	TwoDMatrix matD(n,n);
	matD.zeros();
	matD.place(fypzero,0,0);
	matD.place(fybzero,0,ypart.npred);
	matD.place(fyplus,0,ypart.nys()+ypart.nstat);
	for(int i= 0;i<ypart.nboth;i++)
		matD.get(ypart.ny()+i,ypart.npred+i)= 1.0;
	
	
	/*:6*/
	;
	/*7:*/
	
	TwoDMatrix matE(n,n);
	matE.zeros();
	matE.place(fymins,0,0);
	matE.place(fyszero,0,ypart.nys());
	matE.place(fyfzero,0,ypart.nys()+ypart.nstat+ypart.nboth);
	for(int i= 0;i<ypart.nboth;i++)
		matE.get(ypart.ny()+i,ypart.nys()+ypart.nstat+i)= -1.0;
	matE.mult(-1.0);
	
	/*:7*/
	;
	/*8:*/
	
	TwoDMatrix vsl(n,n);
	TwoDMatrix vsr(n,n);
	int lwork= 100*n+16;
	Vector work(lwork);
	IntSequence bwork(n);
	int info;
	LAPACK_dgges("N","V","S",order_eigs,&n,matE.getData().base(),&n,
		matD.getData().base(),&n,&sdim,alphar.base(),alphai.base(),
		beta.base(),vsl.getData().base(),&n,vsr.getData().base(),&n,
		work.base(),&lwork,&(bwork[0]),&info);
	bk_cond= (sdim==ypart.nys());
	
	
	/*:8*/
	;
	/*9:*/
	
	ConstGeneralMatrix z11(vsr,0,0,ypart.nys(),ypart.nys());
	ConstGeneralMatrix z12(vsr,0,ypart.nys(),ypart.nys(),n-ypart.nys());
	ConstGeneralMatrix z21(vsr,ypart.nys(),0,n-ypart.nys(),ypart.nys());
	ConstGeneralMatrix z22(vsr,ypart.nys(),ypart.nys(),n-ypart.nys(),n-ypart.nys());
	
	/*:9*/
	;
	/*10:*/
	
	GeneralMatrix sfder(z12,"transpose");
	z22.multInvLeftTrans(sfder);
	sfder.mult(-1);
	
	/*:10*/
	;
	/*11:*/
	
	ConstGeneralMatrix s11(matE,0,0,ypart.nys(),ypart.nys());
	ConstGeneralMatrix t11(matD,0,0,ypart.nys(),ypart.nys());
	GeneralMatrix dumm(s11,"transpose");
	z11.multInvLeftTrans(dumm);
	GeneralMatrix preder(dumm,"transpose");
	t11.multInvLeft(preder);
	preder.multLeft(z11);
	
	/*:11*/
	;
	/*12:*/
	
	gy.place(preder,ypart.nstat,0);
	GeneralMatrix sder(sfder,0,0,ypart.nstat,ypart.nys());
	gy.place(sder,0,0);
	GeneralMatrix fder(sfder,ypart.nstat+ypart.nboth,0,ypart.nforw,ypart.nys());
	gy.place(fder,ypart.nstat+ypart.nys(),0);
	
	/*:12*/
	;
	/*13:*/
	
	GeneralMatrix bder((const GeneralMatrix&)sfder,ypart.nstat,0,ypart.nboth,ypart.nys());
	GeneralMatrix bder2(preder,ypart.npred,0,ypart.nboth,ypart.nys());
	bder.add(-1,bder2);
	b_error= bder.getData().getMax();
	
	/*:13*/
	;
	
	
	/*:4*/
	;
	/*14:*/
	
	GeneralMatrix matA(ypart.ny(),ypart.ny());
	matA.zeros();
	ConstGeneralMatrix gss(gy,ypart.nstat+ypart.npred,0,ypart.nyss(),ypart.nys());
	GeneralMatrix aux(fyplus,gss);
	matA.place(aux,0,ypart.nstat);
	ConstGeneralMatrix fyzero(fd,0,ypart.nyss(),ypart.ny(),ypart.ny());
	matA.add(1.0,fyzero);
	gu.zeros();
	gu.add(-1.0,fuzero);
	ConstGeneralMatrix(matA).multInvLeft(gu);
	
	/*:14*/
	;
	journalEigs();
	
	if(!gy.isFinite()||!gu.isFinite()){
		throw KordException(__FILE__,__LINE__,
			"NaN or Inf asserted in first order derivatives in FirstOrder::solve");
	}
}

/*:3*/
;
/*15:*/

void FirstOrder::journalEigs()
{
	if(bk_cond){
		JournalRecord jr(journal);
		jr<<"Blanchard-Kahn conditition satisfied, model stable"<<endrec;
	}else{
		JournalRecord jr(journal);
		jr<<"Blanchard-Kahn condition not satisfied, model not stable: sdim="<<sdim
			<<" "<<"npred="<<ypart.nys()<<endrec;
	}
	if(!bk_cond){
		for(int i= 0;i<alphar.length();i++){
			if(i==sdim||i==ypart.nys()){
				JournalRecord jr(journal);
				jr<<"---------------------------------------------------- ";
				if(i==sdim)
					jr<<"sdim";
				else
					jr<<"npred";
				jr<<endrec;
			}
			JournalRecord jr(journal);
			double mod= sqrt(alphar[i]*alphar[i]+alphai[i]*alphai[i]);
			mod= mod/round(100000*std::abs(beta[i]))*100000;
			jr<<i<<"\t("<<alphar[i]<<","<<alphai[i]<<") / "<<beta[i]
				<<"  \t"<<mod<<endrec;
		}
	}
}


/*:15*/
;

/*:1*/
