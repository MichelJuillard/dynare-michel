/*1:*/

#include "kord_exception.h"
#include "approximation.h"
#include "first_order.h"
#include "korder_stoch.h"

/*2:*/

ZAuxContainer::ZAuxContainer(const _Ctype*gss,int ngss,int ng,int ny,int nu)
:StackContainer<FGSTensor> (4,1)
{
	stack_sizes[0]= ngss;stack_sizes[1]= ng;
	stack_sizes[2]= ny;stack_sizes[3]= nu;
	conts[0]= gss;
	calculateOffsets();
}


/*:2*/
;
/*3:*/

ZAuxContainer::itype ZAuxContainer::getType(int i,const Symmetry&s)const
{
	if(i==0)
		if(s[2]> 0)
			return zero;
		else
			return matrix;
		return zero;
}


/*:3*/
;
/*4:*/

Approximation::Approximation(DynamicModel&m,Journal&j,int ns)
:model(m),journal(j),rule_ders(NULL),rule_ders_ss(NULL),fdr(NULL),udr(NULL),
ypart(model.nstat(),model.npred(),model.nboth(),model.nforw()),
mom(UNormalMoments(model.order(),model.getVcov())),nvs(4),steps(ns),
ss(ypart.ny(),steps+1)
{
	nvs[0]= ypart.nys();nvs[1]= model.nexog();
	nvs[2]= model.nexog();nvs[3]= 1;
	
	ss.nans();
}

/*:4*/
;
/*5:*/

Approximation::~Approximation()
{
	if(rule_ders_ss)delete rule_ders_ss;
	if(rule_ders)delete rule_ders;
	if(fdr)delete fdr;
	if(udr)delete udr;
}

/*:5*/
;
/*6:*/

const FoldDecisionRule&Approximation::getFoldDecisionRule()const
{
	KORD_RAISE_IF(fdr==NULL,
		"Folded decision rule has not been created in Approximation::getFoldDecisionRule");
	return*fdr;
}


/*:6*/
;
/*7:*/

const UnfoldDecisionRule&Approximation::getUnfoldDecisionRule()const
{
	KORD_RAISE_IF(udr==NULL,
		"Unfolded decision rule has not been created in Approximation::getUnfoldDecisionRule");
	return*udr;
}


/*:7*/
;
/*8:*/

void Approximation::approxAtSteady()
{
	model.calcDerivativesAtSteady();
	FirstOrder fo(model.nstat(),model.npred(),model.nboth(),model.nforw(),
		model.nexog(),*(model.getModelDerivatives().get(Symmetry(1))),
		journal);
	KORD_RAISE_IF_X(!fo.isStable(),
		"The model is not Blanchard-Kahn stable",
		KORD_MD_NOT_STABLE);
	
	if(model.order()>=2){
		KOrder korder(model.nstat(),model.npred(),model.nboth(),model.nforw(),
			model.getModelDerivatives(),fo.getGy(),fo.getGu(),
			model.getVcov(),journal);
		korder.switchToFolded();
		for(int k= 2;k<=model.order();k++)
			korder.performStep<KOrder::fold> (k);
		
		saveRuleDerivs(korder.getFoldDers());
	}else{
		FirstOrderDerivs<KOrder::fold> fo_ders(fo);
		saveRuleDerivs(fo_ders);
	}
	//check(0.0);
}

/*:8*/
;
/*9:*/

void Approximation::walkStochSteady()
{
	/*10:*/
	
	model.solveDeterministicSteady();
	approxAtSteady();
	Vector steady0(ss,0);
	steady0= (const Vector&)model.getSteady();
	
	/*:10*/
	;
	double sigma_so_far= 0.0;
	double dsigma= (steps==0)?0.0:1.0/steps;
	for(int i= 1;i<=steps;i++){
		JournalRecordPair pa(journal);
		pa<<"Approximation about stochastic steady for sigma="<<sigma_so_far+dsigma<<endrec;
		
		Vector last_steady((const Vector&)model.getSteady());
		
		/*11:*/
		
		DRFixPoint<KOrder::fold> fp(*rule_ders,ypart,model.getSteady(),dsigma);
		bool converged= fp.calcFixPoint(DecisionRule::horner,model.getSteady());
		JournalRecord rec(journal);
		rec<<"Fix point calcs: iter="<<fp.getNumIter()<<", newton_iter="
			<<fp.getNewtonTotalIter()<<", last_newton_iter="<<fp.getNewtonLastIter()<<".";
		if(converged)
			rec<<" Converged."<<endrec;
		else{
			rec<<" Not converged!!"<<endrec;
			KORD_RAISE_X("Fix point calculation not converged",KORD_FP_NOT_CONV);
		}
		Vector steadyi(ss,i);
		steadyi= (const Vector&)model.getSteady();
		
		/*:11*/
		;
		/*12:*/
		
		Vector dy((const Vector&)model.getSteady());
		dy.add(-1.0,last_steady);
		
		StochForwardDerivs<KOrder::fold> hh(ypart,model.nexog(),*rule_ders_ss,mom,dy,
			dsigma,sigma_so_far);
		JournalRecord rec1(journal);
		rec1<<"Calculation of g** expectations done"<<endrec;
		
		
		/*:12*/
		;
		/*13:*/
		
		model.calcDerivativesAtSteady();
		KOrderStoch korder_stoch(ypart,model.nexog(),model.getModelDerivatives(),
			hh,journal);
		for(int d= 1;d<=model.order();d++){
			korder_stoch.performStep<KOrder::fold> (d);
		}
		saveRuleDerivs(korder_stoch.getFoldDers());
		
		
		/*:13*/
		;
		
		//check(sigma_so_far+dsigma);
		sigma_so_far+= dsigma;
	}
	
	/*14:*/
	
	if(fdr){
		delete fdr;
		fdr= NULL;
	}
	if(udr){
		delete udr;
		udr= NULL;
	}
	
	fdr= new FoldDecisionRule(*rule_ders,ypart,model.nexog(),
		model.getSteady(),1.0-sigma_so_far);
	if(steps==0){
		/*15:*/
		
		DRFixPoint<KOrder::fold> fp(*rule_ders,ypart,model.getSteady(),1.0);
		bool converged= fp.calcFixPoint(DecisionRule::horner,model.getSteady());
		JournalRecord rec(journal);
		rec<<"Fix point calcs: iter="<<fp.getNumIter()<<", newton_iter="
			<<fp.getNewtonTotalIter()<<", last_newton_iter="<<fp.getNewtonLastIter()<<".";
		if(converged)
			rec<<" Converged."<<endrec;
		else{
			rec<<" Not converged!!"<<endrec;
			KORD_RAISE_X("Fix point calculation not converged",KORD_FP_NOT_CONV);
		}
		
		{
			JournalRecordPair recp(journal);
			recp<<"Centralizing about fix-point."<<endrec;
			FoldDecisionRule*dr_backup= fdr;
			fdr= new FoldDecisionRule(*dr_backup,model.getSteady());
			delete dr_backup;
		}
		
		
		/*:15*/
		;
	}
	
	
	/*:14*/
	;
}

/*:9*/
;
/*16:*/

void Approximation::saveRuleDerivs(const FGSContainer&g)
{
	if(rule_ders){
		delete rule_ders;
		delete rule_ders_ss;
	}
	rule_ders= new FGSContainer(g);
	rule_ders_ss= new FGSContainer(4);
	for(FGSContainer::iterator run= (*rule_ders).begin();run!=(*rule_ders).end();++run){
		FGSTensor*ten= new FGSTensor(ypart.nstat+ypart.npred,ypart.nyss(),*((*run).second));
		rule_ders_ss->insert(ten);
	}
}

/*:16*/
;
/*17:*/

void Approximation::calcStochShift(Vector&out,double at_sigma)const
{
	KORD_RAISE_IF(out.length()!=ypart.ny(),
		"Wrong length of output vector for Approximation::calcStochShift");
	out.zeros();
	
	ZAuxContainer zaux(rule_ders_ss,ypart.nyss(),ypart.ny(),
		ypart.nys(),model.nexog());
	
	int dfac= 1;
	for(int d= 1;d<=rule_ders->getMaxDim();d++,dfac*= d){
		if(KOrder::is_even(d)){
			Symmetry sym(0,d,0,0);
			/*18:*/
			
			FGSTensor*ten= new FGSTensor(ypart.ny(),TensorDimens(sym,nvs));
			ten->zeros();
			for(int l= 1;l<=d;l++){
				const FSSparseTensor*f= model.getModelDerivatives().get(Symmetry(l));
				zaux.multAndAdd(*f,*ten);
			}
			
			/*:18*/
			;
			/*19:*/
			
			FGSTensor*tmp= new FGSTensor(ypart.ny(),TensorDimens(Symmetry(0,0,0,0),nvs));
			tmp->zeros();
			ten->contractAndAdd(1,*tmp,*(mom.get(Symmetry(d))));
			
			out.add(pow(at_sigma,d)/dfac,tmp->getData());
			delete ten;
			delete tmp;
			
			
			/*:19*/
			;
		}
	}
}

/*:17*/
;
/*20:*/

void Approximation::check(double at_sigma)const
{
	Vector stoch_shift(ypart.ny());
	Vector system_resid(ypart.ny());
	Vector xx(model.nexog());
	xx.zeros();
	model.evaluateSystem(system_resid,model.getSteady(),xx);
	calcStochShift(stoch_shift,at_sigma);
	stoch_shift.add(1.0,system_resid);
	JournalRecord rec1(journal);
	rec1<<"Error of current approximation for shocks at sigma "<<at_sigma
		<<" is "<<stoch_shift.getMax()<<endrec;
	calcStochShift(stoch_shift,1.0);
	stoch_shift.add(1.0,system_resid);
	JournalRecord rec2(journal);
	rec2<<"Error of current approximation for full shocks is "<<stoch_shift.getMax()<<endrec;
}

/*:20*/
;
/*21:*/

TwoDMatrix*Approximation::calcYCov()const
{
	const TwoDMatrix&gy= *(rule_ders->get(Symmetry(1,0,0,0)));
	const TwoDMatrix&gu= *(rule_ders->get(Symmetry(0,1,0,0)));
	TwoDMatrix G(model.numeq(),model.numeq());
	G.zeros();
	G.place(gy,0,model.nstat());
	TwoDMatrix B((const TwoDMatrix&)G);
	B.mult(-1.0);
	TwoDMatrix C(G,"transpose");
	TwoDMatrix A(model.numeq(),model.numeq());
	A.zeros();
	for(int i= 0;i<model.numeq();i++)
		A.get(i,i)= 1.0;
	
	TwoDMatrix guSigma(gu,model.getVcov());
	TwoDMatrix guTrans(gu,"transpose");
	TwoDMatrix*X= new TwoDMatrix(guSigma,guTrans);
	
	GeneralSylvester gs(1,model.numeq(),model.numeq(),0,
		A.base(),B.base(),C.base(),X->base());
	gs.solve();
	
	return X;
}

/*:21*/
;

/*:1*/
