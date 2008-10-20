/*1:*/

#ifndef DECISION_RULE_H
#define DECISION_RULE_H

#include "kord_exception.h"
#include "korder.h"
#include "normal_conjugate.h"
#include "mersenne_twister.h"

/*2:*/

class ShockRealization{
public:
	virtual~ShockRealization(){}
	virtual void get(int n,Vector&out)= 0;
	virtual int numShocks()const= 0;
};

/*:2*/
;
/*3:*/

class DecisionRule{
public:
	enum emethod{horner,trad};
	virtual~DecisionRule(){}
	virtual TwoDMatrix*simulate(emethod em,int np,const Vector&ystart,
		ShockRealization&sr)const= 0;
	virtual void eval(emethod em,Vector&out,const ConstVector&v)const= 0;
	virtual void evaluate(emethod em,Vector&out,const ConstVector&ys,
		const ConstVector&u)const= 0;
	virtual void writeMat4(FILE*fd,const char*prefix)const= 0;
	virtual DecisionRule*centralizedClone(const Vector&fixpoint)const= 0;
	virtual const Vector&getSteady()const= 0;
	virtual int nexog()const= 0;
	virtual const PartitionY&getYPart()const= 0;
};

/*:3*/
;
/*4:*/

template<int t> 
class DecisionRuleImpl:public ctraits<t> ::Tpol,public DecisionRule{
protected:
	typedef typename ctraits<t> ::Tpol _Tparent;
	const Vector ysteady;
	const PartitionY ypart;
	const int nu;
public:
	DecisionRuleImpl(const _Tparent&pol,const PartitionY&yp,int nuu,
		const Vector&ys)
		:ctraits<t> ::Tpol(pol),ysteady(ys),ypart(yp),nu(nuu){}
	DecisionRuleImpl(_Tparent&pol,const PartitionY&yp,int nuu,
		const Vector&ys)
		:ctraits<t> ::Tpol(0,yp.ny(),pol),ysteady(ys),ypart(yp),
		nu(nuu){}
	DecisionRuleImpl(const _Tg&g,const PartitionY&yp,int nuu,
		const Vector&ys,double sigma)
		:ctraits<t> ::Tpol(yp.ny(),yp.nys()+nuu),ysteady(ys),ypart(yp),nu(nuu)
	{fillTensors(g,sigma);}
	DecisionRuleImpl(const DecisionRuleImpl<t> &dr,const ConstVector&fixpoint)
		:ctraits<t> ::Tpol(dr.ypart.ny(),dr.ypart.nys()+dr.nu),
		ysteady(fixpoint),ypart(dr.ypart),nu(dr.nu)
	{centralize(dr);}
	const Vector&getSteady()const
	{return ysteady;}
	/*8:*/
	
	TwoDMatrix*simulate(emethod em,int np,const Vector&ystart,
		ShockRealization&sr)const
	{
		KORD_RAISE_IF(ysteady.length()!=ystart.length(),
			"Start and steady lengths differ in DecisionRuleImpl::simulate");
		TwoDMatrix*res= new TwoDMatrix(ypart.ny(),np);
		
		/*9:*/
		
		Vector dyu(ypart.nys()+nu);
		ConstVector ystart_pred(ystart,ypart.nstat,ypart.nys());
		ConstVector ysteady_pred(ysteady,ypart.nstat,ypart.nys());
		Vector dy(dyu,0,ypart.nys());
		Vector u(dyu,ypart.nys(),nu);
		
		
		/*:9*/
		;
		/*10:*/
		
		dy= ystart_pred;
		dy.add(-1.0,ysteady_pred);
		sr.get(0,u);
		Vector out(*res,0);
		eval(em,out,dyu);
		
		/*:10*/
		;
		/*11:*/
		
		for(int i= 1;i<np;i++){
			ConstVector ym(*res,i-1);
			ConstVector dym(ym,ypart.nstat,ypart.nys());
			dy= dym;
			sr.get(i,u);
			Vector out(*res,i);
			eval(em,out,dyu);
			if(!out.isFinite()){
				if(i+1<np){
					TwoDMatrix rest(*res,i+1,np-i-1);
					rest.zeros();
				}
				return res;
			}
		}
		
		/*:11*/
		;
		/*12:*/
		
		for(int i= 0;i<res->ncols();i++){
			Vector col(*res,i);
			col.add(1.0,ysteady);
		}
		
		
		/*:12*/
		;
		return res;
	}
	
	/*:8*/
	;
	/*13:*/
	
	void evaluate(emethod em,Vector&out,const ConstVector&ys,
		const ConstVector&u)const
	{
		KORD_RAISE_IF(ys.length()!=ypart.nys()||u.length()!=nu,
			"Wrong dimensions of input vectors in DecisionRuleImpl::evaluate");
		KORD_RAISE_IF(out.length()!=ypart.ny(),
			"Wrong dimension of output vector in DecisionRuleImpl::evaluate");
		ConstVector ysteady_pred(ysteady,ypart.nstat,ypart.nys());
		Vector ys_u(ypart.nys()+nu);
		Vector ys_u1(ys_u,0,ypart.nys());
		ys_u1= ys;
		ys_u1.add(-1.0,ysteady_pred);
		Vector ys_u2(ys_u,ypart.nys(),nu);
		ys_u2= u;
		eval(em,out,ys_u);
		out.add(1.0,ysteady);
	}
	
	/*:13*/
	;
	/*14:*/
	
	DecisionRule*centralizedClone(const Vector&fixpoint)const
	{
		return new DecisionRuleImpl<t> (*this,fixpoint);
	}
	
	/*:14*/
	;
	/*16:*/
	
	void writeMat4(FILE*fd,const char*prefix)const
	{
		ctraits<t> ::Tpol::writeMat4(fd,prefix);
		TwoDMatrix dum(ysteady.length(),1);
		dum.getData()= ysteady;
		char tmp[100];
		sprintf(tmp,"%s_ss",prefix);
		ConstTwoDMatrix(dum).writeMat4(fd,tmp);
	}
	
	/*:16*/
	;
	int nexog()const
	{return nu;}
	const PartitionY&getYPart()const
	{return ypart;}
protected:
	/*5:*/
	
	void fillTensors(const _Tg&g,double sigma)
	{
		IntSequence tns(2);
		tns[0]= ypart.nys();tns[1]= nu;
		int dfact= 1;
		for(int d= 0;d<=g.getMaxDim();d++,dfact*= d){
			_Ttensym*g_yud= new _Ttensym(ypart.ny(),ypart.nys()+nu,d);
			g_yud->zeros();
			/*6:*/
			
			for(int i= 0;i<=d;i++){
				int j= d-i;
				int kfact= 1;
				_Ttensor tmp(ypart.ny(),
					TensorDimens(Symmetry(i,j),tns));
				tmp.zeros();
				for(int k= 0;k+d<=g.getMaxDim();k++,kfact*= k){
					Symmetry sym(i,j,0,k);
					if(g.check(sym)){
						double mult= pow(sigma,k)/dfact/kfact;
						tmp.add(mult,*(g.get(sym)));
					}
				}
				g_yud->addSubTensor(tmp);
			}
			
			/*:6*/
			;
			insert(g_yud);
		}
	}
	
	/*:5*/
	;
	/*7:*/
	
	void centralize(const DecisionRuleImpl&dr)
	{
		Vector dstate(ypart.nys()+nu);
		dstate.zeros();
		Vector dstate_star(dstate,0,ypart.nys());
		ConstVector newsteady_star(ysteady,ypart.nstat,ypart.nys());
		ConstVector oldsteady_star(dr.ysteady,ypart.nstat,ypart.nys());
		dstate_star.add(1.0,newsteady_star);
		dstate_star.add(-1.0,oldsteady_star);
		
		_Tpol pol(dr);
		for(int d= 1;d<=dr.getMaxDim();d++){
			pol.derivative(d-1);
			_Ttensym*der= pol.evalPartially(d,dstate);
			insert(der);
		}
	}
	
	/*:7*/
	;
	/*15:*/
	
	void eval(emethod em,Vector&out,const ConstVector&v)const
	{
		if(em==DecisionRule::horner)
			_Tparent::evalHorner(out,v);
		else
			_Tparent::evalTrad(out,v);
	}
	
	/*:15*/
	;
};

/*:4*/
;
/*17:*/

class UnfoldDecisionRule;
class FoldDecisionRule:public DecisionRuleImpl<KOrder::fold> {
	friend class UnfoldDecisionRule;
public:
	FoldDecisionRule(const ctraits<KOrder::fold> ::Tpol&pol,const PartitionY&yp,int nuu,
		const Vector&ys)
		:DecisionRuleImpl<KOrder::fold> (pol,yp,nuu,ys){}
	FoldDecisionRule(ctraits<KOrder::fold> ::Tpol&pol,const PartitionY&yp,int nuu,
		const Vector&ys)
		:DecisionRuleImpl<KOrder::fold> (pol,yp,nuu,ys){}
	FoldDecisionRule(const ctraits<KOrder::fold> ::Tg&g,const PartitionY&yp,int nuu,
		const Vector&ys,double sigma)
		:DecisionRuleImpl<KOrder::fold> (g,yp,nuu,ys,sigma){}
	FoldDecisionRule(const DecisionRuleImpl<KOrder::fold> &dr,const ConstVector&fixpoint)
		:DecisionRuleImpl<KOrder::fold> (dr,fixpoint){}
	FoldDecisionRule(const UnfoldDecisionRule&udr);
};

/*:17*/
;
/*18:*/

class UnfoldDecisionRule:public DecisionRuleImpl<KOrder::unfold> {
	friend class FoldDecisionRule;
public:
	UnfoldDecisionRule(const ctraits<KOrder::unfold> ::Tpol&pol,const PartitionY&yp,int nuu,
		const Vector&ys)
		:DecisionRuleImpl<KOrder::unfold> (pol,yp,nuu,ys){}
	UnfoldDecisionRule(ctraits<KOrder::unfold> ::Tpol&pol,const PartitionY&yp,int nuu,
		const Vector&ys)
		:DecisionRuleImpl<KOrder::unfold> (pol,yp,nuu,ys){}
	UnfoldDecisionRule(const ctraits<KOrder::unfold> ::Tg&g,const PartitionY&yp,int nuu,
		const Vector&ys,double sigma)
		:DecisionRuleImpl<KOrder::unfold> (g,yp,nuu,ys,sigma){}
	UnfoldDecisionRule(const DecisionRuleImpl<KOrder::unfold> &dr,const ConstVector&fixpoint)
		:DecisionRuleImpl<KOrder::unfold> (dr,fixpoint){}
	UnfoldDecisionRule(const FoldDecisionRule&udr);
};


/*:18*/
;
/*19:*/

template<int t> 
class DRFixPoint:public ctraits<t> ::Tpol{
	typedef typename ctraits<t> ::Tpol _Tparent;
	static int max_iter;
	static int max_newton_iter;
	static int newton_pause;
	static double tol;
	const Vector ysteady;
	const PartitionY ypart;
	_Tparent*bigf;
	_Tparent*bigfder;
public:
	typedef typename DecisionRule::emethod emethod;
	/*20:*/
	
	DRFixPoint(const _Tg&g,const PartitionY&yp,
		const Vector&ys,double sigma)
		:ctraits<t> ::Tpol(yp.ny(),yp.nys()),
		ysteady(ys),ypart(yp),bigf(NULL),bigfder(NULL)
	{
		fillTensors(g,sigma);
		_Tparent yspol(ypart.nstat,ypart.nys(),*this);
		bigf= new _Tparent((const _Tparent&)yspol);
		_Ttensym*frst= bigf->get(Symmetry(1));
		for(int i= 0;i<ypart.nys();i++)
			frst->get(i,i)= frst->get(i,i)-1;
		bigfder= new _Tparent(*bigf,0);
	}
	
	/*:20*/
	;
	/*21:*/
	
	virtual~DRFixPoint()
	{
		if(bigf)
			delete bigf;
		if(bigfder)
			delete bigfder;
	}
	
	/*:21*/
	;
	/*25:*/
	
	bool calcFixPoint(emethod em,Vector&out)
	{
		KORD_RAISE_IF(out.length()!=ypart.ny(),
			"Wrong length of out in DRFixPoint::calcFixPoint");
		
		Vector delta(ypart.nys());
		Vector ystar(ypart.nys());
		ystar.zeros();
		
		iter= 0;
		newton_iter_last= 0;
		newton_iter_total= 0;
		bool converged= false;
		do{
			if((iter/newton_pause)*newton_pause==iter)
				converged= solveNewton(ystar);
			if(!converged){
				bigf->evalHorner(delta,ystar);
				KORD_RAISE_IF_X(!delta.isFinite(),
					"NaN or Inf asserted in DRFixPoint::calcFixPoint",
					KORD_FP_NOT_FINITE);
				ystar.add(1.0,delta);
				converged= delta.getNorm()<tol;
			}
			iter++;
		}while(iter<max_iter&&!converged);
		
		if(converged){
			_Tparent::evalHorner(out,ystar);
			out.add(1.0,ysteady);
		}
		
		return converged;
	}
	
	
	/*:25*/
	;
	int getNumIter()const
	{return iter;}
	int getNewtonLastIter()const
	{return newton_iter_last;}
	int getNewtonTotalIter()const
	{return newton_iter_total;}
protected:
	/*22:*/
	
	void fillTensors(const _Tg&g,double sigma)
	{
		int dfact= 1;
		for(int d= 0;d<=g.getMaxDim();d++,dfact*= d){
			_Ttensym*g_yd= new _Ttensym(ypart.ny(),ypart.nys(),d);
			g_yd->zeros();
			int kfact= 1;
			for(int k= 0;d+k<=g.getMaxDim();k++,kfact*= k){
				if(g.check(Symmetry(d,0,0,k))){
					const _Ttensor*ten= g.get(Symmetry(d,0,0,k));
					double mult= pow(sigma,k)/dfact/kfact;
					g_yd->add(mult,*ten);
				}
			}
			insert(g_yd);
		}
	}
	
	/*:22*/
	;
	/*23:*/
	
	bool solveNewton(Vector&y)
	{
		const double urelax_threshold= 1.e-5;
		Vector sol((const Vector&)y);
		Vector delta(y.length());
		newton_iter_last= 0;
		bool delta_finite= true;
		double flastnorm= 0.0;
		double fnorm= 0.0;
		bool converged= false;
		double urelax= 1.0;
		
		do{
			_Ttensym*jacob= bigfder->evalPartially(1,sol);
			bigf->evalHorner(delta,sol);
			if(newton_iter_last==0)
				flastnorm= delta.getNorm();
			delta_finite= delta.isFinite();
			if(delta_finite){
				ConstTwoDMatrix(*jacob).multInvLeft(delta);
				/*24:*/
				
				bool urelax_found= false;
				urelax= 1.0;
				while(!urelax_found&&urelax> urelax_threshold){
					Vector soltmp((const Vector&)sol);
					soltmp.add(-urelax,delta);
					Vector f(sol.length());
					bigf->evalHorner(f,soltmp);
					fnorm= f.getNorm();
					if(fnorm<=flastnorm)
						urelax_found= true;
					else
						urelax*= std::min(0.5,flastnorm/fnorm);
				}
				
				
				/*:24*/
				;
				sol.add(-urelax,delta);
				delta_finite= delta.isFinite();
			}
			delete jacob;
			newton_iter_last++;
			converged= delta_finite&&fnorm<tol;
			flastnorm= fnorm;
		}while(!converged&&newton_iter_last<max_newton_iter
			&&urelax> urelax_threshold);
		
		newton_iter_total+= newton_iter_last;
		if(!converged)
			newton_iter_last= 0;
		y= (const Vector&)sol;
		return converged;
	}
	
	/*:23*/
	;
private:
	int iter;
	int newton_iter_last;
	int newton_iter_total;
};


/*:19*/
;
/*26:*/

class ExplicitShockRealization;
class SimResults{
protected:
	int num_y;
	int num_per;
	vector<TwoDMatrix*> data;
	vector<ExplicitShockRealization*> shocks;
public:
	SimResults(int ny,int nper)
		:num_y(ny),num_per(nper){}
	virtual~SimResults();
	void simulate(int num_sim,const DecisionRule&dr,const Vector&start,
		const TwoDMatrix&vcov,Journal&journal);
	void simulate(int num_sim,const DecisionRule&dr,const Vector&start,
		const TwoDMatrix&vcov);
	int getNumPer()const
	{return num_per;}
	int getNumSets()const
	{return(int)data.size();}
	const TwoDMatrix&getData(int i)const
	{return*(data[i]);}
	const ExplicitShockRealization&getShocks(int i)const
	{return*(shocks[i]);}
	bool addDataSet(TwoDMatrix*d,ExplicitShockRealization*sr);
	void writeMat4(const char*base,const char*lname)const;
	void writeMat4(FILE*fd,const char*lname)const;
};

/*:26*/
;
/*27:*/

class SimResultsStats:public SimResults{
protected:
	Vector mean;
	TwoDMatrix vcov;
public:
	SimResultsStats(int ny,int nper)
		:SimResults(ny,nper),mean(ny),vcov(ny,ny){}
	void simulate(int num_sim,const DecisionRule&dr,const Vector&start,
		const TwoDMatrix&vcov,Journal&journal);
	void writeMat4(FILE*fd,const char*lname)const;
protected:
	void calcMean();
	void calcVcov();
};

/*:27*/
;
/*28:*/

class SimResultsDynamicStats:public SimResults{
protected:
	TwoDMatrix mean;
	TwoDMatrix variance;
public:
	SimResultsDynamicStats(int ny,int nper)
		:SimResults(ny,nper),mean(ny,nper),variance(ny,nper){}
	void simulate(int num_sim,const DecisionRule&dr,const Vector&start,
		const TwoDMatrix&vcov,Journal&journal);
	void writeMat4(FILE*fd,const char*lname)const;
protected:
	void calcMean();
	void calcVariance();
};


/*:28*/
;
/*29:*/

class SimulationIRFWorker;
class SimResultsIRF:public SimResults{
	friend class SimulationIRFWorker;
protected:
	const SimResults&control;
	int ishock;
	double imp;
	TwoDMatrix means;
	TwoDMatrix variances;
public:
	SimResultsIRF(const SimResults&cntl,int ny,int nper,int i,double impulse)
		:SimResults(ny,nper),control(cntl),
		ishock(i),imp(impulse),
		means(ny,nper),variances(ny,nper){}
	void simulate(const DecisionRule&dr,const Vector&start,
		Journal&journal);
	void simulate(const DecisionRule&dr,const Vector&start);
	void writeMat4(FILE*fd,const char*lname)const;
protected:
	void calcMeans();
	void calcVariances();
};

/*:29*/
;
/*30:*/

class RTSimulationWorker;
class RTSimResultsStats{
	friend class RTSimulationWorker;
protected:
	Vector mean;
	TwoDMatrix vcov;
	int num_per;
	NormalConj nc;
	int incomplete_simulations;
	int thrown_periods;
public:
	RTSimResultsStats(int ny,int nper)
		:mean(ny),vcov(ny,ny),
		num_per(nper),nc(ny),
		incomplete_simulations(0),thrown_periods(0){}
	void simulate(int num_sim,const DecisionRule&dr,const Vector&start,
		const TwoDMatrix&vcov,Journal&journal);
	void simulate(int num_sim,const DecisionRule&dr,const Vector&start,
		const TwoDMatrix&vcov);
	void writeMat4(FILE*fd,const char*lname);
};

/*:30*/
;
/*31:*/

class DynamicModel;
class IRFResults{
	vector<SimResultsIRF*> irf_res;
	const DynamicModel&model;
	vector<int> irf_list_ind;
public:
	IRFResults(const DynamicModel&mod,const DecisionRule&dr,
		const SimResults&control,const vector<int> &ili,
		Journal&journal);
	~IRFResults();
	void writeMat4(FILE*fd,const char*prefix)const;
};

/*:31*/
;
/*32:*/

class SimulationWorker:public THREAD{
protected:
	SimResults&res;
	const DecisionRule&dr;
	DecisionRule::emethod em;
	int np;
	const Vector&st;
	ShockRealization&sr;
public:
	SimulationWorker(SimResults&sim_res,
		const DecisionRule&dec_rule,
		DecisionRule::emethod emet,int num_per,
		const Vector&start,ShockRealization&shock_r)
		:res(sim_res),dr(dec_rule),em(emet),np(num_per),st(start),sr(shock_r){}
	void operator()();
};

/*:32*/
;
/*33:*/

class SimulationIRFWorker:public THREAD{
	SimResultsIRF&res;
	const DecisionRule&dr;
	DecisionRule::emethod em;
	int np;
	const Vector&st;
	int idata;
	int ishock;
	double imp;
public:
	SimulationIRFWorker(SimResultsIRF&sim_res,
		const DecisionRule&dec_rule,
		DecisionRule::emethod emet,int num_per,
		const Vector&start,int id,
		int ishck,double impulse)
		:res(sim_res),dr(dec_rule),em(emet),np(num_per),st(start),
		idata(id),ishock(ishck),imp(impulse){}
	void operator()();
};

/*:33*/
;
/*34:*/

class RTSimulationWorker:public THREAD{
protected:
	RTSimResultsStats&res;
	const DecisionRule&dr;
	DecisionRule::emethod em;
	int np;
	const Vector&ystart;
	ShockRealization&sr;
public:
	RTSimulationWorker(RTSimResultsStats&sim_res,
		const DecisionRule&dec_rule,
		DecisionRule::emethod emet,int num_per,
		const Vector&start,ShockRealization&shock_r)
		:res(sim_res),dr(dec_rule),em(emet),np(num_per),ystart(start),sr(shock_r){}
	void operator()();
};

/*:34*/
;
/*35:*/

class RandomShockRealization:virtual public ShockRealization{
protected:
	MersenneTwister mtwister;
	TwoDMatrix factor;
public:
	RandomShockRealization(const TwoDMatrix&v,unsigned int iseed)
		:mtwister(iseed),factor(v.nrows(),v.nrows())
	{schurFactor(v);}
	RandomShockRealization(const RandomShockRealization&sr)
		:mtwister(sr.mtwister),factor(sr.factor){}
	virtual~RandomShockRealization(){}
	void get(int n,Vector&out);
	int numShocks()const
	{return factor.nrows();}
protected:
	void choleskyFactor(const TwoDMatrix&v);
	void schurFactor(const TwoDMatrix&v);
};

/*:35*/
;
/*36:*/

class ExplicitShockRealization:virtual public ShockRealization{
	TwoDMatrix shocks;
public:
	ExplicitShockRealization(const TwoDMatrix&sh)
		:shocks(sh){}
	ExplicitShockRealization(const ExplicitShockRealization&sr)
		:shocks(sr.shocks){}
	ExplicitShockRealization(ShockRealization&sr,int num_per);
	void get(int n,Vector&out);
	int numShocks()const
	{return shocks.nrows();}
	void addToShock(int ishock,int iper,double val);
	void print()const
	{shocks.print();}
};

/*:36*/
;
/*37:*/

class GenShockRealization:public RandomShockRealization,public ExplicitShockRealization{
public:
	GenShockRealization(const TwoDMatrix&v,const TwoDMatrix&sh,int seed)
		:RandomShockRealization(v,seed),ExplicitShockRealization(sh)
	{
		KORD_RAISE_IF(sh.nrows()!=v.nrows()||v.nrows()!=v.ncols(),
			"Wrong dimension of input matrix in GenShockRealization constructor");
	}
	void get(int n,Vector&out);
	int numShocks()const
	{return RandomShockRealization::numShocks();}
};

/*:37*/
;

#endif

/*:1*/
