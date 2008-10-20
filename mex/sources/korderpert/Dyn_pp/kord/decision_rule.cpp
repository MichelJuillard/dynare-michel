/*1:*/


#include "kord_exception.h"
#include "decision_rule.h"
#include "dynamic_model.h"

#include "SymSchurDecomp.h"
#include "cpplapack.h"

#include <limits> 

template<> 
int DRFixPoint<KOrder::fold> ::max_iter= 10000;
template<> 
int DRFixPoint<KOrder::unfold> ::max_iter= 10000;
template<> 
double DRFixPoint<KOrder::fold> ::tol= 1.e-10;
template<> 
double DRFixPoint<KOrder::unfold> ::tol= 1.e-10;
template<> 
int DRFixPoint<KOrder::fold> ::max_newton_iter= 50;
template<> 
int DRFixPoint<KOrder::unfold> ::max_newton_iter= 50;
template<> 
int DRFixPoint<KOrder::fold> ::newton_pause= 100;
template<> 
int DRFixPoint<KOrder::unfold> ::newton_pause= 100;

/*2:*/

FoldDecisionRule::FoldDecisionRule(const UnfoldDecisionRule&udr)
:DecisionRuleImpl<KOrder::fold> (ctraits<KOrder::fold> ::Tpol(udr.nrows(),udr.nvars()),
								 udr.ypart,udr.nu,udr.ysteady)
{
	for(ctraits<KOrder::unfold> ::Tpol::const_iterator it= udr.begin();
	it!=udr.end();++it){
		insert(new ctraits<KOrder::fold> ::Ttensym(*((*it).second)));
	}
}

/*:2*/
;
/*3:*/

UnfoldDecisionRule::UnfoldDecisionRule(const FoldDecisionRule&fdr)
:DecisionRuleImpl<KOrder::unfold> (ctraits<KOrder::unfold> ::Tpol(fdr.nrows(),fdr.nvars()),
								   fdr.ypart,fdr.nu,fdr.ysteady)
{
	for(ctraits<KOrder::fold> ::Tpol::const_iterator it= fdr.begin();
	it!=fdr.end();++it){
		insert(new ctraits<KOrder::unfold> ::Ttensym(*((*it).second)));
	}
}

/*:3*/
;
/*4:*/

SimResults::~SimResults()
{
	for(int i= 0;i<getNumSets();i++){
		delete data[i];
		delete shocks[i];
	}
}

/*:4*/
;
/*5:*/

void SimResults::simulate(int num_sim,const DecisionRule&dr,const Vector&start,
						  const TwoDMatrix&vcov,Journal&journal)
{
	JournalRecordPair paa(journal);
	paa<<"Performing "<<num_sim<<" stochastic simulations for "
		<<num_per<<" periods"<<endrec;
	simulate(num_sim,dr,start,vcov);
	int thrown= num_sim-data.size();
	if(thrown> 0){
		JournalRecord rec(journal);
		rec<<"I had to throw "<<thrown<<" simulations away due to Nan or Inf"<<endrec;
	}
}

/*:5*/
;
/*6:*/

void SimResults::simulate(int num_sim,const DecisionRule&dr,const Vector&start,
						  const TwoDMatrix&vcov)
{
	std::vector<RandomShockRealization> rsrs;
	rsrs.reserve(num_sim);
	
	THREAD_GROUP gr;
	for(int i= 0;i<num_sim;i++){
		RandomShockRealization sr(vcov,system_random_generator.int_uniform());
		rsrs.push_back(sr);
		THREAD*worker= new
			SimulationWorker(*this,dr,DecisionRule::horner,
			num_per,start,rsrs.back());
		gr.insert(worker);
	}
	gr.run();
}

/*:6*/
;
/*7:*/

bool SimResults::addDataSet(TwoDMatrix*d,ExplicitShockRealization*sr)
{
	KORD_RAISE_IF(d->nrows()!=num_y,
		"Incompatible number of rows for SimResults::addDataSets");
	KORD_RAISE_IF(d->ncols()!=num_per,
		"Incompatible number of cols for SimResults::addDataSets");
	if(d->isFinite()){
		data.push_back(d);
		shocks.push_back(sr);
		return true;
	}else{
		delete d;
		delete sr;
		return false;
	}
}

/*:7*/
;
/*8:*/

void SimResults::writeMat4(const char*base,const char*lname)const
{
	char matfile_name[100];
	sprintf(matfile_name,"%s.mat",base);
	FILE*out;
	if(NULL!=(out= fopen(matfile_name,"wb"))){
		writeMat4(out,lname);
		fclose(out);
	}
}

/*:8*/
;
/*9:*/

void SimResults::writeMat4(FILE*fd,const char*lname)const
{
	char tmp[100];
	for(int i= 0;i<getNumSets();i++){
		if(getNumSets()> 1)
			sprintf(tmp,"%s_data%d",lname,i+1);
		else
			sprintf(tmp,"%s_data",lname);
		ConstTwoDMatrix m(*(data[i]));
		m.writeMat4(fd,tmp);
	}
}

/*:9*/
;
/*10:*/

void SimResultsStats::simulate(int num_sim,const DecisionRule&dr,
							   const Vector&start,
							   const TwoDMatrix&vcov,Journal&journal)
{
	SimResults::simulate(num_sim,dr,start,vcov,journal);
	{
		JournalRecordPair paa(journal);
		paa<<"Calculating means from the simulations."<<endrec;
		calcMean();
	}
	{
		JournalRecordPair paa(journal);
		paa<<"Calculating covariances from the simulations."<<endrec;
		calcVcov();
	}
}


/*:10*/
;
/*11:*/

void SimResultsStats::writeMat4(FILE*fd,const char*lname)const
{
	char tmp[100];
	sprintf(tmp,"%s_mean",lname);
	ConstTwoDMatrix m(num_y,1,mean.base());
	m.writeMat4(fd,tmp);
	sprintf(tmp,"%s_vcov",lname);
	ConstTwoDMatrix(vcov).writeMat4(fd,tmp);
}

/*:11*/
;
/*12:*/

void SimResultsStats::calcMean()
{
	mean.zeros();
	if(data.size()*num_per> 0){
		double mult= 1.0/data.size()/num_per;
		for(unsigned int i= 0;i<data.size();i++){
			for(int j= 0;j<num_per;j++){
				ConstVector col(*data[i],j);
				mean.add(mult,col);
			}
		}
	}
}

/*:12*/
;
/*13:*/

void SimResultsStats::calcVcov()
{
	if(data.size()*num_per> 1){
		vcov.zeros();
		double mult= 1.0/(data.size()*num_per-1);
		for(unsigned int i= 0;i<data.size();i++){
			const TwoDMatrix&d= *(data[i]);
			for(int j= 0;j<num_per;j++){
				for(int m= 0;m<num_y;m++){
					for(int n= m;n<num_y;n++){
						double s= (d.get(m,j)-mean[m])*(d.get(n,j)-mean[n]);
						vcov.get(m,n)+= mult*s;
						if(m!=n)
							vcov.get(n,m)+= mult*s;
					}
				}
			}
		}
	}else{
		vcov.infs();
	}
}

/*:13*/
;
/*14:*/

void SimResultsDynamicStats::simulate(int num_sim,const DecisionRule&dr,
									  const Vector&start,
									  const TwoDMatrix&vcov,Journal&journal)
{
	SimResults::simulate(num_sim,dr,start,vcov,journal);
	{
		JournalRecordPair paa(journal);
		paa<<"Calculating means of the conditional simulations."<<endrec;
		calcMean();
	}
	{
		JournalRecordPair paa(journal);
		paa<<"Calculating variances of the conditional simulations."<<endrec;
		calcVariance();
	}
}

/*:14*/
;
/*15:*/

void SimResultsDynamicStats::writeMat4(FILE*fd,const char*lname)const
{
	char tmp[100];
	sprintf(tmp,"%s_cond_mean",lname);
	ConstTwoDMatrix(mean).writeMat4(fd,tmp);
	sprintf(tmp,"%s_cond_variance",lname);
	ConstTwoDMatrix(variance).writeMat4(fd,tmp);
}

/*:15*/
;
/*16:*/

void SimResultsDynamicStats::calcMean()
{
	mean.zeros();
	if(data.size()> 0){
		double mult= 1.0/data.size();
		for(int j= 0;j<num_per;j++){
			Vector meanj(mean,j);
			for(unsigned int i= 0;i<data.size();i++){
				ConstVector col(*data[i],j);
				meanj.add(mult,col);
			}
		}
	}
}

/*:16*/
;
/*17:*/

void SimResultsDynamicStats::calcVariance()
{
	if(data.size()> 1){
		variance.zeros();
		double mult= 1.0/(data.size()-1);
		for(int j= 0;j<num_per;j++){
			ConstVector meanj(mean,j);
			Vector varj(variance,j);
			for(int i= 0;i<(int)data.size();i++){
				Vector col(ConstVector((*data[i]),j));
				col.add(-1.0,meanj);
				for(int k= 0;k<col.length();k++)
					col[k]= col[k]*col[k];
				varj.add(mult,col);
			}
		}
	}else{
		variance.infs();
	}
}


/*:17*/
;
/*18:*/

void SimResultsIRF::simulate(const DecisionRule&dr,const Vector&start,
							 Journal&journal)
{
	JournalRecordPair paa(journal);
	paa<<"Performing "<<control.getNumSets()<<" IRF simulations for "
		<<num_per<<" periods; shock="<<ishock<<", impulse="<<imp<<endrec;
	simulate(dr,start);
	int thrown= control.getNumSets()-data.size();
	if(thrown> 0){
		JournalRecord rec(journal);
		rec<<"I had to throw "<<thrown
			<<" simulations away due to Nan or Inf"<<endrec;
	}
	calcMeans();
	calcVariances();
}

/*:18*/
;
/*19:*/

void SimResultsIRF::simulate(const DecisionRule&dr,const Vector&start)
{
	THREAD_GROUP gr;
	for(int idata= 0;idata<control.getNumSets();idata++){
		THREAD*worker= new
			SimulationIRFWorker(*this,dr,DecisionRule::horner,
			num_per,start,idata,ishock,imp);
		gr.insert(worker);
	}
	gr.run();
}

/*:19*/
;
/*20:*/

void SimResultsIRF::calcMeans()
{
	means.zeros();
	if(data.size()> 0){
		for(unsigned int i= 0;i<data.size();i++)
			means.add(1.0,*(data[i]));
		means.mult(1.0/data.size());
	}
}

/*:20*/
;
/*21:*/

void SimResultsIRF::calcVariances()
{
	if(data.size()> 1){
		variances.zeros();
		for(unsigned int i= 0;i<data.size();i++){
			TwoDMatrix d((const TwoDMatrix&)(*(data[i])));
			d.add(-1.0,means);
			for(int j= 0;j<d.nrows();j++)
				for(int k= 0;k<d.ncols();k++)
					variances.get(j,k)+= d.get(j,k)*d.get(j,k);
				d.mult(1.0/(data.size()-1));
		}
	}else{
		variances.infs();
	}
}

/*:21*/
;
/*22:*/

void SimResultsIRF::writeMat4(FILE*fd,const char*lname)const
{
	char tmp[100];
	sprintf(tmp,"%s_mean",lname);
	means.writeMat4(fd,tmp);
	sprintf(tmp,"%s_var",lname);
	variances.writeMat4(fd,tmp);
}

/*:22*/
;
/*23:*/

void RTSimResultsStats::simulate(int num_sim,const DecisionRule&dr,const Vector&start,
								 const TwoDMatrix&v,Journal&journal)
{
	JournalRecordPair paa(journal);
	paa<<"Performing "<<num_sim<<" real-time stochastic simulations for "
		<<num_per<<" periods"<<endrec;
	simulate(num_sim,dr,start,v);
	mean= nc.getMean();
	mean.add(1.0,dr.getSteady());
	nc.getVariance(vcov);
	if(thrown_periods> 0){
		JournalRecord rec(journal);
		rec<<"I had to throw "<<thrown_periods<<" periods away due to Nan or Inf"<<endrec;
		JournalRecord rec1(journal);
		rec1<<"This affected "<<incomplete_simulations<<" out of "
			<<num_sim<<" simulations"<<endrec;
	}
}

/*:23*/
;
/*24:*/

void RTSimResultsStats::simulate(int num_sim,const DecisionRule&dr,const Vector&start,
								 const TwoDMatrix&vcov)
{
	std::vector<RandomShockRealization> rsrs;
	rsrs.reserve(num_sim);
	
	THREAD_GROUP gr;
	for(int i= 0;i<num_sim;i++){
		RandomShockRealization sr(vcov,system_random_generator.int_uniform());
		rsrs.push_back(sr);
		THREAD*worker= new
			RTSimulationWorker(*this,dr,DecisionRule::horner,
			num_per,start,rsrs.back());
		gr.insert(worker);
	}
	gr.run();
}

/*:24*/
;
/*25:*/

void RTSimResultsStats::writeMat4(FILE*fd,const char*lname)
{
	char tmp[100];
	sprintf(tmp,"%s_rt_mean",lname);
	ConstTwoDMatrix m(nc.getDim(),1,mean.base());
	m.writeMat4(fd,tmp);
	sprintf(tmp,"%s_rt_vcov",lname);
	ConstTwoDMatrix(vcov).writeMat4(fd,tmp);
}

/*:25*/
;
/*26:*/

IRFResults::IRFResults(const DynamicModel&mod,const DecisionRule&dr,
					   const SimResults&control,const vector<int> &ili,
					   Journal&journal)
					   :model(mod),irf_list_ind(ili)
{
	int num_per= control.getNumPer();
	JournalRecordPair pa(journal);
	pa<<"Calculating IRFs against control for "<<(int)irf_list_ind.size()<<" shocks and for "
		<<num_per<<" periods"<<endrec;
	const TwoDMatrix&vcov= mod.getVcov();
	for(unsigned int ii= 0;ii<irf_list_ind.size();ii++){
		int ishock= irf_list_ind[ii];
		double stderror= sqrt(vcov.get(ishock,ishock));
		irf_res.push_back(new SimResultsIRF(control,model.numeq(),num_per,
			ishock,stderror));
		irf_res.push_back(new SimResultsIRF(control,model.numeq(),num_per,
			ishock,-stderror));
	}
	
	for(unsigned int ii= 0;ii<irf_list_ind.size();ii++){
		irf_res[2*ii]->simulate(dr,model.getSteady(),journal);
		irf_res[2*ii+1]->simulate(dr,model.getSteady(),journal);
	}
}

/*:26*/
;
/*27:*/

IRFResults::~IRFResults()
{
	for(unsigned int i= 0;i<irf_res.size();i++)
		delete irf_res[i];
}

/*:27*/
;
/*28:*/

void IRFResults::writeMat4(FILE*fd,const char*prefix)const
{
	for(unsigned int i= 0;i<irf_list_ind.size();i++){
		char tmp[100];
		int ishock= irf_list_ind[i];
		const char*shockname= model.getExogNames().getName(ishock);
		sprintf(tmp,"%s_irfp_%s",prefix,shockname);
		irf_res[2*i]->writeMat4(fd,tmp);
		sprintf(tmp,"%s_irfm_%s",prefix,shockname);
		irf_res[2*i+1]->writeMat4(fd,tmp);
	}
}

/*:28*/
;
/*29:*/

void SimulationWorker::operator()()
{
	ExplicitShockRealization*esr= new ExplicitShockRealization(sr,np);
	TwoDMatrix*m= dr.simulate(em,np,st,*esr);
	{
		SYNCHRO syn(&res,"simulation");
		res.addDataSet(m,esr);
	}
}

/*:29*/
;
/*30:*/

void SimulationIRFWorker::operator()()
{
	ExplicitShockRealization*esr= 
		new ExplicitShockRealization(res.control.getShocks(idata));
	esr->addToShock(ishock,0,imp);
	TwoDMatrix*m= dr.simulate(em,np,st,*esr);
	m->add(-1.0,res.control.getData(idata));
	{
		SYNCHRO syn(&res,"simulation");
		res.addDataSet(m,esr);
	}
}

/*:30*/
;
/*31:*/

void RTSimulationWorker::operator()()
{
	NormalConj nc(res.nc.getDim());
	const PartitionY&ypart= dr.getYPart();
	int nu= dr.nexog();
	const Vector&ysteady= dr.getSteady();
	
	/*32:*/
	
	Vector dyu(ypart.nys()+nu);
	ConstVector ystart_pred(ystart,ypart.nstat,ypart.nys());
	ConstVector ysteady_pred(ysteady,ypart.nstat,ypart.nys());
	Vector dy(dyu,0,ypart.nys());
	Vector u(dyu,ypart.nys(),nu);
	Vector y(nc.getDim());
	ConstVector ypred(y,ypart.nstat,ypart.nys());
	
	/*:32*/
	;
	/*33:*/
	
	int ip= 0;
	dy= ystart_pred;
	dy.add(-1.0,ysteady_pred);
	sr.get(ip,u);
	dr.eval(em,y,dyu);
	nc.update(y);
	
	/*:33*/
	;
	/*34:*/
	
	while(y.isFinite()&&ip<res.num_per){
		ip++;
		dy= ypred;
		sr.get(ip,u);
		dr.eval(em,y,dyu);
		nc.update(y);
	}
	
	/*:34*/
	;
	{
		SYNCHRO syn(&res,"rtsimulation");
		res.nc.update(nc);
		if(res.num_per-ip> 0){
			res.incomplete_simulations++;
			res.thrown_periods+= res.num_per-ip;
		}
	}
}

/*:31*/
;
/*35:*/

void RandomShockRealization::choleskyFactor(const TwoDMatrix&v)
{
	factor= v;
	int rows= factor.nrows();
	for(int i= 0;i<rows;i++)
		for(int j= i+1;j<rows;j++)
			factor.get(i,j)= 0.0;
		int info;
		
		LAPACK_dpotrf("L",&rows,factor.base(),&rows,&info);
		KORD_RAISE_IF(info!=0,
			"Info!=0 in RandomShockRealization::choleskyFactor");
}

/*:35*/
;
/*36:*/

void RandomShockRealization::schurFactor(const TwoDMatrix&v)
{
	SymSchurDecomp ssd(v);
	ssd.getFactor(factor);
}

/*:36*/
;
/*37:*/

void RandomShockRealization::get(int n,Vector&out)
{
	KORD_RAISE_IF(out.length()!=numShocks(),
		"Wrong length of out vector in RandomShockRealization::get");
	Vector d(out.length());
	for(int i= 0;i<d.length();i++){
		d[i]= mtwister.normal();
	}
	out.zeros();
	factor.multaVec(out,ConstVector(d));
}

/*:37*/
;
/*38:*/

ExplicitShockRealization::ExplicitShockRealization(ShockRealization&sr,
												   int num_per)
												   :shocks(sr.numShocks(),num_per)
{
	for(int j= 0;j<num_per;j++){
		Vector jcol(shocks,j);
		sr.get(j,jcol);
	}
}

/*:38*/
;
/*39:*/

void ExplicitShockRealization::get(int n,Vector&out)
{
	KORD_RAISE_IF(out.length()!=numShocks(),
		"Wrong length of out vector in ExplicitShockRealization::get");
	int i= n%shocks.ncols();
	ConstVector icol(shocks,i);
	out= icol;
}

/*:39*/
;
/*40:*/

void ExplicitShockRealization::addToShock(int ishock,int iper,double val)
{
	KORD_RAISE_IF(ishock<0||ishock> numShocks(),
		"Wrong index of shock in ExplicitShockRealization::addToShock");
	int j= iper%shocks.ncols();
	shocks.get(ishock,j)+= val;
}


/*:40*/
;
/*41:*/

void GenShockRealization::get(int n,Vector&out)
{
	KORD_RAISE_IF(out.length()!=numShocks(),
		"Wrong length of out vector in GenShockRealization::get");
	ExplicitShockRealization::get(n,out);
	Vector r(numShocks());
	RandomShockRealization::get(n,r);
	for(int j= 0;j<numShocks();j++)
		if(!isfinite(out[j]))
			out[j]= r[j];
}


/*:41*/
;

/*:1*/
