#define _Ttensor TYPENAME ctraits<t> ::Ttensor
#define _Ttensym TYPENAME ctraits<t> ::Ttensym
#define _Tg TYPENAME ctraits<t> ::Tg
#define _Tgs TYPENAME ctraits<t> ::Tgs
#define _Tgss TYPENAME ctraits<t> ::Tgss
#define _TG TYPENAME ctraits<t> ::TG
#define _TZstack TYPENAME ctraits<t> ::TZstack
#define _TGstack TYPENAME ctraits<t> ::TGstack
#define _TZXstack TYPENAME ctraits<t> ::TZXstack
#define _TGXstack TYPENAME ctraits<t> ::TGXstack
#define _Tm TYPENAME ctraits<t> ::Tm
#define _Tpol TYPENAME ctraits<t> ::Tpol \
	\
	
/*1:*/

#ifndef KORDER_H
#define KORDER_H

#include "int_sequence.h"
#include "fs_tensor.h"
#include "gs_tensor.h"
#include "t_container.h"
#include "stack_container.h"
#include "normal_moments.h"
#include "t_polynomial.h"
#include "faa_di_bruno.h"
#include "journal.h"

#include "kord_exception.h"
#include "GeneralSylvester.h"

#include <cmath> 

#define TYPENAME typename

/*2:*/

class FoldedZXContainer;
class UnfoldedZXContainer;
class FoldedGXContainer;
class UnfoldedGXContainer;

template<bool condition,class Then,class Else> 
struct IF{
	typedef Then RET;
};

template<class Then,class Else> 
struct IF<false,Then,Else> {
	typedef Else RET;
};

template<int type> 
class ctraits{
public:
	enum{fold,unfold};
	typedef TYPENAME IF<type==fold,FGSTensor,UGSTensor> ::RET Ttensor;
	typedef TYPENAME IF<type==fold,FFSTensor,UFSTensor> ::RET Ttensym;
	typedef TYPENAME IF<type==fold,FGSContainer,UGSContainer> ::RET Tg;
	typedef TYPENAME IF<type==fold,FGSContainer,UGSContainer> ::RET Tgs;
	typedef TYPENAME IF<type==fold,FGSContainer,UGSContainer> ::RET Tgss;
	typedef TYPENAME IF<type==fold,FGSContainer,UGSContainer> ::RET TG;
	typedef TYPENAME IF<type==fold,FoldedZContainer,UnfoldedZContainer> ::RET TZstack;
	typedef TYPENAME IF<type==fold,FoldedGContainer,UnfoldedGContainer> ::RET TGstack;
	typedef TYPENAME IF<type==fold,FNormalMoments,UNormalMoments> ::RET Tm;
	typedef TYPENAME IF<type==fold,FTensorPolynomial,UTensorPolynomial> ::RET Tpol;
	typedef TYPENAME IF<type==fold,FoldedZXContainer,UnfoldedZXContainer> ::RET TZXstack;
	typedef TYPENAME IF<type==fold,FoldedGXContainer,UnfoldedGXContainer> ::RET TGXstack;
};


/*:2*/
;
/*3:*/

struct PartitionY{
	const int nstat;
	const int npred;
	const int nboth;
	const int nforw;
	PartitionY(int num_stat,int num_pred,
		int num_both,int num_forw)
		:nstat(num_stat),npred(num_pred),
		nboth(num_both),nforw(num_forw)
	{}
	int ny()const
	{return nstat+npred+nboth+nforw;}
	int nys()const
	{return npred+nboth;}
	int nyss()const
	{return nboth+nforw;}
};


/*:3*/
;
/*4:*/

class PLUMatrix:public TwoDMatrix{
public:
	PLUMatrix(int n)
		:TwoDMatrix(n,n),
		inv(nrows()*ncols()),
		ipiv(new int[nrows()]){}
	PLUMatrix(const PLUMatrix&plu);
	virtual~PLUMatrix()
	{delete[]ipiv;}
	void multInv(TwoDMatrix&m)const;
private:
	Vector inv;
	int*ipiv;
protected:
	void calcPLU();
};

/*:4*/
;
/*5:*/

class MatrixA:public PLUMatrix{
public:
	MatrixA(const FSSparseTensor&f,const IntSequence&ss,
		const TwoDMatrix&gy,const PartitionY&ypart);
};

/*:5*/
;
/*6:*/

class MatrixS:public PLUMatrix{
public:
	MatrixS(const FSSparseTensor&f,const IntSequence&ss,
		const TwoDMatrix&gy,const PartitionY&ypart);
};


/*:6*/
;
/*7:*/

class MatrixB:public TwoDMatrix{
public:
	MatrixB(const FSSparseTensor&f,const IntSequence&ss)
		:TwoDMatrix(FGSTensor(f,ss,IntSequence(1,0),
		TensorDimens(ss,IntSequence(1,0))))
	{}
};

/*:7*/
;
/*8:*/

class KOrder{
protected:
	const PartitionY ypart;
	const int ny;
	const int nu;
	const int maxk;
	IntSequence nvs;
	/*30:*/
	
	UGSContainer _ug;
	FGSContainer _fg;
	UGSContainer _ugs;
	FGSContainer _fgs;
	UGSContainer _ugss;
	FGSContainer _fgss;
	UGSContainer _uG;
	FGSContainer _fG;
	UnfoldedZContainer _uZstack;
	FoldedZContainer _fZstack;
	UnfoldedGContainer _uGstack;
	FoldedGContainer _fGstack;
	UNormalMoments _um;
	FNormalMoments _fm;
	const TensorContainer<FSSparseTensor> &f;
	
	/*:30*/
	;
	const MatrixA matA;
	const MatrixS matS;
	const MatrixB matB;
	/*31:*/
	
	template<int t> _Tg&g();
	template<int t> const _Tg&g()const;
	template<int t> _Tgs&gs();
	template<int t> const _Tgs&gs()const;
	template<int t> _Tgss&gss();
	template<int t> const _Tgss&gss()const;
	template<int t> _TG&G();
	template<int t> const _TG&G()const;
	template<int t> _TZstack&Zstack();
	template<int t> const _TZstack&Zstack()const;
	template<int t> _TGstack&Gstack();
	template<int t> const _TGstack&Gstack()const;
	template<int t> _Tm&m();
	template<int t> const _Tm&m()const;
	
	
	/*:31*/
	;
	Journal&journal;
public:
	KOrder(int num_stat,int num_pred,int num_both,int num_forw,
		const TensorContainer<FSSparseTensor> &fcont,
		const TwoDMatrix&gy,const TwoDMatrix&gu,const TwoDMatrix&v,
		Journal&jr);
	enum{fold,unfold};
	/*24:*/
	
	template<int t> 
		void performStep(int order)
	{
		KORD_RAISE_IF(order-1!=g<t> ().getMaxDim(),
			"Wrong order for KOrder::performStep");
		JournalRecordPair pa(journal);
		pa<<"Performing step for order = "<<order<<endrec;
		
		recover_y<t> (order);
		
		for(int i= 0;i<order;i++){
			recover_yu<t> (i,order-i);
		}
		
		for(int j= 1;j<order;j++){
			for(int i= j-1;i>=1;i--){
				recover_yus<t> (order-j,i,j-i);
			}
			recover_ys<t> (order-j,j);
		}
		
		for(int i= order-1;i>=1;i--){
			recover_yus<t> (0,i,order-i);
		}
		recover_s<t> (order);
	}
	
	/*:24*/
	;
	/*25:*/
	
	template<int t> 
		double check(int dim)const
	{
		KORD_RAISE_IF(dim> g<t> ().getMaxDim(),
			"Wrong dimension for KOrder::check");
		JournalRecordPair pa(journal);
		pa<<"Checking residuals for order = "<<dim<<endrec;
		
		double maxerror= 0.0;
		
		/*26:*/
		
		for(int i= 0;i<=dim;i++){
			Symmetry sym(dim-i,i,0,0);
			_Ttensor*r= faaDiBrunoZ<t> (sym);
			double err= r->getData().getMax();
			JournalRecord(journal)<<"\terror for symmetry "<<sym<<"\tis "<<err<<endrec;
			if(err> maxerror)
				maxerror= err;
			delete r;
		}
		
		/*:26*/
		;
		/*27:*/
		
		SymmetrySet ss(dim,3);
		for(symiterator si(ss);!si.isEnd();++si){
			int i= (*si)[0];
			int j= (*si)[1];
			int k= (*si)[2];
			if(i+j> 0&&k> 0){
				Symmetry sym(i,j,0,k);
				_Ttensor*r= faaDiBrunoZ<t> (sym);
				_Ttensor*D_ijk= calcD_ijk<t> (i,j,k);
				r->add(1.0,*D_ijk);
				delete D_ijk;
				_Ttensor*E_ijk= calcE_ijk<t> (i,j,k);
				r->add(1.0,*E_ijk);
				delete E_ijk;
				double err= r->getData().getMax();
				JournalRecord(journal)<<"\terror for symmetry "<<sym<<"\tis "<<err<<endrec;
				delete r;
			}
		}
		
		
		/*:27*/
		;
		/*28:*/
		
		_Ttensor*r= faaDiBrunoZ<t> (Symmetry(0,0,0,dim));
		_Ttensor*D_k= calcD_k<t> (dim);
		r->add(1.0,*D_k);
		delete D_k;
		_Ttensor*E_k= calcE_k<t> (dim);
		r->add(1.0,*E_k);
		delete E_k;
		double err= r->getData().getMax();
		Symmetry sym(0,0,0,dim);
		JournalRecord(journal)<<"\terror for symmetry "<<sym<<"\tis "<<err<<endrec;
		if(err> maxerror)
			maxerror= err;
		delete r;
		
		/*:28*/
		;
		
		return maxerror;
	}
	
	
	/*:25*/
	;
	/*29:*/
	
	template<int t> 
		Vector*calcStochShift(int order,double sigma)const
	{
		Vector*res= new Vector(ny);
		res->zeros();
		int jfac= 1;
		for(int j= 1;j<=order;j++,jfac*= j)
			if(is_even(j)){
				_Ttensor*ten= calcD_k<t> (j);
				res->add(std::pow(sigma,j)/jfac,ten->getData());
				delete ten;
			}
			return res;
	}
	
	
	/*:29*/
	;
	void switchToFolded();
	const PartitionY&getPartY()const
	{return ypart;}
	const FGSContainer&getFoldDers()const
	{return _fg;}
	const UGSContainer&getUnfoldDers()const
	{return _ug;}
	static bool is_even(int i)
	{return(i/2)*2==i;}
protected:
	/*9:*/
	
	template<int t> 
		void insertDerivative(_Ttensor*der)
	{
		g<t> ().insert(der);
		gs<t> ().insert(new _Ttensor(ypart.nstat,ypart.nys(),*der));
		gss<t> ().insert(new _Ttensor(ypart.nstat+ypart.npred,
			ypart.nyss(),*der));
	}
	
	
	/*:9*/
	;
	template<int t> 
		void sylvesterSolve(_Ttensor&der)const;
	
	/*10:*/
	
	template<int t> 
		_Ttensor*faaDiBrunoZ(const Symmetry&sym)const
	{
		JournalRecordPair pa(journal);
		pa<<"Faa Di Bruno Z container for "<<sym<<endrec;
		_Ttensor*res= new _Ttensor(ny,TensorDimens(sym,nvs));
		FaaDiBruno bruno(journal);
		bruno.calculate(Zstack<t> (),f,*res);
		return res;
	}
	
	/*:10*/
	;
	/*11:*/
	
	template<int t> 
		_Ttensor*faaDiBrunoG(const Symmetry&sym)const
	{
		JournalRecordPair pa(journal);
		pa<<"Faa Di Bruno G container for "<<sym<<endrec;
		TensorDimens tdims(sym,nvs);
		_Ttensor*res= new _Ttensor(ypart.nyss(),tdims);
		FaaDiBruno bruno(journal);
		bruno.calculate(Gstack<t> (),gss<t> (),*res);
		return res;
	}
	
	/*:11*/
	;
	
	/*12:*/
	
	template<int t> 
		void recover_y(int i)
	{
		Symmetry sym(i,0,0,0);
		JournalRecordPair pa(journal);
		pa<<"Recovering symmetry "<<sym<<endrec;
		
		_Ttensor*G_yi= faaDiBrunoG<t> (sym);
		G<t> ().insert(G_yi);
		
		_Ttensor*g_yi= faaDiBrunoZ<t> (sym);
		g_yi->mult(-1.0);
		
		sylvesterSolve<t> (*g_yi);
		
		insertDerivative<t> (g_yi);
		
		_Ttensor*gss_y= gss<t> ().get(Symmetry(1,0,0,0));
		gs<t> ().multAndAdd(*gss_y,*G_yi);
		_Ttensor*gss_yi= gss<t> ().get(sym);
		gs<t> ().multAndAdd(*gss_yi,*G_yi);
	}
	
	
	/*:12*/
	;
	/*13:*/
	
	template<int t> 
		void recover_yu(int i,int j)
	{
		Symmetry sym(i,j,0,0);
		JournalRecordPair pa(journal);
		pa<<"Recovering symmetry "<<sym<<endrec;
		
		_Ttensor*G_yiuj= faaDiBrunoG<t> (sym);
		G<t> ().insert(G_yiuj);
		
		_Ttensor*g_yiuj= faaDiBrunoZ<t> (sym);
		g_yiuj->mult(-1.0);
		matA.multInv(*g_yiuj);
		insertDerivative<t> (g_yiuj);
		
		gs<t> ().multAndAdd(*(gss<t> ().get(Symmetry(1,0,0,0))),*G_yiuj);
	}
	
	/*:13*/
	;
	/*14:*/
	
	template<int t> 
		void recover_ys(int i,int j)
	{
		Symmetry sym(i,0,0,j);
		JournalRecordPair pa(journal);
		pa<<"Recovering symmetry "<<sym<<endrec;
		
		fillG<t> (i,0,j);
		
		if(is_even(j)){
			_Ttensor*G_yisj= faaDiBrunoG<t> (sym);
			G<t> ().insert(G_yisj);
			
			_Ttensor*g_yisj= faaDiBrunoZ<t> (sym);
			
			{
				_Ttensor*D_ij= calcD_ik<t> (i,j);
				g_yisj->add(1.0,*D_ij);
				delete D_ij;
			}
			
			if(j>=3){
				_Ttensor*E_ij= calcE_ik<t> (i,j);
				g_yisj->add(1.0,*E_ij);
				delete E_ij;
			}
			
			g_yisj->mult(-1.0);
			
			sylvesterSolve<t> (*g_yisj);
			
			insertDerivative<t> (g_yisj);
			
			Gstack<t> ().multAndAdd(1,gss<t> (),*G_yisj);
			Gstack<t> ().multAndAdd(i+j,gss<t> (),*G_yisj);
		}
	}
	
	/*:14*/
	;
	/*15:*/
	
	template<int t> 
		void recover_yus(int i,int j,int k)
	{
		Symmetry sym(i,j,0,k);
		JournalRecordPair pa(journal);
		pa<<"Recovering symmetry "<<sym<<endrec;
		
		fillG<t> (i,j,k);
		
		if(is_even(k)){
			_Ttensor*G_yiujsk= faaDiBrunoG<t> (sym);
			G<t> ().insert(G_yiujsk);
			
			_Ttensor*g_yiujsk= faaDiBrunoZ<t> (sym);
			
			{
				_Ttensor*D_ijk= calcD_ijk<t> (i,j,k);
				g_yiujsk->add(1.0,*D_ijk);
				delete D_ijk;
			}
			
			if(k>=3){
				_Ttensor*E_ijk= calcE_ijk<t> (i,j,k);
				g_yiujsk->add(1.0,*E_ijk);
				delete E_ijk;
			}
			
			g_yiujsk->mult(-1.0);
			
			matA.multInv(*g_yiujsk);
			insertDerivative<t> (g_yiujsk);
			
			Gstack<t> ().multAndAdd(1,gss<t> (),*G_yiujsk);
		}
	}
	
	/*:15*/
	;
	/*16:*/
	
	template<int t> 
		void recover_s(int i)
	{
		Symmetry sym(0,0,0,i);
		JournalRecordPair pa(journal);
		pa<<"Recovering symmetry "<<sym<<endrec;
		
		fillG<t> (0,0,i);
		
		if(is_even(i)){
			_Ttensor*G_si= faaDiBrunoG<t> (sym);
			G<t> ().insert(G_si);
			
			_Ttensor*g_si= faaDiBrunoZ<t> (sym);
			
			{
				_Ttensor*D_i= calcD_k<t> (i);
				g_si->add(1.0,*D_i);
				delete D_i;
			}
			
			if(i>=3){
				_Ttensor*E_i= calcE_k<t> (i);
				g_si->add(1.0,*E_i);
				delete E_i;
			}
			
			g_si->mult(-1.0);
			
			
			matS.multInv(*g_si);
			insertDerivative<t> (g_si);
			
			Gstack<t> ().multAndAdd(1,gss<t> (),*G_si);
			Gstack<t> ().multAndAdd(i,gss<t> (),*G_si);
		}
	}
	
	/*:16*/
	;
	/*17:*/
	
	template<int t> 
		void fillG(int i,int j,int k)
	{
		for(int m= 1;m<=k;m++){
			if(is_even(k-m)){
				_Ttensor*G_yiujupms= faaDiBrunoG<t> (Symmetry(i,j,m,k-m));
				G<t> ().insert(G_yiujupms);
			}
		}
	}
	
	
	/*:17*/
	;
	
	/*18:*/
	
	template<int t> 
		_Ttensor*calcD_ijk(int i,int j,int k)const
	{
		_Ttensor*res= new _Ttensor(ny,TensorDimens(Symmetry(i,j,0,0),nvs));
		res->zeros();
		if(is_even(k)){
			_Ttensor*tmp= faaDiBrunoZ<t> (Symmetry(i,j,k,0));
			tmp->contractAndAdd(2,*res,*(m<t> ().get(Symmetry(k))));
			delete tmp;
		}
		return res;
	}
	
	
	/*:18*/
	;
	/*20:*/
	
	template<int t> 
		_Ttensor*calcD_ik(int i,int k)const
	{
		return calcD_ijk<t> (i,0,k);
	}
	
	/*:20*/
	;
	/*21:*/
	
	template<int t> 
		_Ttensor*calcD_k(int k)const
	{
		return calcD_ijk<t> (0,0,k);
	}
	
	/*:21*/
	;
	
	/*19:*/
	
	template<int t> 
		_Ttensor*calcE_ijk(int i,int j,int k)const
	{
		_Ttensor*res= new _Ttensor(ny,TensorDimens(Symmetry(i,j,0,0),nvs));
		res->zeros();
		for(int n= 2;n<=k-1;n+= 2){
			_Ttensor*tmp= faaDiBrunoZ<t> (Symmetry(i,j,n,k-n));
			tmp->mult((double)(Tensor::noverk(k,n)));
			tmp->contractAndAdd(2,*res,*(m<t> ().get(Symmetry(n))));
			delete tmp;
		}
		return res;
	}
	
	/*:19*/
	;
	/*22:*/
	
	template<int t> 
		_Ttensor*calcE_ik(int i,int k)const
	{
		return calcE_ijk<t> (i,0,k);
	}
	
	/*:22*/
	;
	/*23:*/
	
	template<int t> 
		_Ttensor*calcE_k(int k)const
	{
		return calcE_ijk<t> (0,0,k);
	}
	
	/*:23*/
	;
};



/*:8*/
;


#endif

/*:1*/
