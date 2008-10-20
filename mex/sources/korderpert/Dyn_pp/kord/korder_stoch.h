/*1:*/

#include "korder.h"
#include "faa_di_bruno.h"
#include "journal.h"


/*2:*/

template<int t> 
class IntegDerivs:public ctraits<t> ::Tgss{
public:
	/*3:*/
	
	IntegDerivs(int r,const IntSequence&nvs,const _Tgss&g,const _Tm&mom,
		double at_sigma)
		:ctraits<t> ::Tgss(4)
	{
		int maxd= g.getMaxDim();
		for(int d= 1;d<=maxd;d++){
			for(int i= 0;i<=d;i++){
				int p= d-i;
				Symmetry sym(i,0,0,p);
				_Ttensor*ten= new _Ttensor(r,TensorDimens(sym,nvs));
				/*4:*/
				
				ten->zeros();
				for(int n= 0;n<=p;n++){
					int k= p-n;
					int povern= Tensor::noverk(p,n);
					int mfac= 1;
					for(int m= 0;i+m+n+k<=maxd;m++,mfac*= m){
						double mult= (pow(at_sigma,m)*povern)/mfac;
						Symmetry sym_mn(i,m+n,0,k);
						if(m+n==0&&g.check(sym_mn))
							ten->add(mult,*(g.get(sym_mn)));
						if(m+n> 0&&KOrder::is_even(m+n)&&g.check(sym_mn)){
							_Ttensor gtmp(*(g.get(sym_mn)));
							gtmp.mult(mult);
							gtmp.contractAndAdd(1,*ten,*(mom.get(Symmetry(m+n))));
						}
					}
				}
				
				/*:4*/
				;
				insert(ten);
			}
		}
	}
	
	/*:3*/
	;
};

/*:2*/
;
/*5:*/

template<int t> 
class StochForwardDerivs:public ctraits<t> ::Tgss{
public:
	/*6:*/
	
	StochForwardDerivs(const PartitionY&ypart,int nu,
		const _Tgss&g,const _Tm&m,
		const Vector&ydelta,double sdelta,
		double at_sigma)
		:ctraits<t> ::Tgss(4)
	{
		int maxd= g.getMaxDim();
		int r= ypart.nyss();
		
		/*7:*/
		
		IntSequence nvs(4);
		nvs[0]= ypart.nys();nvs[1]= 0;nvs[2]= 0;nvs[3]= 1;
		IntegDerivs<t> g_int(r,nvs,g,m,at_sigma);
		
		/*:7*/
		;
		/*8:*/
		
		_Tpol g_int_sym(r,ypart.nys()+1);
		for(int d= 1;d<=maxd;d++){
			_Ttensym*ten= new _Ttensym(r,ypart.nys()+1,d);
			ten->zeros();
			for(int i= 0;i<=d;i++){
				int k= d-i;
				if(g_int.check(Symmetry(i,0,0,k)))
					ten->addSubTensor(*(g_int.get(Symmetry(i,0,0,k))));
			}
			g_int_sym.insert(ten);
		}
		
		/*:8*/
		;
		/*9:*/
		
		Vector delta(ypart.nys()+1);
		Vector dy(delta,0,ypart.nys());
		ConstVector dy_in(ydelta,ypart.nstat,ypart.nys());
		dy= dy_in;
		delta[ypart.nys()]= sdelta;
		_Tpol g_int_cent(r,ypart.nys()+1);
		for(int d= 1;d<=maxd;d++){
			g_int_sym.derivative(d-1);
			_Ttensym*der= g_int_sym.evalPartially(d,delta);
			g_int_cent.insert(der);
		}
		
		/*:9*/
		;
		/*10:*/
		
		IntSequence ss(4);
		ss[0]= ypart.nys();ss[1]= 0;ss[2]= 0;ss[3]= 1;
		IntSequence pp(4);
		pp[0]= 0;pp[1]= 1;pp[2]= 2;pp[3]= 3;
		IntSequence true_nvs(nvs);
		true_nvs[1]= nu;true_nvs[2]= nu;
		for(int d= 1;d<=maxd;d++){
			if(g_int_cent.check(Symmetry(d))){
				for(int i= 0;i<=d;i++){
					Symmetry sym(i,0,0,d-i);
					IntSequence coor(sym,pp);
					_Ttensor*ten= new _Ttensor(*(g_int_cent.get(Symmetry(d))),ss,coor,
						TensorDimens(sym,true_nvs));
					insert(ten);
				}
			}
		}
		
		
		/*:10*/
		;
	}
	
	/*:6*/
	;
};

/*:5*/
;
/*11:*/

template<class _Ttype> 
class GXContainer:public GContainer<_Ttype> {
public:
	typedef StackContainerInterface<_Ttype> _Stype;
	typedef typename StackContainer<_Ttype> ::_Ctype _Ctype;
	typedef typename StackContainer<_Ttype> ::itype itype;
	GXContainer(const _Ctype*gs,int ngs,int nu)
		:GContainer<_Ttype> (gs,ngs,nu){}
	/*12:*/
	
	itype getType(int i,const Symmetry&s)const
	{
		if(i==0)
			if(s[2]> 0)
				return _Stype::zero;
			else
				return _Stype::matrix;
			if(i==1)
				return _Stype::zero;
			if(i==2)
				return _Stype::zero;
			if(i==3)
				if(s==Symmetry(0,0,0,1))
					return _Stype::unit;
				else
					return _Stype::zero;
				
				KORD_RAISE("Wrong stack index in GXContainer::getType");
	}
	
	
	/*:12*/
	;
};

/*:11*/
;
/*13:*/

template<class _Ttype> 
class ZXContainer:public ZContainer<_Ttype> {
public:
	typedef StackContainerInterface<_Ttype> _Stype;
	typedef typename StackContainer<_Ttype> ::_Ctype _Ctype;
	typedef typename StackContainer<_Ttype> ::itype itype;
	ZXContainer(const _Ctype*gss,int ngss,const _Ctype*g,int ng,int ny,int nu)
		:ZContainer<_Ttype> (gss,ngss,g,ng,ny,nu){}
	/*14:*/
	
	itype getType(int i,const Symmetry&s)const
	{
		if(i==0)
			if(s[2]> 0)
				return _Stype::zero;
			else
				return _Stype::matrix;
			if(i==1)
				if(s[2]> 0)
					return _Stype::zero;
				else
					return _Stype::matrix;
				if(i==2)
					if(s==Symmetry(1,0,0,0))
						return _Stype::unit;
					else
						return _Stype::zero;
					if(i==3)
						if(s==Symmetry(0,1,0,0))
							return _Stype::unit;
						else
							return _Stype::zero;
						
						KORD_RAISE("Wrong stack index in ZXContainer::getType");
	}
	
	/*:14*/
	;
};

/*:13*/
;
/*15:*/

class UnfoldedGXContainer:public GXContainer<UGSTensor> ,public UnfoldedStackContainer{
public:
	typedef TensorContainer<UGSTensor> _Ctype;
	UnfoldedGXContainer(const _Ctype*gs,int ngs,int nu)
		:GXContainer<UGSTensor> (gs,ngs,nu){}
};

/*:15*/
;
/*16:*/

class FoldedGXContainer:public GXContainer<FGSTensor> ,public FoldedStackContainer{
public:
	typedef TensorContainer<FGSTensor> _Ctype;
	FoldedGXContainer(const _Ctype*gs,int ngs,int nu)
		:GXContainer<FGSTensor> (gs,ngs,nu){}
};

/*:16*/
;
/*17:*/

class UnfoldedZXContainer:public ZXContainer<UGSTensor> ,public UnfoldedStackContainer{
public:
	typedef TensorContainer<UGSTensor> _Ctype;
	UnfoldedZXContainer(const _Ctype*gss,int ngss,const _Ctype*g,int ng,int ny,int nu)
		:ZXContainer<UGSTensor> (gss,ngss,g,ng,ny,nu){}
};

/*:17*/
;
/*18:*/

class FoldedZXContainer:public ZXContainer<FGSTensor> ,public FoldedStackContainer{
public:
	typedef TensorContainer<FGSTensor> _Ctype;
	FoldedZXContainer(const _Ctype*gss,int ngss,const _Ctype*g,int ng,int ny,int nu)
		:ZXContainer<FGSTensor> (gss,ngss,g,ng,ny,nu){}
};

/*:18*/
;
/*19:*/

class MatrixAA:public PLUMatrix{
public:
	MatrixAA(const FSSparseTensor&f,const IntSequence&ss,
		const TwoDMatrix&gyss,const PartitionY&ypart);
};


/*:19*/
;
/*20:*/

class KOrderStoch{
protected:
	IntSequence nvs;
	PartitionY ypart;
	Journal&journal;
	UGSContainer _ug;
	FGSContainer _fg;
	UGSContainer _ugs;
	FGSContainer _fgs;
	UGSContainer _uG;
	FGSContainer _fG;
	const UGSContainer*_uh;
	const FGSContainer*_fh;
	UnfoldedZXContainer _uZstack;
	FoldedZXContainer _fZstack;
	UnfoldedGXContainer _uGstack;
	FoldedGXContainer _fGstack;
	const TensorContainer<FSSparseTensor> &f;
	MatrixAA matA;
public:
	KOrderStoch(const PartitionY&ypart,int nu,const TensorContainer<FSSparseTensor> &fcont,
		const FGSContainer&hh,Journal&jr);
	KOrderStoch(const PartitionY&ypart,int nu,const TensorContainer<FSSparseTensor> &fcont,
		const UGSContainer&hh,Journal&jr);
	/*23:*/
	
	template<int t> 
		void performStep(int order)
	{
		int maxd= g<t> ().getMaxDim();
		KORD_RAISE_IF(order-1!=maxd&&(order!=1||maxd!=-1),
			"Wrong order for KOrderStoch::performStep");
		SymmetrySet ss(order,4);
		for(symiterator si(ss);!si.isEnd();++si){
			if((*si)[2]==0){
				JournalRecordPair pa(journal);
				pa<<"Recovering symmetry "<<*si<<endrec;
				
				_Ttensor*G_sym= faaDiBrunoG<t> (*si);
				G<t> ().insert(G_sym);
				
				_Ttensor*g_sym= faaDiBrunoZ<t> (*si);
				g_sym->mult(-1.0);
				matA.multInv(*g_sym);
				g<t> ().insert(g_sym);
				gs<t> ().insert(new _Ttensor(ypart.nstat,ypart.nys(),*g_sym));
				
				Gstack<t> ().multAndAdd(1,h<t> (),*G_sym);
			}
		}
	}
	
	/*:23*/
	;
	const FGSContainer&getFoldDers()const
	{return _fg;}
	const UGSContainer&getUnfoldDers()const
	{return _ug;}
protected:
	/*21:*/
	
	template<int t> 
		_Ttensor*faaDiBrunoZ(const Symmetry&sym)const
	{
		JournalRecordPair pa(journal);
		pa<<"Faa Di Bruno ZX container for "<<sym<<endrec;
		_Ttensor*res= new _Ttensor(ypart.ny(),TensorDimens(sym,nvs));
		FaaDiBruno bruno(journal);
		bruno.calculate(Zstack<t> (),f,*res);
		return res;
	}
	
	/*:21*/
	;
	/*22:*/
	
	template<int t> 
		_Ttensor*faaDiBrunoG(const Symmetry&sym)const
	{
		JournalRecordPair pa(journal);
		pa<<"Faa Di Bruno GX container for "<<sym<<endrec;
		TensorDimens tdims(sym,nvs);
		_Ttensor*res= new _Ttensor(ypart.nyss(),tdims);
		FaaDiBruno bruno(journal);
		bruno.calculate(Gstack<t> (),h<t> (),*res);
		return res;
	}
	
	/*:22*/
	;
	/*24:*/
	
	template<int t> _Tg&g();
	template<int t> const _Tg&g()const;
	template<int t> _Tgs&gs();
	template<int t> const _Tgs&gs()const;
	template<int t> const _Tgss&h()const;
	template<int t> _TG&G();
	template<int t> const _TG&G()const;
	template<int t> _TZXstack&Zstack();
	template<int t> const _TZXstack&Zstack()const;
	template<int t> _TGXstack&Gstack();
	template<int t> const _TGXstack&Gstack()const;
	
	
	/*:24*/
	;
};

/*:20*/
;

/*:1*/
