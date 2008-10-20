/*1:*/

#include "faa_di_bruno.h"
#include "fine_container.h"

#include <math.h> 

double FaaDiBruno::magic_mult= 1.5;
/*2:*/

void FaaDiBruno::calculate(const StackContainer<FGSTensor> &cont,
						   const TensorContainer<FSSparseTensor> &f,
						   FGSTensor&out)
{
	out.zeros();
	for(int l= 1;l<=out.dimen();l++){
		int mem_mb,p_size_mb;
		int max= estimRefinment(out.getDims(),out.nrows(),l,mem_mb,p_size_mb);
		FoldedFineContainer fine_cont(cont,max);
		fine_cont.multAndAdd(l,f,out);
		JournalRecord recc(journal);
		recc<<"dim="<<l<<" avmem="<<mem_mb<<" tmpmem="<<p_size_mb<<" max="<<max
			<<" stacks="<<cont.numStacks()<<"->"<<fine_cont.numStacks()<<endrec;
	}
}

/*:2*/
;
/*3:*/

void FaaDiBruno::calculate(const FoldedStackContainer&cont,const FGSContainer&g,
						   FGSTensor&out)
{
	out.zeros();
	for(int l= 1;l<=out.dimen();l++){
		long int mem= SystemResources::availableMemory();
		cont.multAndAdd(l,g,out);
		JournalRecord rec(journal);
		int mem_mb= mem/1024/1024;
		rec<<"dim="<<l<<" avmem="<<mem_mb<<endrec;
	}
}

/*:3*/
;
/*4:*/

void FaaDiBruno::calculate(const StackContainer<UGSTensor> &cont,
						   const TensorContainer<FSSparseTensor> &f,
						   UGSTensor&out)
{
	out.zeros();
	for(int l= 1;l<=out.dimen();l++){
		int mem_mb,p_size_mb;
		int max= estimRefinment(out.getDims(),out.nrows(),l,mem_mb,p_size_mb);
		UnfoldedFineContainer fine_cont(cont,max);
		fine_cont.multAndAdd(l,f,out);
		JournalRecord recc(journal);
		recc<<"dim="<<l<<" avmem="<<mem_mb<<" tmpmem="<<p_size_mb<<" max="<<max
			<<" stacks="<<cont.numStacks()<<"->"<<fine_cont.numStacks()<<endrec;
	}
}

/*:4*/
;
/*5:*/

void FaaDiBruno::calculate(const UnfoldedStackContainer&cont,const UGSContainer&g,
						   UGSTensor&out)
{
	out.zeros();
	for(int l= 1;l<=out.dimen();l++){
		long int mem= SystemResources::availableMemory();
		cont.multAndAdd(l,g,out);
		JournalRecord rec(journal);
		int mem_mb= mem/1024/1024;
		rec<<"dim="<<l<<" avmem="<<mem_mb<<endrec;
	}
}

/*:5*/
;
/*6:*/

int FaaDiBruno::estimRefinment(const TensorDimens&tdims,int nr,int l,
							   int&avmem_mb,int&tmpmem_mb)
{
	int nthreads= THREAD_GROUP::max_parallel_threads;
	long int per_size1= tdims.calcUnfoldMaxOffset();
	long int per_size2= (long int)pow((double)tdims.getNVS().getMax(),l);
	double lambda= 0.0;
	long int per_size= sizeof(double)*nr
		*(long int)(lambda*per_size1+(1-lambda)*per_size2);
	long int mem= SystemResources::availableMemory();
	int max= 0;
	double num_cols= ((double)(mem-magic_mult*nthreads*per_size))
		/nthreads/sizeof(double)/nr;
	if(num_cols> 0){
		double maxd= pow(num_cols,((double)1)/l);
		max= (int)floor(maxd);
	}
	if(max==0){
		max= 10;
		JournalRecord rec(journal);
		rec<<"dim="<<l<<" run out of memory, imposing max="<<max;
		if(nthreads> 1)
			rec<<" (decrease number of threads)";
		rec<<endrec;
	}
	avmem_mb= mem/1024/1024;
	tmpmem_mb= (nthreads*per_size)/1024/1024;
	return max;
}


/*:6*/
;

/*:1*/
