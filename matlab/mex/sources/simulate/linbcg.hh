/*
 * Copyright (C) 2007-2008 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LINBCG_HH_INCLUDED
#define LINBCG_HH_INCLUDED

#include <iostream>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <complex>
#include <map>
#include <string>
#include <time.h>
#include "mex.h"

using namespace std;

const int debile=10;
const double tol=1e-10;

typedef double DP;


template <class T>
class NRVec {
private:
	int nn;	// size of array. upper index is nn-1
	T *v;
public:
	NRVec();
	explicit NRVec(int n);		// Zero-based array
	NRVec(const T &a, int n);	//initialize to constant value
	NRVec(const T *a, int n);	// Initialize to array
	NRVec(const NRVec &rhs);	// Copy constructor
	NRVec & operator=(const NRVec &rhs);	//assignment
	NRVec & operator=(const T &a);	//assign a to every element
  //	NRVec & operator=(const T* rhs);
	inline T & operator[](const int i);	//i'th element
	inline const T & operator[](const int i) const;
	inline int size() const;
	inline void Format(T* rhs, int n);
	inline void Format(const NRVec rhs);
	inline void Format(int n);
	inline void Print();
	~NRVec();
};

template <class T>
NRVec<T>::NRVec() : nn(0), v(0) {}

template <class T>
NRVec<T>::NRVec(int n) : nn(n), v(/*new T[n]*/(T*)mxMalloc(sizeof(T)*n)) {}

template <class T>
NRVec<T>::NRVec(const T& a, int n) : nn(n), v(/*new T[n]*/(T*)mxMalloc(sizeof(T)*n))
{
	for(int i=0; i<n; i++)
		v[i] = a;
}

template <class T>
NRVec<T>::NRVec(const T *a, int n) : nn(n), v(/*new T[n]*/(T*)mxMalloc(sizeof(T)*n))
{
	for(int i=0; i<n; i++)
		v[i] = *a++;
}

template <class T>
NRVec<T>::NRVec(const NRVec<T> &rhs) : nn(rhs.nn), v(/*new T[nn]*/(T*)mxMalloc(sizeof(T)*rhs.nn))
{
	for(int i=0; i<nn; i++)
		v[i] = rhs[i];
}

template <class T>
void NRVec<T>::Format(T* rhs, int n)
{
   if(v)
     mxFree(v);
   v=NULL;
   nn=n;
   v=(T*)mxMalloc(sizeof(T)*n);
   for(int i=0; i<nn; i++)
		 v[i] = rhs[i];
}

template <class T>
void NRVec<T>::Format(const NRVec rhs)
{
   if(v)
     mxFree(v);
   v=NULL;
   nn=rhs.nn;
   v=(T*)mxMalloc(sizeof(T)*nn);
   for(int i=0; i<nn; i++)
		 v[i] = rhs[i];
}


template <class T>
void NRVec<T>::Format(int n)
{
   if(v)
     mxFree(v);
   v=NULL;
   nn=n;
   v=(T*)mxMalloc(sizeof(T)*n);
}


template <class T>
NRVec<T> & NRVec<T>::operator=(const NRVec<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
	if (this != &rhs)
	  {
		  if (nn != rhs.nn)
		    {
			    if (v != 0)
			      mxFree(v);
			    nn=rhs.nn;
			    v= (T*)mxMalloc(sizeof(T)*nn);
		    }
		  for (int i=0; i<nn; i++)
			  v[i]=rhs[i];
	  }
	return *this;
}

template <class T>
void NRVec<T>::Print()
{
  int i;
  for(i=0;i<nn;i++)
    mexPrintf("%f\n",v[i]);
}

/*template <class T>
NRVec<T> & NRVec<T>::operator=(const T* rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
	if (this != &rhs)
	  {
       if (v != 0)
			    mxFree(v);
       v= (T*)mxMalloc(sizeof(T)*nn);
	  }
  for (int i=0; i<nn; i++)
	  v[i]=rhs[i];
	return *this;
}*/

template <class T>
NRVec<T> & NRVec<T>::operator=(const T &a)	//assign a to every element
{
	for (int i=0; i<nn; i++)
		v[i]=a;
	return *this;
}

template <class T>
inline T & NRVec<T>::operator[](const int i)	//subscripting
{
  if(i>=nn || i<0)
    mexPrintf("Error: out of range in Vect[] operator\n");
	return v[i];
}

template <class T>
inline const T & NRVec<T>::operator[](const int i) const	//subscripting
{
  if(i>=nn || i<0)
    mexPrintf("Error: out of range in Vect[] operator\n");
	return v[i];
}

template <class T>
inline int NRVec<T>::size() const
{
	return nn;
}

extern double  *u;
template <class T>
NRVec<T>::~NRVec()
{
	if (v != 0)
		mxFree(v);
}


template <class T>
class NRMat {
private:
	int nn;
	int mm;
	T **v;
public:
	NRMat();
	NRMat(int n, int m);			// Zero-based array
	NRMat(const T &a, int n, int m);	//Initialize to constant
	NRMat(const T *a, int n, int m);	// Initialize to array
	NRMat(const NRMat &rhs);		// Copy constructor
	NRMat & operator=(const NRMat &rhs);	//assignment
	NRMat & operator=(const T &a);		//assign a to every element
	inline T* operator[](const int i);	//subscripting: pointer to row i
	inline const T* operator[](const int i) const;
	inline int nrows() const;
	inline int ncols() const;
	inline void Format(const int n,const int m);
	void Print();
	~NRMat();
};

template <class T>
NRMat<T>::NRMat() : nn(0), mm(0), v(0) {}

template <class T>
NRMat<T>::NRMat(int n, int m) :
                               nn(n),
                               mm(m),
                               v((T**)mxMalloc(sizeof(T*)*n))
{
	v[0] = /*new T[m*n]*/(T*)mxMalloc(sizeof(T)*m*n);
	for (int i=1; i< n; i++)
		v[i] = v[i-1] + m;
}

template <class T>
NRMat<T>::NRMat(const T &a, int n, int m) : nn(n), mm(m), v(/*new T*[n]*/(T**)mxMalloc(sizeof(T*)*n))
{
	int i,j;
	v[0] = /*new T[m*n]*/(T*)mxMalloc(sizeof(T)*m*n);
	for (i=1; i< n; i++)
		v[i] = v[i-1] + m;
	for (i=0; i< n; i++)
		for (j=0; j<m; j++)
			v[i][j] = a;
}

template <class T>
NRMat<T>::NRMat(const T *a, int n, int m) : nn(n), mm(m), v(/*new T*[n]*/(T**)mxMalloc(sizeof(T*)*n))
{
	int i,j;
	v[0] = /*new T[m*n]*/(T*)mxMalloc(sizeof(T)*m*n);
	for (i=1; i< n; i++)
		v[i] = v[i-1] + m;
	for (i=0; i< n; i++)
		for (j=0; j<m; j++)
			v[i][j] = *a++;
}

template <class T>
NRMat<T>::NRMat(const NRMat &rhs) : nn(rhs.nn), mm(rhs.mm), v(/*new T*[nn]*/(T**)mxMalloc(sizeof(T*)*nn))
{
	int i,j;
	v[0] = /*new T[mm*nn]*/(T*)mxMalloc(sizeof(T)*mm*nn);
	for (i=1; i< nn; i++)
		v[i] = v[i-1] + mm;
	for (i=0; i< nn; i++)
		for (j=0; j<mm; j++)
			v[i][j] = rhs[i][j];
}

template <class T>
NRMat<T> & NRMat<T>::operator=(const NRMat<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
	if (this != &rhs) {
		int i,j;
		if (nn != rhs.nn || mm != rhs.mm) {
			if (v != 0)
			  {
			    mxFree(v[0]);
			    mxFree(v);
				  /*delete[] (v[0]);
				  delete[] (v);*/
			  }
			nn=rhs.nn;
			mm=rhs.mm;
			//v = new T*[nn];
			v = (T**)mxMalloc(sizeof(T*)*nn);
			v[0] = /*new T[mm*nn]*/(T*)mxMalloc(sizeof(T)*mm*nn);
		}
		for (i=1; i< nn; i++)
			v[i] = v[i-1] + mm;
		for (i=0; i< nn; i++)
			for (j=0; j<mm; j++)
				v[i][j] = rhs[i][j];
	}
	return *this;
}

template <class T>
NRMat<T> & NRMat<T>::operator=(const T &a)	//assign a to every element
{
	for (int i=0; i< nn; i++)
		for (int j=0; j<mm; j++)
			v[i][j] = a;
	return *this;
}

template <class T>
inline T* NRMat<T>::operator[](const int i)	//subscripting: pointer to row i
{
  if(i>=nn || i<0)
    mexPrintf("Error: out of range in Mat[] operator\n");
	return v[i];
}

template <class T>
inline const T* NRMat<T>::operator[](const int i) const
{
  if(i>=nn || i<0)
    mexPrintf("Error: out of range in Mat[] operator\n");
	return v[i];
}

template <class T>
inline int NRMat<T>::nrows() const
{
	return nn;
}

template <class T>
inline int NRMat<T>::ncols() const
{
	return mm;
}

template <class T>
inline void NRMat<T>::Format(const int n,const int m)
{
  if (v != 0)
    {
		  /*delete[] (v[0]);
		  delete[] (v);*/
		  mxFree(v[0]);
      mxFree(v);
	  }
	nn=n;
	mm=m;
	v = (T**)mxMalloc(sizeof(T*)*n);//new T*[n];
	v[0] = (T*)mxMalloc(sizeof(T)*m*n);//new T[m*n];
	for (int i=1; i< n; i++)
		v[i] = v[i-1] + m;
}



template <class T>
NRMat<T>::~NRMat()
{
	if (v != 0)
	  {
	    mxFree(v[0]);
      mxFree(v);
		  /*delete[] (v[0]);
		  delete[] (v);*/
	  }
}


template <class T>
void NRMat<T>::Print()
{
  int i,j;
  mexPrintf("nn=%d, mm=%d\n",nn,mm);
  for(i=0;i<nn;i++)
    {
      for(j=0;j<mm;j++)
        mexPrintf("%f ",v[i][j]);
      mexPrintf("\n");
    }
}


typedef const NRVec<int> Vec_I_INT;
typedef NRVec<int> Vec_INT, Vec_O_INT, Vec_IO_INT;
typedef const NRVec<DP> Vec_I_DP;
typedef NRVec<DP> Vec_DP, Vec_O_DP, Vec_IO_DP;

typedef const NRMat<DP> Mat_I_DP;
typedef NRMat<DP> Mat_DP, Mat_O_DP, Mat_IO_DP;




class LinBCG
{
  Vec_INT ija_p, ijat_p;
  Vec_DP sa_p, sat_p;
  private:


  clock_t time00;
  string filename;
  double res1, res2, slowc, slowc_save, *ya, *direction, max_res;
  int iter;

  void sprs_swap_line_copy(map<std::pair<int, int>, double> &list, int &pos_2_blck, int begin, int end);
  void sprs_swap_line_exchange(map<std::pair<int, int>, double> &list, int &pos_2_blck, int LS, int LD);

  public:
  LinBCG();
  void SolveLinear(int periods, int y_kmin, int y_kmax, int Size, std::map<std::pair<std::pair<int, int> ,int>, int> IM, int* index_vara, int* index_equa, int y_size, double* y, bool print_it, bool cvg, Mat_DP &a, Vec_IO_INT &indx);
  void Preconditioner(int periods, int y_kmin, int y_kmax, int Size, std::map<std::pair<std::pair<int, int> ,int>, int> IM, int* index_vara, int* index_equa, int y_size, double* y, bool print_it, int type, Mat_O_DP &a, Vec_O_INT &indx);
  void asolve(Vec_I_DP &b, Vec_O_DP &x, const int itrnsp, Mat_DP &a, Vec_IO_INT &indx, const int periods);
  void atimes(Vec_I_DP &x, Vec_O_DP &r, const int itrnsp);
  DP snrm(Vec_I_DP &sx, const int itol);
  void sprsax(Vec_I_DP &sa, Vec_I_INT &ija, Vec_I_DP &x, Vec_O_DP &b);
  void sprstx(Vec_I_DP &sa, Vec_I_INT &ija, Vec_I_DP &x, Vec_O_DP &b);
  void sprsin(double *a, int NbRow, int NbCol, const DP thresh);
  void sprsin(Vec_DP s, Vec_INT ij);
  void sprsin(DP *s, int *ij, int size);
  void sprsin(Mat_DP &a, const DP thresh);
  void sprsout(Mat_DP &a);
  void BiCG(Vec_I_DP &b, Vec_IO_DP &x, const int itol, const DP tol,
              const int itmax, int &iter, DP &err, Mat_DP &a, Vec_IO_INT &indx, const int periods);
  void BiCGStab(Vec_I_DP &b, Vec_IO_DP &x, const int itol, const DP tol,
	            const int itmax, int &iter, DP &err, Mat_DP &a, Vec_IO_INT &indx, const int periods);
  void sprsprt();
  void sprs_sprt();
  void sprs_col_index();
  void sprs_swap_line(int L0, int L1);
  void sprsludcmp(Mat_DP &a, Vec_O_INT &indx, DP &d);
  void lubksb_blck(Mat_I_DP &a, Vec_I_INT &indx, Vec_IO_DP &b, const int periods);
  void lubksb_blck_trp(Mat_I_DP &a, Vec_I_INT &indx, Vec_IO_DP &b, const int periods);
  void lubksb(Mat_I_DP &a, Vec_I_INT &indx, Vec_IO_DP &b);
  void lubksb_trp(Mat_I_DP &a, Vec_I_INT &indx, Vec_IO_DP &b);
  void sprs_LU(Vec_O_DP &sL, Vec_O_INT &ijL, Vec_O_DP &sU, Vec_O_INT &ijU);
  void Initialize(string filename_arg, double res1_arg, double res2_arg, double max_res_arg, double slowc_arg, double *ya_arg, double *direction_arg, int iter_arg);
};


#endif
