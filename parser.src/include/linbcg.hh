#ifndef LINBCG_HH_INCLUDED
#define LINBCG_HH_INCLUDED

#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include "mex.h"
using namespace std;

const int debile=10;

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
	inline T & operator[](const int i);	//i'th element
	inline const T & operator[](const int i) const;
	inline int size() const;
	~NRVec();
};

template <class T>
NRVec<T>::NRVec() : nn(0), v(0) {}

template <class T>
NRVec<T>::NRVec(int n) : nn(n), v(new T[n]) {}

template <class T>
NRVec<T>::NRVec(const T& a, int n) : nn(n), v(new T[n])
{
	for(int i=0; i<n; i++)
		v[i] = a;
}

template <class T>
NRVec<T>::NRVec(const T *a, int n) : nn(n), v(new T[n])
{
	for(int i=0; i<n; i++)
		v[i] = *a++;
}

template <class T>
NRVec<T>::NRVec(const NRVec<T> &rhs) : nn(rhs.nn), v(new T[nn])
{
	for(int i=0; i<nn; i++)
		v[i] = rhs[i];
}

template <class T>
NRVec<T> & NRVec<T>::operator=(const NRVec<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
	if (this != &rhs)
	{
		if (nn != rhs.nn) {
			if (v != 0) delete [] (v);
			nn=rhs.nn;
			v= new T[nn];
		}
		for (int i=0; i<nn; i++)
			v[i]=rhs[i];
	}
	return *this;
}

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
	return v[i];
}

template <class T>
inline const T & NRVec<T>::operator[](const int i) const	//subscripting
{
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
		delete[] (v);
}

typedef const NRVec<int> Vec_I_INT;
typedef NRVec<int> Vec_INT, Vec_O_INT, Vec_IO_INT;
typedef const NRVec<DP> Vec_I_DP;
typedef NRVec<DP> Vec_DP, Vec_O_DP, Vec_IO_DP;

class LinBCG
{
  Vec_INT *ija_p;
  Vec_DP *sa_p;
  public:
  LinBCG();
  void Init(int periods, int y_kmin, int y_kmax, int Size, std::map<std::pair<std::pair<int, int> ,int>, int> IM, int* index_vara, int* index_equa);
  void asolve(Vec_I_DP &b, Vec_O_DP &x, const int itrnsp);
  void atimes(Vec_I_DP &x, Vec_O_DP &r, const int itrnsp);
  DP snrm(Vec_I_DP &sx, const int itol);
  void sprsax(Vec_I_DP &sa, Vec_I_INT &ija, Vec_I_DP &x, Vec_O_DP &b);
  void sprstx(Vec_I_DP &sa, Vec_I_INT &ija, Vec_I_DP &x, Vec_O_DP &b);
  void linbcg(Vec_I_DP &b, Vec_IO_DP &x, const int itol, const DP tol,
              const int itmax, int &iter, DP &err);
};
#endif
