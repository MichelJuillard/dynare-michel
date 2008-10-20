/*1:*/

#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "equivalence.h"
#include "int_sequence.h"

#include <list> 
#include <vector> 

/*2:*/

class Symmetry:public IntSequence{
public:
/*3:*/

Symmetry(int len,const char*dummy)
:IntSequence(len,0){}
Symmetry(int i1)
:IntSequence(1,i1){}
Symmetry(int i1,int i2)
:IntSequence(2){operator[](0)= i1;operator[](1)= i2;}
Symmetry(int i1,int i2,int i3)
:IntSequence(3)
{
operator[](0)= i1;
operator[](1)= i2;
operator[](2)= i3;
}
Symmetry(int i1,int i2,int i3,int i4)
:IntSequence(4)
{
operator[](0)= i1;
operator[](1)= i2;
operator[](2)= i3;
operator[](3)= i4;
}
Symmetry(const Symmetry&s)
:IntSequence(s){}
Symmetry(const Symmetry&s,const OrdSequence&cl)
:IntSequence(s,cl.getData()){}
Symmetry(Symmetry&s,int len)
:IntSequence(s,s.size()-len,s.size()){}
Symmetry(const IntSequence&s);

/*:3*/
;
int num()const
{return size();}
int dimen()const
{return sum();}
int findClass(int i)const;
bool isFull()const;
};

/*:2*/
;
/*4:*/

class SymmetrySet{
Symmetry run;
int dim;
public:
SymmetrySet(int d,int length)
:run(length,""),dim(d){}
SymmetrySet(SymmetrySet&s,int d)
:run(s.run,s.size()-1),dim(d){}
int dimen()const
{return dim;}
const Symmetry&sym()const
{return run;}
Symmetry&sym()
{return run;}
int size()const
{return run.size();}
};

/*:4*/
;
/*5:*/

class symiterator{
SymmetrySet&s;
symiterator*subit;
SymmetrySet*subs;
bool end_flag;
public:
symiterator(SymmetrySet&ss);
~symiterator();
symiterator&operator++();
bool isEnd()const
{return end_flag;}
const Symmetry&operator*()const
{return s.sym();}
};


/*:5*/
;
/*6:*/

class InducedSymmetries:public vector<Symmetry> {
public:
InducedSymmetries(const Equivalence&e,const Symmetry&s);
InducedSymmetries(const Equivalence&e,const Permutation&p,const Symmetry&s);
void print()const;
};



/*:6*/
;

#endif

/*:1*/
