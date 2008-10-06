/*1:*/
#line 28 "./int_sequence.hweb"

#ifndef INT_SEQUENCE_H
#define INT_SEQUENCE_H


#include <string.h> 
#include <vector> 

using namespace std;

/*2:*/
#line 46 "./int_sequence.hweb"

class Symmetry;
class IntSequence{
int*data;
int length;
bool destroy;
public:
/*3:*/
#line 83 "./int_sequence.hweb"

IntSequence(int l)
:data(new int[l]),length(l),destroy(true){}
IntSequence(int l,int n)
:data(new int[l]),length(l),destroy(true)
{for(int i= 0;i<length;i++)data[i]= n;}
IntSequence(const IntSequence&s)
:data(new int[s.length]),length(s.length),destroy(true)
{memcpy(data,s.data,length*sizeof(int));}
IntSequence(IntSequence&s,int i1,int i2)
:data(s.data+i1),length(i2-i1),destroy(false){}
IntSequence(const IntSequence&s,int i1,int i2)
:data(new int[i2-i1]),length(i2-i1),destroy(true)
{memcpy(data,s.data+i1,sizeof(int)*length);}
IntSequence(const Symmetry&sy,const vector<int> &se);
IntSequence(const Symmetry&sy,const IntSequence&se);
IntSequence(int i,const IntSequence&s);
IntSequence(int i,const IntSequence&s,int pos);
IntSequence(int l,const int*d)
:data(new int[l]),length(l),destroy(true)
{memcpy(data,d,sizeof(int)*length);}


/*:3*/
#line 53 "./int_sequence.hweb"
;
/*4:*/
#line 107 "./int_sequence.hweb"

const IntSequence&operator= (const IntSequence&s);
virtual~IntSequence()
{if(destroy)delete[]data;}
bool operator==(const IntSequence&s)const;
bool operator!=(const IntSequence&s)const
{return!operator==(s);}
int&operator[](int i)
{return data[i];}
int operator[](int i)const
{return data[i];}
int size()const
{return length;}

/*:4*/
#line 54 "./int_sequence.hweb"
;
/*5:*/
#line 124 "./int_sequence.hweb"

bool operator<(const IntSequence&s)const;
bool operator<=(const IntSequence&s)const
{return(operator==(s)||operator<(s));}
bool lessEq(const IntSequence&s)const;
bool less(const IntSequence&s)const;


/*:5*/
#line 55 "./int_sequence.hweb"
;
void sort();
void monotone();
void pmonotone(const Symmetry&s);
int sum()const;
int mult(int i1,int i2)const;
int mult()const
{return mult(0,length);}
void add(int i);
void add(int f,const IntSequence&s);
int getPrefixLength()const;
int getNumDistinct()const;
int getMax()const;
bool isPositive()const;
bool isConstant()const;
bool isSorted()const;
void print()const;
};

/*:2*/
#line 38 "./int_sequence.hweb"
;

#endif

/*:1*/
