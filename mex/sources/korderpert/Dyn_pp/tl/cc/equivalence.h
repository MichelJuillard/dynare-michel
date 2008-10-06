/*1:*/
#line 46 "./equivalence.hweb"

#ifndef EQUIVALENCE_H
#define EQUIVALENCE_H

#include "int_sequence.h"

#include <vector> 
#include <list> 

using namespace std;

/*2:*/
#line 72 "./equivalence.hweb"

class OrdSequence{
vector<int> data;
public:
OrdSequence():data(){}
OrdSequence(const OrdSequence&s):data(s.data){}
const OrdSequence&operator= (const OrdSequence&s)
{data= s.data;return*this;}
bool operator==(const OrdSequence&s)const;
int operator[](int i)const;
bool operator<(const OrdSequence&s)const;
const vector<int> &getData()const
{return data;}
int length()const{return data.size();}
void add(int i);
void add(const OrdSequence&s);
bool has(int i)const;
void print(const char*prefix)const;
private:
double average()const;
};


/*:2*/
#line 57 "./equivalence.hweb"
;
/*3:*/
#line 101 "./equivalence.hweb"

class Permutation;
class Equivalence{
private:
int n;
list<OrdSequence> classes;
public:
typedef list<OrdSequence> ::const_iterator const_seqit;
typedef list<OrdSequence> ::iterator seqit;

/*6:*/
#line 178 "./equivalence.hweb"

Equivalence(int num);
Equivalence(int num,const char*dummy);
Equivalence(const Equivalence&e);
Equivalence(const Equivalence&e,int i1,int i2);

/*:6*/
#line 111 "./equivalence.hweb"
;
const Equivalence&operator= (const Equivalence&e);
bool operator==(const Equivalence&e)const;
bool operator!=(const Equivalence&e)const
{return!operator==(e);}
int getN()const{return n;}
int numClasses()const{return classes.size();}
void trace(IntSequence&out,int n)const;
void trace(IntSequence&out)const
{trace(out,numClasses());}
void trace(IntSequence&out,const Permutation&per)const;
void print(const char*prefix)const;
/*7:*/
#line 185 "./equivalence.hweb"

seqit begin(){return classes.begin();}
const_seqit begin()const{return classes.begin();}
seqit end(){return classes.end();}
const_seqit end()const{return classes.end();}

/*:7*/
#line 123 "./equivalence.hweb"
;
const_seqit find(int i)const;
seqit find(int i);
protected:
/*8:*/
#line 198 "./equivalence.hweb"

const_seqit findHaving(int i)const;
seqit findHaving(int i);
void insert(const OrdSequence&s);

/*:8*/
#line 127 "./equivalence.hweb"
;
};

/*:3*/
#line 58 "./equivalence.hweb"
;
/*4:*/
#line 137 "./equivalence.hweb"

class EquivalenceSet{
int n;
list<Equivalence> equis;
public:
typedef list<Equivalence> ::const_iterator const_iterator;
EquivalenceSet(int num);
void print(const char*prefix)const;
const_iterator begin()const
{return equis.begin();}
const_iterator end()const
{return equis.end();}
private:
bool has(const Equivalence&e)const;
void addParents(const Equivalence&e,list<Equivalence> &added);
};

/*:4*/
#line 59 "./equivalence.hweb"
;
/*5:*/
#line 161 "./equivalence.hweb"

class EquivalenceBundle{
vector<EquivalenceSet*> bundle;
public:
EquivalenceBundle(int nmax);
~EquivalenceBundle();
const EquivalenceSet&get(int n)const;
void generateUpTo(int nmax);
};

/*:5*/
#line 60 "./equivalence.hweb"
;

#endif


/*:1*/
