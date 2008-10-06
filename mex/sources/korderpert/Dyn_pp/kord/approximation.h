/*1:*/
#line 53 "./approximation.hweb"

#ifndef APPROXIMATION_H
#define APPROXIMATION_H

#include "dynamic_model.h"
#include "decision_rule.h"
#include "korder.h"
#include "journal.h"

/*2:*/
#line 80 "./approximation.hweb"

class ZAuxContainer:public StackContainer<FGSTensor> ,public FoldedStackContainer{
public:
typedef StackContainer<FGSTensor> ::_Ctype _Ctype;
typedef StackContainer<FGSTensor> ::itype itype;
ZAuxContainer(const _Ctype*gss,int ngss,int ng,int ny,int nu);
itype getType(int i,const Symmetry&s)const;
};



/*:2*/
#line 62 "./approximation.hweb"
;
/*3:*/
#line 115 "./approximation.hweb"

class Approximation{
DynamicModel&model;
Journal&journal;
FGSContainer*rule_ders;
FGSContainer*rule_ders_ss;
FoldDecisionRule*fdr;
UnfoldDecisionRule*udr;
const PartitionY ypart;
const FNormalMoments mom;
IntSequence nvs;
int steps;
TwoDMatrix ss;
public:
Approximation(DynamicModel&m,Journal&j,int ns);
virtual~Approximation();

const FoldDecisionRule&getFoldDecisionRule()const;
const UnfoldDecisionRule&getUnfoldDecisionRule()const;
const TwoDMatrix&getSS()const
{return ss;}
const DynamicModel&getModel()const
{return model;}

void walkStochSteady();
TwoDMatrix*calcYCov()const;
protected:
void approxAtSteady();
void calcStochShift(Vector&out,double at_sigma)const;
void saveRuleDerivs(const FGSContainer&g);
void check(double at_sigma)const;
};


/*:3*/
#line 63 "./approximation.hweb"
;

#endif


/*:1*/
