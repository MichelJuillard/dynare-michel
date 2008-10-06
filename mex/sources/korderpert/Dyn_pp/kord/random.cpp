/*1:*/
#line 5 "./random.cweb"


#include "random.h"

#include <stdlib.h> 
#include <limits> 
#include <cmath> 

/*2:*/
#line 20 "./random.cweb"

int RandomGenerator::int_uniform()
{
double s= std::numeric_limits<int> ::max()*uniform();
return(int)s;
}

/*:2*/
#line 13 "./random.cweb"
;
/*3:*/
#line 28 "./random.cweb"

double RandomGenerator::normal()
{
double x1,x2;
double w;
do{
x1= 2*uniform()-1;
x2= 2*uniform()-1;
w= x1*x1+x2*x2;
}while(w>=1.0||w<1.0e-30);
return x1*std::sqrt((-2.0*std::log(w))/w);
}

/*:3*/
#line 14 "./random.cweb"
;
SystemRandomGenerator system_random_generator;
/*4:*/
#line 42 "./random.cweb"

double SystemRandomGenerator::uniform()
{
#ifndef __MINGW32__
return drand48();
#else
return((double)rand())/RAND_MAX;
#endif
}

/*:4*/
#line 16 "./random.cweb"
;
/*5:*/
#line 53 "./random.cweb"

void SystemRandomGenerator::initSeed(int seed)
{
#ifndef __MINGW32__
srand48(seed);
#else
srand(seed);
#endif
}

/*:5*/
#line 17 "./random.cweb"
;

/*:1*/
