/*1:*/


#include "random.h"

#include <stdlib.h> 
#include <limits> 
#include <cmath> 

/*2:*/

int RandomGenerator::int_uniform()
{
	double s= std::numeric_limits<int> ::max()*uniform();
	return(int)s;
}

/*:2*/
;
/*3:*/

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
;
SystemRandomGenerator system_random_generator;
/*4:*/

double SystemRandomGenerator::uniform()
{
#ifndef __MINGW32__
	return drand48();
#else
	return((double)rand())/RAND_MAX;
#endif
}

/*:4*/
;
/*5:*/

void SystemRandomGenerator::initSeed(int seed)
{
#ifndef __MINGW32__
	srand48(seed);
#else
	srand(seed);
#endif
}

/*:5*/
;

/*:1*/
