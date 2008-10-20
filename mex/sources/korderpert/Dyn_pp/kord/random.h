/*1:*/

#ifndef RANDOM_H
#define RANDOM_H

/*2:*/

class RandomGenerator{
public:
	virtual double uniform()= 0;
	int int_uniform();
	double normal();
};

/*:2*/
;
/*3:*/

class SystemRandomGenerator:public RandomGenerator{
public:
	double uniform();
	void initSeed(int seed);
};

/*:3*/
;
extern SystemRandomGenerator system_random_generator;

#endif

/*:1*/
