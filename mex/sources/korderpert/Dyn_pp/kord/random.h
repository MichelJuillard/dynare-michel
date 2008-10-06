/*1:*/
#line 9 "./random.hweb"

#ifndef RANDOM_H
#define RANDOM_H

/*2:*/
#line 22 "./random.hweb"

class RandomGenerator{
public:
virtual double uniform()= 0;
int int_uniform();
double normal();
};

/*:2*/
#line 13 "./random.hweb"
;
/*3:*/
#line 32 "./random.hweb"

class SystemRandomGenerator:public RandomGenerator{
public:
double uniform();
void initSeed(int seed);
};

/*:3*/
#line 14 "./random.hweb"
;
extern SystemRandomGenerator system_random_generator;

#endif

/*:1*/
