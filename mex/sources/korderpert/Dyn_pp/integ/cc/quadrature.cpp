/*1:*/

#include "quadrature.h"
#include "precalc_quadrature.dat"

#include <cmath> 

/*2:*/

void OneDPrecalcQuadrature::calcOffsets()
{
offsets[0]= 0;
for(int i= 1;i<num_levels;i++)
offsets[i]= offsets[i-1]+num_points[i-1];
}

/*:2*/
;
/*3:*/

GaussHermite::GaussHermite()
:OneDPrecalcQuadrature(gh_num_levels,gh_num_points,gh_weights,gh_points){}

/*:3*/
;
/*4:*/

GaussLegendre::GaussLegendre()
:OneDPrecalcQuadrature(gl_num_levels,gl_num_points,gl_weights,gl_points){}

/*:4*/
;
/*5:*/

double NormalICDF::get(double x)
{
double xx= (2*normal_icdf_end-1)*std::abs(x-0.5);
int i= (int)floor(xx/normal_icdf_step);
double xx1= normal_icdf_step*i;
double yy1= normal_icdf_data[i];
double y;
if(i<normal_icdf_num-1){
double yy2= normal_icdf_data[i+1];
y= yy1+(yy2-yy1)*(xx-xx1)/normal_icdf_step;
}else{
y= yy1;
}
if(x> 0.5)
return y;
else
return-y;
}


/*:5*/
;

/*:1*/
