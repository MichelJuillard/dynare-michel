// Computes the steady state used in initval block of ramst_homotopy.mod

@#include "common.mod"

bet=0.05;

initval;
x = 1;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
end;

steady;
