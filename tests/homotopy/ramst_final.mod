// Computes the steady state which should be arrived at in ramst_homotopy.mod

@#include "common.mod"

bet = 0.1;

initval;
x = 2;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
end;

steady;
