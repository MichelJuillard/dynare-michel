// Contains the test for homotopy.
// The values in initval block are obtained from ramst_initial.mod 
// Result of the computation should be the same than in ramst_final.mod

@include "common.mod"

initval;
k = 12.75;
c = 1.5;
x = 1;
end;

homotopy_setup;
bet, 0.05, 0.1;
x, 2;
end;

steady(homotopy_mode = 2, homotopy_steps = 50);
//steady(homotopy_mode = 2, homotopy_steps = 50);
//steady(homotopy_mode = 3, homotopy_steps = 50);
