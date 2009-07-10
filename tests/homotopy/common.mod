var c k;
varexo x;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
aa=0.5;
bet=0.05;

model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
end;

initval;
x = 1;
k = ((delt+bet)/(aa*x*alph))^(1/(alph-1));
c = aa*x*k^alph-delt*k;
end;

homotopy_setup;
bet, 0.05, 0.1;
x, 2;
end;

steady(homotopy_mode = @{homotopy_mode}, homotopy_steps = 50);

if abs(oo_.steady_state(1)/(aa*oo_.exo_steady_state(1)*oo_.steady_state(2)^alph-delt*oo_.steady_state(2)) - 1) > 1e-4
   error('Error in homotopy for c')
end

if abs(oo_.steady_state(2)/((delt+get_param_by_name('bet'))/(aa*oo_.exo_steady_state(1)*alph))^(1/(alph-1)) - 1) > 1e-4
   error('Error in homotopy for k')
end
