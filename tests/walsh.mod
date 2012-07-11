@#define Blocks = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15","16", "17", "18", "19", "20", "21", "22", "23", "24", "25","26", "27", "28", "29", "30", "31", "32", "33", "34", "35","36", "37", "38", "39", "40", "41", "42", "43", "44", "45","46", "47", "48", "49", "50", "51", "52", "53", "54", "55","56", "57", "58", "59", "60", "61", "62", "63", "64", "65","66", "67", "68", "69", "70", "71", "72", "73", "74", "75","76", "77", "78", "79", "80", "81", "82", "83", "84", "85","86", "87", "88", "89", "90", "91", "92", "93", "94", "95","96", "97", "98", "99", "100", "101", "102", "103", "104", "105", "106", "107", "108", "109", "110", "111", "112", "113", "114", "115","116", "117", "118", "119", "120", "121", "122", "123", "124", "125","126", "127", "128", "129", "130", "131", "132", "133", "134", "135","136", "137", "138", "139", "140", "141", "142", "143", "144", "145","146", "147", "148", "149", "150", "151", "152", "153", "154", "155","156", "157", "158", "159", "160", "161", "162", "163", "164", "165","166", "167", "168", "169", "170", "171", "172", "173", "174", "175","176", "177", "178", "179", "180", "181", "182", "183", "184", "185","186", "187", "188", "189", "190", "191", "192", "193", "194", "195","196", "197", "198", "199", "200"]

@#define blocks = Blocks[1:20]

@#for block in blocks
  var y_@{block} c_@{block} k_@{block} m_@{block} n_@{block} R_@{block} pi_@{block} z_@{block} u_@{block};
@#endfor

@#for block in blocks
varexo	e_@{block} sigma_@{block};
@#endfor

// sigma stands for phi in the eq 2.37 p.69

parameters alpha beta delta gamm phi1 eta a b rho  phi2 Psi thetass;
//phi1 stands for capital phi in eq.2.68 and 2.69
//phi2 stands for lowercase phi in eq. 2.66

alpha = 0.36;
beta = 0.989;
gamm = 0.5;
delta = 0.019;
phi1 = 2;
phi2 = 0;
eta = 1;
a = 0.95;
b = 2.56;
rho = 0.95;
Psi = 1.47630583;
thetass = 1.0125;

model(use_dll);

@#for block in blocks
(a*exp(c_@{block})^(1-b)+(1-a)*exp(m_@{block})^(1-b))^((b-phi1)/(1-b))*a*exp(c_@{block})^(-b) = (a*exp(c_@{block})^(1-b)+(1-a)*exp(m_@{block})^(1-b))^((b-phi1)/(1-b))*(1-a)*exp(m_@{block})^(-b)+beta*(a*exp(c_@{block}(+1))^(1-b)+(1-a)*exp(m_@{block}(+1))^(1-b))^((b-phi1)/(1-b))*a*exp(c_@{block}(+1))^(-b)/(1+pi_@{block}(+1));

Psi*(1-exp(n_@{block}))^(-eta)/(a*exp(c_@{block})^(-b)*(a*exp(c_@{block})^(1-b) + (1-a)*exp(m_@{block})^(1-b))^((b-phi1)/(1-b))) = (1-alpha)*exp(y_@{block})/exp(n_@{block});

(a*exp(c_@{block})^(1-b)+(1-a)*exp(m_@{block})^(1-b))^((b-phi1)/(1-b))*a*exp(c_@{block})^(-b) = beta*exp(R_@{block}(+1))*(a*exp(c_@{block}(+1))^(1-b)+(1-a)*exp(m_@{block}(+1))^(1-b))^((b-phi1)/(1-b))*a*exp(c_@{block}(+1))^(-b);

exp(R_@{block}) = alpha*exp(y_@{block})/exp(k_@{block}(-1)) + 1-delta;

exp(k_@{block}) = (1-delta)*exp(k_@{block}(-1))+exp(y_@{block})-exp(c_@{block});

exp(y_@{block}) = exp(z_@{block})*exp(k_@{block}(-1))^alpha*exp(n_@{block})^(1-alpha);

exp(m_@{block}) = exp(m_@{block}(-1))*(u_@{block}+thetass)/(1+pi_@{block});

z_@{block} = rho*z_@{block}(-1) + e_@{block};

u_@{block} = gamm*u_@{block}(-1) + phi2*z_@{block}(-1) + sigma_@{block};

@#endfor

end;

shocks;
@#for block in blocks
  var e_@{block}; stderr 0.007;
  var sigma_@{block};stderr 0.0089;
@#endfor
end;

steady_state_model;
// solving in levels
// calibrating n = 1/3 and recovering the value of Psi
// adapting solution Walsh (2003) p. 84
en = 1/3;
eR = 1/beta;
y_k = (1/alpha)*(1/beta-1+delta);
ek = en*y_k^(-1/(1-alpha));
ec = ek*(y_k-delta);
em = ec*(a/(1-a))^(-1/b)*((thetass-beta)/thetass)^(-1/b);
ey = ek*y_k;
Xss = a*ec^(1-b)*(1+(a/(1-a))^(-1/b)*((thetass-beta)/thetass)^((b-1)/b));
Psi = (1-alpha)*(ey/en)*Xss^((b-phi1)/(1-b))*a*ec^(-b)*(1-en)^eta;
@#for block in blocks
  pi_@{block} = thetass-1;
  n_@{block} = log(en);
  k_@{block} = log(ek);
  m_@{block} = log(em);
  c_@{block} = log(ec);
  y_@{block} = log(ey);
  R_@{block} = log(eR);
@#endfor
end;

steady;

loop_size = 100;

tic;
for i=1:loop_size
  stoch_simul(order=1,irf=0,noprint,nomoments,dr=default);
  dr = oo_.dr;
end;
e1 = toc;

tic;
for i=1:loop_size
  stoch_simul(order=1,irf=0,noprint,nomoments,dr=cycle_reduction);
  DR = oo_.dr;
end;
e2 = toc;

tic;
for i=1:loop_size
  stoch_simul(order=1,irf=0,noprint,nomoments,dr=cycle_reduction,dr_cycle_reduction_tol=1e-6);
  DR_ = oo_.dr;
end;
e3 = toc;

tic;
for i=1:loop_size
  stoch_simul(order=1,irf=0,noprint,nomoments,dr=cycle_reduction,dr_cycle_reduction_tol=1e-5);
  DR__ = oo_.dr;
end;
e4 = toc;

tic;
for i=1:loop_size
  stoch_simul(order=1,irf=0,noprint,nomoments,dr=cycle_reduction,dr_cycle_reduction_tol=1e-4);
  DR___ = oo_.dr;
end;
e5 = toc;

tic;
for i=1:loop_size
  stoch_simul(order=1,irf=0,noprint,nomoments,dr=cycle_reduction,dr_cycle_reduction_tol=1e-3);
  DR____ = oo_.dr;
end;
e6 = toc;

tic;
for i=1:loop_size
  stoch_simul(order=1,irf=0,noprint,nomoments,dr=cycle_reduction,dr_cycle_reduction_tol=1e-2);
  DR_____ = oo_.dr;
end;
e7 = toc;

tic;
for i=1:loop_size
  stoch_simul(order=1,irf=0,noprint,nomoments,dr=cycle_reduction,dr_cycle_reduction_tol=200);
  DR______ = oo_.dr;
end;
e8 = toc;

disp(['Elapsed time (QZ algorithm): ' num2str(e1)])
disp(['Elapsed time (Cycle reduction algorithm, tol = 1e-7): ' num2str(e2)])
disp(['Elapsed time (Cycle reduction algorithm, tol = 1e-6): ' num2str(e3)])
disp(['Elapsed time (Cycle reduction algorithm, tol = 1e-5): ' num2str(e4)])
disp(['Elapsed time (Cycle reduction algorithm, tol = 1e-4): ' num2str(e5)])
disp(['Elapsed time (Cycle reduction algorithm, tol = 1e-3): ' num2str(e6)])
disp(['Elapsed time (Cycle reduction algorithm, tol = 1e-2): ' num2str(e7)])
disp(['Elapsed time (Cycle reduction algorithm, tol = 1e-1): ' num2str(e8)])

if e1<e2
   disp(['cycle reduction algorithm is slower than QZ algorithm in solving the model (' int2str(M_.endo_nbr) ' equations)!'])
end

dx = max(max(dr.ghx-DR.ghx));
du = max(max(dr.ghu-DR.ghu));
dx_ = max(max(dr.ghx-DR_.ghx));
du_ = max(max(dr.ghu-DR_.ghu));
dx__ = max(max(dr.ghx-DR__.ghx));
du__ = max(max(dr.ghu-DR__.ghu));
dx___ = max(max(dr.ghx-DR___.ghx));
du___ = max(max(dr.ghu-DR___.ghu));
dx____ = max(max(dr.ghx-DR____.ghx));
du____ = max(max(dr.ghu-DR____.ghu));
dx_____ = max(max(dr.ghx-DR_____.ghx));
du_____ = max(max(dr.ghu-DR_____.ghu));
dx______ = max(max(dr.ghx-DR______.ghx));
du______ = max(max(dr.ghu-DR______.ghu));


if dx>1e-8
   error('QZ and Cycle reduction algorithms don''t return the same results (tolerance parameter is 1e-7)!')
end

if du>1e-8
   error('QZ and Cycle reduction algorithms don''t return the same results (tolerance parameter is 1e-7)!')
end

if dx_>1e-8 || du_>1e-8
   disp('QZ and Cycle reduction algorithms don''t return the same results (tolerance parameter is 1e-6)!')
else
   [dx, dx_]
   [du, du_]
end

if dx__>1e-8 || du__>1e-8
   disp('QZ and Cycle reduction algorithms don''t return the same results (tolerance parameter is 1e-5)!')
else
   [dx, dx_, dx__]
   [du, du_, du__]
end

if dx___>1e-8 || du___>1e-8
   disp('QZ and Cycle reduction algorithms don''t return the same results (tolerance parameter is 1e-4)!')
else
   [dx, dx_, dx__, dx___]
   [du, du_, du__, du___]
end

if dx____>1e-8 || du____>1e-8
   disp('QZ and Cycle reduction algorithms don''t return the same results (tolerance parameter is 1e-3)!')
else
   [dx, dx_, dx__, dx___,dx____]
   [du, du_, du__, du___,du____]
end

if dx_____>1e-8 || du_____>1e-8
   disp('QZ and Cycle reduction algorithms don''t return the same results (tolerance parameter is 1e-2)!')
else
   [dx, dx_, dx__, dx___,dx____,dx_____]
   [du, du_, du__, du___,du____,du_____]
end

if dx_____>1e-8 || du_____>1e-8
   disp('QZ and Cycle reduction algorithms don''t return the same results (tolerance parameter is 1e-1)!')
else
   [dx, dx_, dx__, dx___,dx____,dx_____,dx______]
   [du, du_, du__, du___,du____,du_____,du______]
end