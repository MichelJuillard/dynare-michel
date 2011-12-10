// comments.mod File to ensure that C/C++/Matlab-style comments are filtered out by
// Dynare and that native Matlab statements (both single and multi-line) pass through to comments.m
var y, c, k, a, h, b;
varexo e, u;

parameters beta, rho, alpha, delta, theta, psi, tau;

alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;

phi   = 0.1;

dsge_prior_weight = .8;
aaa = 8343553 / 3737 * 90
bbb = 8343553 / 3737 * 39            %comment
ccc = [ 3 1 4; 49 23 2 ; 0 90 0] ;           // comment
ddd = [ 4 2 3; rho tau beta; 3 1 43]      /*comment*/
eee = {alpha}/*the comment across several lines ****
*****/
fff = [ 1 delta 4 ; ...       /* COMMENT
                                     keeps on
                                     going*/
     1     0 4 ; ...     // comment */
         phi 9 4 ]        % comment

disp(' %% This is not a comment %% ')

model;
c*theta*h^(1+psi)=(1-alpha)*y;
k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
    *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));
y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
k = exp(b)*(y-c)+(1-delta)*k(-1);
a = rho*a(-1)+tau*b(-1) + e;
b = tau*a(-1)+rho*b(-1) + u;
end;

initval;
y = 1.08068253095672;
c = 0.80359242014163;
h = 0.29175631001732;
k = 11.08360443260358;
a = 0;
b = 0;
e = 0;
u = 0;
end;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
var e, u = phi*0.009*0.009;
end;

stoch_simul;
