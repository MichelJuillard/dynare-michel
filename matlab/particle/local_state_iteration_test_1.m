function [T,t,LOG] = local_state_iteration_test_1()
try
addpath ../matlab

n = 2;
q = 3;

yhat = zeros(n,1);
epsilon = zeros(q,1);
ghx = rand(n,n);
ghu = rand(n,q);
constant = ones(n,1);
ghxx = rand(n,n*n);
ghuu = rand(n,q*q);
ghxu = rand(n,n*q);
yhat_ = zeros(n,1);
ss = ones(n,1);

y1 = local_state_iteration(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu);
[y2,y2_] = local_state_iteration(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,yhat_,ss);

t(1) = dyn_assert(y1,ones(n,1));
t(2) = dyn_assert(y2,ones(n,1));
t(3) = dyn_assert(y2_,ones(n,1));
T = all(t);
LOG = NaN;
catch exception
LOG = getReport(exception,'extended');
T = NaN;
t = NaN;
end
