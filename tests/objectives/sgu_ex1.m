clear all
global scalv_ ex_ recur_ recurs_ ys_ y_ exe_ lgy_ lgx_ lgr_ dsmpl_ endval_
 ...
global endo_nbr exo_nbr iy_  ykmin_  ykmax_  xkmin_  xkmax_ zkmin_ zkmax_ iter_
 ...
global dynatol_ slowc_ maxit_ valf_ ys0_ recurs0_ timing_ ct_ gstep_ Sigma_e_ fname_ lgx_orig_ord_
dsmpl_=0;
dynatol_=0.00001;
maxit_=10;
slowc_=1;
timing_=0;
ct_=0;
gstep_=1e-2;
endval_=0;rplottype_=0;
fname_ = 'sgu_ex1';
logname_ = 'sgu_ex1.log';
diary off;
warning off;
delete sgu_ex1.log;
warning on;
warning backtrace;
diary sgu_ex1.log;

iter_ = 20000;




global alpha beta delta gamma rho 


beta = 0.95;
delta = 1;
alpha = 0.3;
rho = 0;
gamma = 2;

lgy_ = 'a';
lgy_ = str2mat(lgy_,'c');
lgy_ = str2mat(lgy_,'k');
lgx_ = 'e';
lgx_orig_ord_ = [1];
endo_nbr = 3;
exo_nbr = 1;
recur_nbr = 0;
iy_ = [ 1 0 2];
temp = [ 3 4 5];
iy_ = [ iy_ ; temp ];
temp = [ 6 7 0];
iy_ = [ iy_ ; temp ];
ykmin_ = 1;
ykmax_ = 1;
xkmin_ = 0;
xkmax_ = 0;
zkmin_ = 0;
zkmax_ = 0;


% INITVAL 
valf_ = 0;
endval_=0;
ys_ = zeros(3,1);
exe_ = zeros(1,1);
ys0_ = 0;
ex0_ = 0;
recurs0_ = 0;
ys_(3)=0;
ys_(2)=0;
ys_(1)=0;
exe_(1)=0;
if exo_nbr > 0;
  ex_ =ones(iter_ + xkmin_ + xkmax_,1) * exe_';
end;


Sigma_e_ = 1;

var_list_ = [];
options.ar = 0;
options.dr_algo = 0;
options.simul_algo = 0;
options.nocorr = 1;
options.drop = 100;
options.linear = 0;
options.nofunctions = 0;
options.nomoments = 1;
options.irf = 0;
options.order = 2;
options.replic = 0;
stoch_simul(options,var_list_);


global dr_
dr_obj_ = dr_;

save sgu_ex1 dr_obj_;

diary off
