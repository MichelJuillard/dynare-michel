% 
% Status : main Dynare file 
% 
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ exedet_ exdet_ recur_ recurs_ 
global options_ endval_
global ys0_ recurs0_ ex0_ ct_
options_ = [];
M_.fname = 'bidon';
% 
% Some global variables initialisation
% 
global_initialization;
diary off;
warning off;

delete bidon.log;
warning on;
warning backtrace;
logname_ = 'bidon.log';
diary 'bidon.log';
erase_compiled_function('bidon_static');
erase_compiled_function('bidon_dynamic');
M_.exo_names = 'g_bar';
M_.exo_names_tex = '';
M_.endo_names = 'c';
M_.endo_names_tex = '';
M_.endo_names = strvcat(M_.endo_names,'i');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'');
M_.endo_names = strvcat(M_.endo_names,'g');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'');
M_.endo_names = strvcat(M_.endo_names,'infl');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'');
M_.endo_names = strvcat(M_.endo_names,'y');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'');
M_.endo_names = strvcat(M_.endo_names,'k');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'');
M_.endo_names = strvcat(M_.endo_names,'l');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'');
M_.endo_names = strvcat(M_.endo_names,'r');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'');
M_.endo_names = strvcat(M_.endo_names,'p1');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'');
M_.endo_names = strvcat(M_.endo_names,'q1');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'');
M_.endo_names = strvcat(M_.endo_names,'p2');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'');
M_.endo_names = strvcat(M_.endo_names,'q2');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'');
M_.endo_names = strvcat(M_.endo_names,'r0');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'');
M_.param_names = 'a';
M_.param_names_tex = '';
M_.param_names = strvcat(M_.param_names,'b');
M_.param_names_tex = strvcat(M_.param_names_tex,'');
M_.param_names = strvcat(M_.param_names,'d');
M_.param_names_tex = strvcat(M_.param_names_tex,'');
M_.param_names = strvcat(M_.param_names,'e');
M_.param_names_tex = strvcat(M_.param_names_tex,'');
M_.param_names = strvcat(M_.param_names,'f');
M_.param_names_tex = strvcat(M_.param_names_tex,'');
M_.param_names = strvcat(M_.param_names,'h');
M_.param_names_tex = strvcat(M_.param_names_tex,'');
M_.param_names = strvcat(M_.param_names,'j');
M_.param_names_tex = strvcat(M_.param_names_tex,'');
M_.param_names = strvcat(M_.param_names,'m');
M_.param_names_tex = strvcat(M_.param_names_tex,'');
M_.param_names = strvcat(M_.param_names,'n');
M_.param_names_tex = strvcat(M_.param_names_tex,'');
M_.param_names = strvcat(M_.param_names,'o');
M_.param_names_tex = strvcat(M_.param_names_tex,'');
M_.param_names = strvcat(M_.param_names,'p');
M_.param_names_tex = strvcat(M_.param_names_tex,'');
M_.param_names = strvcat(M_.param_names,'q');
M_.param_names_tex = strvcat(M_.param_names_tex,'');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.Sigma_e = zeros(1, 1);
M_.endo_nbr = 13;
M_.recur_nbr = 0;
M_.param_nbr = 12;
M_.lead_lag_incidence = [
	 0 2 5 18 0;
	 0 0 6 0 0;
	 0 0 7 0 0;
	 0 0 8 0 0;
	 1 3 9 19 20;
	 0 4 10 0 0;
	 0 0 11 0 0;
	 0 0 12 0 0;
	 0 0 13 0 0;
	 0 0 14 0 0;
	 0 0 15 0 0;
	 0 0 16 0 0;
	 0 0 17 0 0;]';
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 10;
M_.maximum_lead = 2;
M_.maximum_endo_lag = 2;
M_.maximum_endo_lead = 2;
oo_.steady_state = zeros(13, 1);
M_.maximum_exo_lag = 10;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = zeros(12, 1);
M_.params( 1 ) = (0.4*0.6);
a = M_.params( 1 );
M_.params( 2 ) = (0.3*0.6);
b = M_.params( 2 );
M_.params( 3 ) = 0.1;
d = M_.params( 3 );
M_.params( 4 ) = 0.15;
e = M_.params( 4 );
M_.params( 5 ) = 1;
f = M_.params( 5 );
M_.params( 6 ) = 0.15;
h = M_.params( 6 );
M_.params( 7 ) = 1;
j = M_.params( 7 );
M_.params( 8 ) = 1;
m = M_.params( 8 );
M_.params( 9 ) = 1;
n = M_.params( 9 );
M_.params( 10 ) = 1;
o = M_.params( 10 );
% 
% INITVAL instructions 
% 
options_.initval_file = 0;
endval_=0;
oo_.exo_steady_state( 1 ) = 0.15;
oo_.steady_state( 1 ) = 0.7;
oo_.steady_state( 2 ) = 0.15;
oo_.steady_state( 3 ) = 0.15;
oo_.steady_state( 5 ) = 1;
oo_.steady_state( 6 ) = 1;
oo_.steady_state( 7 ) = 1;
oo_.steady_state( 4 ) = 0.02;
oo_.steady_state( 8 ) = 0;
oo_.steady_state( 13 ) = oo_.steady_state(8);
oo_.steady_state( 9 ) = (2/3);
oo_.steady_state( 10 ) = (3.1/3);
oo_.steady_state( 12 ) = (4/9);
oo_.steady_state( 11 ) = (8/9);
oo_.y_simul=[oo_.steady_state*ones(1,M_.maximum_lag)];
if M_.exo_nbr > 0;
	oo_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];
end;
if M_.exo_det_nbr > 0;
	oo_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];
end;
options_.solve_algo = 2;
steady;
options_.periods = 100;
options_.simul = 1;
if ~ options_.initval_file
  make_y_;
  make_ex_;
end
disp('compiling...');
mex bidon_dynamic.cc;
oo_.endo_simul=bidon_dynamic;
var_list_=[];
var_list_ = strvcat(var_list_,'c');
rplot(var_list_);
var_list_=[];
var_list_ = strvcat(var_list_,'y');
rplot(var_list_);
save('bidon_results', 'oo_');
diary off

disp(['Total computing time : ' sec2hms(round(toc)) ]);
