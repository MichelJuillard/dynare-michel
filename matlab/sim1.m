% Copyright (C) 2001 Michel Juillard
%
function sim1
% function sim1
% performs deterministic simulations with lead or lag on one period
%
% INPUTS
%   ...
% OUTPUTS
%   ...
% ALGORITHM
%   Laffargue, Boucekkine, Juillard (LBJ)
%   see Juillard (1996) Dynare: A program for the resolution and
%   simulation of dynamic models with forward variables through the use
%   of a relaxation algorithm. CEPREMAP. Couverture Orange. 9602.
%
% SPECIAL REQUIREMENTS
%   None.
%  
%  
% part of DYNARE, copyright S. Adjemian, M. Juillard (1996-2006)
% Gnu Public License.

global M_ options_ oo_
global  iyp iyf ct_ M_ it_ c

lead_lag_incidence = M_.lead_lag_incidence;

ny = size(oo_.endo_simul,1) ;
nyp = nnz(lead_lag_incidence(1,:)) ;
nyf = nnz(lead_lag_incidence(3,:)) ;
nrs = ny+nyp+nyf+1 ;
nrc = nyf+1 ;
iyf = find(lead_lag_incidence(3,:)>0) ;
iyp = find(lead_lag_incidence(1,:)>0) ;
isp = [1:nyp] ;
is = [nyp+1:ny+nyp] ;
isf = iyf+nyp ;
isf1 = [nyp+ny+1:nyf+nyp+ny+1] ;
stop = 0 ;
iz = [1:ny+nyp+nyf];

disp (['-----------------------------------------------------']) ;
disp (['MODEL SIMULATION :']) ;
fprintf('\n') ;

it_init = M_.maximum_lag+1 ;

h1 = clock ;
for iter = 1:options_.maxit
	h2 = clock ;

	if ct_ == 0
		c = zeros(ny*options_.periods,nrc) ;
	else
		c = zeros(ny*(options_.periods+1),nrc) ;
	end

	it_ = it_init ;
	z = [oo_.endo_simul(iyp,it_-1) ; oo_.endo_simul(:,it_) ; oo_.endo_simul(iyf,it_+1)] ;
	[d1,M_.jacobia] = feval([M_.fname '_dynamic'],z,oo_.exo_simul);
	M_.jacobia = [M_.jacobia(:,iz) -d1] ;
	ic = [1:ny] ;
	icp = iyp ;
	c (ic,:) = M_.jacobia(:,is)\M_.jacobia(:,isf1) ;
	for it_ = it_init+(1:options_.periods-1)
		z = [oo_.endo_simul(iyp,it_-1) ; oo_.endo_simul(:,it_) ; oo_.endo_simul(iyf,it_+1)] ;
		[d1,M_.jacobia] = feval([M_.fname '_dynamic'],z,oo_.exo_simul);
		M_.jacobia = [M_.jacobia(:,iz) -d1] ;
		M_.jacobia(:,[isf nrs]) = M_.jacobia(:,[isf nrs])-M_.jacobia(:,isp)*c(icp,:) ;
		ic = ic + ny ;
		icp = icp + ny ;
		c (ic,:) = M_.jacobia(:,is)\M_.jacobia(:,isf1) ;
	end

	if ct_ == 1
		s = eye(ny) ;
		s(:,isf) = s(:,isf)+c(ic,1:nyf) ;
		ic = ic + ny ;
		c(ic,nrc) = s\c(:,nrc) ;
		c = bksup1(ny,nrc) ;
		c = reshape(c,ny,options_.periods+1) ;
		oo_.endo_simul(:,it_init+(0:options_.periods)) = oo_.endo_simul(:,it_init+(0:options_.periods))+options_.slowc*c ;
	else
		c = bksup1(ny,nrc) ;
		c = reshape(c,ny,options_.periods) ;
		oo_.endo_simul(:,it_init+(0:options_.periods-1)) = oo_.endo_simul(:,it_init+(0:options_.periods-1))+options_.slowc*c ;
	end

	err = max(max(abs(c./options_.scalv')));
	disp([num2str(iter) ' -	err = ' num2str(err)]) ;
	disp(['	Time of iteration 	:' num2str(etime(clock,h2))]) ;

	if err < options_.dynatol
		stop = 1 ;
		fprintf('\n') ;
		disp(['	Total time of simulation 	:' num2str(etime(clock,h1))]) ;
		fprintf('\n') ;
		disp(['	Convergency obtained.']) ;
		fprintf('\n') ;
		break
	end
end

if ~ stop
	fprintf('\n') ;
	disp(['	Total time of simulation 	:' num2str(etime(clock,h1))]) ;
	fprintf('\n') ;
	disp(['WARNING : maximum number of iterations is reached (modify options_.maxit).']) ;
	fprintf('\n') ;
end
disp (['-----------------------------------------------------']) ;
return ;

% 08/24/01 MJ added start_simul