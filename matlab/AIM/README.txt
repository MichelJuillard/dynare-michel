AIM Subsystem cosists of 

A function dynAIMsolver1(jacobia_,M_,dr) which maps Dynare to Gary Anderson AIM package subsystem and derives the solution for gy=dr.hgx and gu=dr.hgu from the AIM outputs (where "1" in the title is for 1st order solver).

A subset of routines from Gary Anderson AIM package starting with SP... needed to compute and solve system passed on and returned by dynAIMsolver1.
( see  http://www.federalreserve.gov/Pubs/oss/oss4/aimindex.html  )

The path to the AIM directory,if exists, is added by dynare_config.m using addpath.

DR1 tries to invoke AIM if options_.aim_solver == 1 is set and, if not check only, and if 1st order only: 
    if (options_.aim_solver == 1) && (task == 0) && (options_.order == 1) 

for start, options_.aim_solver = 0 is set by default in global_initialization.m so that system uses mjdgges by default. 

If AIM is to be used, options_.aim_solver = 1 needs to be set either in the model <>.mod file, before invoking, estimate and/or stoch_simul, or by issuing appropriate command for estimate and/or stoch_simul. 

NOTE: in the current implementation, as of July 2008, handling of exceptions is rather fundamental and, in particular, when Blanchard and Kahn conditions are not met, only a large penalty value 1.0e+8 is being set.

Hence, system may not coverge ot resluts may not be accurate if there were many messages like 
Error in AIM: aimcode=4 : Aim: too few big roots
or
Error in AIM: aimcode=3 : Aim: too many big roots
especially close to the point of convergence.

However, if other exceptions occur and aimcode (see codes below) is higher than 5, the system resets options_.aim_solver = 0 and tries to use mjdgges instead.


APPENDIX

% AIM System is given as a sum: 
% i.e. for i=-$...+&   SUM(Hi*xt+i)= £*zt, t = 0, . . . ,?
% and its input as single array of matrices: [H-$...  Hi ... H+&]
% and its solution as xt=SUM( Bi*xt+i) + @*£*zt for i=-$...-1 
% with the output in form bb=[B-$...  Bi ... B-1] and @=inv(Ho+H1*B-1) 
% Dynare jacobian = [fy'-$...  fy'i ... fy'+&  fu'] 
% where [fy'-$...  fy'i ... fy'+&]=[H-$...  Hi ... H+&] and fu'= £


%function [dr,aimcode]=dynAIMsolver1(jacobia_,M_,dr)
% INPUTS
%   jacobia_   [matrix]           1st order derivative of the model 
%   dr         [matlab structure] Decision rules for stochastic simulations.
%   M_         [matlab structure] Definition of the model.           
%    
% OUTPUTS
%   dr         [matlab structure] Decision rules for stochastic simulations.
%   aimcode    [integer]          1: the model defines variables uniquely
%   aimcode is resolved in AIMerr as
%      (c==1)  e='Aim: unique solution.';
%      (c==2)  e='Aim: roots not correctly computed by real_schur.';
%      (c==3)  e='Aim: too many big roots.';
%      (c==35) e='Aim: too many big roots, and q(:,right) is singular.';
%      (c==4)  e='Aim: too few big roots.';
%      (c==45) e='Aim: too few big roots, and q(:,right) is singular.';
%      (c==5)  e='Aim: q(:,right) is singular.';
%      (c==61) e='Aim: too many exact shiftrights.';
%      (c==62) e='Aim: too many numeric shiftrights.';
%      else    e='Aimerr: return code not properly specified';
%    
% SPECIAL REQUIREMENTS
% Dynare use: 
%       1) the lognormal block in DR1 is being invoked for some models and changing
%       values of ghx and ghy. We need to return the AIM output
%       values before that block and run the block with the current returned values
%       of gy (i.e. dr.ghx) and gu (dr.ghu) if it is needed even when the AIM is used  
%       (it does not depend on mjdgges output).
%       
%       2) for forward looking models, passing into dynAIMsolver 
%		aa={Q'|1}*jacobia_ can produce ~ one order closer  
%       results to the Dynare solutiion then when if plain jacobia_ is passed, 
%       i.e. diff < e-14 for aa and diff < *e-13 for jacobia_ if Q' is used.  
%
% GP July 2008  
% part of DYNARE, copyright Dynare Team (1996-2008)
% Gnu Public License.
