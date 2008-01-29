function dsample(s1,s2)
		
% function dsample(s1,s2)
% This optional command permits to reduce the number of periods considered in following output commands.
% If only one argument is provided, output is from period 1 to the period specified in the DSAMPLE command. 
% If two arguments are present output is done for the interval between the two periods.
% DSAMPLE without arguments reset the sample to the one specified by PERIODS
%
% INPUTS
%    s1:      first period
%    s2:      last period
%    
% OUTPUTS
%    none
%    
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2001-2008)
% Gnu Public License.


global options_

options_.smpl = zeros(2,1) ;

if s1 > options_.periods | s2 > options_.periods
  t = ['DYNARE dsample error: one of the arguments is larger than the one' ...
       ' specified in PERIODS'];
  error(t);
end

if nargin == 0
	options_.smpl(1) = 1 ;
	options_.smpl(2) = options_.periods ;
elseif nargin == 1
	options_.smpl(1) = 1 ;
	options_.smpl(2) = s1 ;
else
	options_.smpl(1) = s1 ;
	options_.smpl(2) = s2 ;
end

% 02/23/01 MJ added error checking