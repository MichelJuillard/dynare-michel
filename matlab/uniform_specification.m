function [m,s,p1,p2] = uniform_specification(m,s,p3,p4)

% function [m,s,p1,p2] = uniform_specification(m,s,p3,p4)
% Specification of the uniform density function parameters
%
% INPUTS
%    m:      mean
%    s:      standard deviation 
%    p3:     lower bound 
%    p4:     upper bound 

% OUTPUTS
%    m:      mean
%    s:      standard deviation 
%    p1:     lower bound 
%    p2:     upper bound 
%        
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2004-2008)
% Gnu Public License.


    if ~(isnan(p3) | isnan(p4))
      p1 = p3;
      p2 = p4;
      m = (p3+p4)/2;
      s = (p4-p3)/(sqrt(12));
    else
      p1 = m-s*sqrt(3);
      p2 = m+s*sqrt(3);
    end
