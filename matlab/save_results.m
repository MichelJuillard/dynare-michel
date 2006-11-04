function save_results(x,s_name,names)
% function save_results(x,s_name,names)
% save results in appropriate structure
% INPUT
%   x: matrix to be saved column by column
%   s_name: name of the structure where to save the results
%   names: names of the individual series
% OUTPUT
%   none
% ALGORITHM
%   none
% SPECIAL REQUIREMENT
%   none
%    
% part of DYNARE, copyright S. Adjemian, M. Juillard (2006)
% Gnu Public License.
  global oo_
  
  for i=1:size(x,2)
    eval([s_name deblank(names(i,:)) '= x(:,i);']);
  end