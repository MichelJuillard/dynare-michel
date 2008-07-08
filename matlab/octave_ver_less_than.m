function r = octave_ver_less_than(verstr)
% function r = octave_ver_less_than(verstr)
%
% Returns 1 if current Octave version is strictly older than
% the one given in argument.
%
% Note that this function will fail under Matlab.
%
% INPUTS
%    verstr: a string of the format 'x.y' or 'x.y.z'
%    
% OUTPUTS
%    r: 0 or 1
%
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2008)
% Gnu Public License.

  cur_verstr = version();
  
  r = get_ver_numeric(cur_verstr) < get_ver_numeric(verstr);


function x = get_ver_numeric(verstr)
  nums = sscanf(verstr, '%d.%d.%d')';
  if length(nums) < 3
    nums(3) = 0;
  end
  x = nums * [1; 0.01; 0.0001 ];
