function warning_config()
% Activates useful warnings
%
% INPUTS
%   none
%             
% OUTPUTS
%   none
%        
% SPECIAL REQUIREMENTS
%   none
%  
% part of DYNARE, copyright Dynare Team (2008)
% Gnu Public License.

    warning on;
    if exist('OCTAVE_VERSION')
        warning('off', 'Octave:separator-insert');
        warning('off', 'Octave:matlab-incompatible');
        warning('off', 'Octave:single-quote-string');
        warning('off', 'Octave:missing-semicolon');
        warning('off', 'Octave:empty-list-elements');
        warning('off', 'Octave:num-to-str');
        warning('off', 'Octave:resize-on-range-error');
        warning('off', 'Octave:str-to-num');
        warning('off', 'Octave:string-concat');
        warning('off', 'Octave:variable-switch-label');
        warning('off', 'Octave:fortran-indexing');
    else
        warning backtrace;
    end
