function moments=compute_model_moments(dr,M_,options,_)
%
% INPUTS
%    dr:        structure describing model solution
%    M_:   structure of Dynare options
%     options_
%    
% OUTPUTS
%    moments: a cell array containing requested moments
%        
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2008)
% Gnu Public License.

    moments = th_autocovariances(dr,options_.varlist,M_,options_);
