function gsa_ = gsa_sdp_dyn(varargin)

% function gsa_ =
% gsa_sdp_dyn(y, x, tvp, nvr, ialias, nvrT, ifig, fname, pnames, yi_anal)
%             1  2  3    4    5       6     7     8      9       10
% INPUTS
% y, the output
% x, the inputs
% tvp: tipe of estimation
%      0 = load previous estimation
%      1 = RW for all parameters
%      2 = IRW for all parameters
%      -1 = RW for all parameters, fixed nvr's (no estimation)
%      -2 = IRW for all parameters, fixed nvr's (no estimation)
%      an array of length(x) specifies the process for each parameter
% nvr  initial nvr value for estimation or fixed nvr is no estimation
% ialias to be used with LPTAU samples < 1024
% nvrT maximum nvr allowable
% ifig: 1 plot figures, 0 no plots
% fname, file name to save the analysis
% pnames: names of the params
% yi_anal, analytical value of the first order HDMR terms (optional)
%
% OUTPUT
%   gsa_.univariate.f
%   gsa_.univariate.fs
%   gsa_.univariate.fses
%   gsa_.univariate.si
%   gsa_.univariate.si_std
%   gsa_.multivariate.f
%   gsa_.multivariate.fs
%   gsa_.multivariate.fses
%   gsa_.multivariate.si
%   gsa_.multivariate.si_std
%   gsa_.multivariate.stat
%   gsa_.x0
%   gsa_.xx
%   gsa_.y
% 
% f, function estimates
% fs, sorted function estimates
% fses, sorted standard error of function estimates
% si, sensitivity indices
% si_std, st. error of sensitivity indices
% stat, euristic tstat for significance of f
% xx, transformed inputs used (rank-transformed)
% x0, original inputs
% y, output
%
% Part of the Sensitivity Analysis Toolbox for DYNARE
%
% Written by Marco Ratto, 2006
% Joint Research Centre, The European Commission,
% (http://eemc.jrc.ec.europa.eu/),
% marco.ratto@jrc.it 
%
% Disclaimer: This software is not subject to copyright protection and is in the public domain. 
% It is an experimental system. The Joint Research Centre of European Commission 
% assumes no responsibility whatsoever for its use by other parties
% and makes no guarantees, expressed or implied, about its quality, reliability, or any other
% characteristic. We would appreciate acknowledgement if the software is used.
% Reference:
% M. Ratto, Global Sensitivity Analysis for Macroeconomic models, MIMEO, 2006.
%

if exist('sdr','file')==6 | exist('sdr','file')==2,
  gsa_ = gsa_sdp_fn(varargin{:});
else
  disp('Download the SDP mapping routines at:')
  disp('http://eemc.jrc.ec.europa.eu/softwareDYNARE-Dowload.htm')
  disp(' ' )
  error('SDP mapping routines missing!')

end
