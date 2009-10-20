function [yt, j0, ir, ic]=teff(T,Nsam,istable)
%
% [yt, j0, ir, ic]=teff(T,Nsam,istable)
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

ndim = (length(size(T)));
if ndim==3,
if nargin==1,
  Nsam=size(T,3);
  istable = [1:Nsam]';
end
tmax=max(T,[],3);
tmin=min(T,[],3);
[ir, ic]=(find( (tmax-tmin)>1.e-8));
j0 = length(ir);
yt=zeros(Nsam, j0);

for j=1:j0,
  y0=squeeze(T(ir(j),ic(j),:));
  %y1=ones(size(lpmat,1),1)*NaN;
  y1=ones(Nsam,1)*NaN;
  y1(istable,1)=y0;
  yt(:,j)=y1;
end

else
tmax=max(T,[],2);
tmin=min(T,[],2);
ir=(find( (tmax-tmin)>1.e-8));
j0 = length(ir);
yt=NaN(Nsam, j0);
yt(istable,:)=T(ir,:)';

  
end
%clear y0 y1;
