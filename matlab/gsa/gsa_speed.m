function [tadj, iff] = gsa_speed(A,B,mf,p),
% [tadj, iff] = gsa_speed(A,B,mf,p),
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

nvar=length(mf);
nstate= size(A,1);
nshock = size(B,2);
nrun=size(B,3);

iff=zeros(nvar,nshock,nrun);
tadj=iff;
disp('Computing speed of adjustement ...')
h = dyn_waitbar(0,'Speed of adjustement...');

for i=1:nrun,
  irf=zeros(nvar,nshock);
  a=squeeze(A(:,:,i));
  b=squeeze(B(:,:,i));
  IFF=inv(eye(nstate)-a)*b;
  iff(:,:,i)=IFF(mf,:);
  IF=IFF-b;
  
  t=0;
  while any(any(irf<0.5))
    t=t+1;
    IFT=((eye(nstate)-a^(t+1))*inv(eye(nstate)-a))*b-b;
    irf=IFT(mf,:)./(IF(mf,:)+eps);
    irf = irf.*(abs(IF(mf,:))>1.e-7)+(abs(IF(mf,:))<=1.e-7);
    %irf=ft(mf,:);
    tt=(irf>0.5).*t;
    tadj(:,:,i)=((tt-tadj(:,:,i))==tt).*tt+tadj(:,:,i);
  end
  dyn_waitbar(i/nrun,h)
end
disp(' ')
disp('.. done !')
dyn_waitbar_close(h)
