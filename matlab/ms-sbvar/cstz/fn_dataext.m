function [xdsube,Brow,Erow] = fn_dataext(Byrqm,Eyrqm,xdatae)
% xdsube = dataext(Byrqm,Eyrqm,xdatae)
%      Extract subset of xdatae, given Byrqm and Eyrqm
%
% Byrqm: [year quarter(month)]: beginning year and period.  If Byqm(2)=0, we get annual data.
% Eyrqm: [year period]:  end year and period.  If Byrqm(2)=0, it must be Eyrqm(2)=0.
% xdatae:  all data (some of which may be NaN) with the first 2 columns indicating years and periods.
%----------
% xdsube:  subset of xdatae.
% Brow:  the row number in xdatee that marks the first row of xdsube.
% Erow:  the row number in xdatee that marks the last row of xdsube.
%
% Tao Zha, April 2000

if (Byrqm(2)==0) && (Eyrqm(2)~=0)
   error('If annual data, make sure both Byrqm(2) and Eyrqm(2) are zero')
end

Brow = min(find(xdatae(:,1)==Byrqm(1)));
if isempty(Brow)
   error('Byrqm is outside the date range indicated by xdatae(:,1:2)!')
end
if Byrqm(2)>0
   nadt=Byrqm(2)-xdatae(Brow,2);
   if nadt<0
      error('Byrqm is outside the date range indicated by xdatae(:,1:2)!')
   end
   Brow=Brow+nadt;
end
%
Erow = min(find(xdatae(:,1)==Eyrqm(1)));
if isempty(Erow)
   error('Eyrqm is outside the date range indicated by xdatae(:,1:2)!')
end
if Eyrqm(2)>0
   nadt2=Eyrqm(2)-xdatae(Erow,2);
   if nadt<0
      error('Eyrqm is outside the date range indicated by xdatae(:,1:2)!')
   end
   Erow=Erow+nadt2;
end
%
xdsube = xdatae(Brow:Erow,:);   % with dates
