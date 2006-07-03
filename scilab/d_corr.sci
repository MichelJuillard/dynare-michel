function []=d_corr(varargin)

if length(varargin) == 0
  for i=1:size(lgy_,1)
    varargin(i) = lgy_(i);
  end
end

m = length(varargin);
k = [];
if m > 0 then
  dyn_disp(m)
  for i=1:m
    s = varargin(i)
    k = [k; grep_exact(lgy_,s)];
  end
end

if k==[] then
  error('One of the variable specified does not exist');
end

y = y_(k,:);
y = y(:,100:$);
n = size(y,2);
y = y-mean(y,2)*ones(1,n);
v = sum(y.*y,2);
c = (y*y')./sqrt(v*v');

dyn_disp('Correlation matrix') 
  mprintf('%8s ','')
  mfprintf(fh_log,'%8s ','')
for i=1:m
  mprintf('%8s ',varargin(i))
  mfprintf(fh_log'%8s ',varargin(i))
end
mprintf('\n');
mfprintf(fh_log,'\n');
for i=1:m
  mprintf('%8s ',varargin(i))
  mfprintf(fh_log,'%8s ',varargin(i))
  for j=1:m
    mprintf('%8.4f ',c(i,j))
    mfprintf(fh_log,'%8.4f ',c(i,j))
  end
  mprintf('\n');
  mfprintf(fh_log,'\n');
end
