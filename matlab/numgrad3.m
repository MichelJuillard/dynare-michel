function [g, badg, f0, f1, f2] = numgrad3(fcn,x,varargin)
% Computes the gradient of the objective function fcn using a three points
% formula if possible.
%
% Adapted from Sims' numgrad routine.
%
% See section 25.3.4 in Abramovitz and Stegun (1972, Tenth Printing, December) Handbook of Mathematical Functions.
% http://www.math.sfu.ca/~cbm/aands/ 
%
% part of DYNARE, copyright Dynare Team (2008)
% Gnu Public License.

f0 = NaN;
f1 = NaN;
f2 = NaN;

delta = 1e-6;
n=length(x);
tvec=delta*eye(n);
g=zeros(n,1);

[f0,cost_flag] = feval(fcn, x, varargin{:});

badg=0;
goog=1;
scale=1;
for i=1:n
   if size(x,1)>size(x,2)
      tvecv=tvec(i,:);
   else
      tvecv=tvec(:,i);
   end
   [f1,cost_flag1] = feval(fcn, x+scale*transpose(tvecv), varargin{:});
   [f2,cost_flag2] = feval(fcn, x-scale*transpose(tvecv), varargin{:});
   if cost_flag1 && cost_flag2
       g0 = (f1 - f2) / (2*scale*delta);
   else
       if cost_flag1
           g0 = (f1-f0) / (scale*delta);
       elseif cost_flag2
           g0 = (f0-f2) / (scale*delta);
       else
           goog=0;
       end
   end
   if goog && abs(g0)< 1e15 
      g(i)=g0;
   else
       disp('bad gradient ------------------------')
       g(i)=0;
       badg=1;
   end
end