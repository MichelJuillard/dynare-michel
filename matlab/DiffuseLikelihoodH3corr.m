function [LIK lik] = DiffuseLikelihoodH3corr(T,R,Q,H,Pinf,Pstar,Y,trend,start)
% Same as DiffuseLikelihoodH3 but allows correlation between the measurement
% errors (this is not a problem with the multivariate approach). 

% Copyright (C) 2004 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global bayestopt_ options_
  
mf  = bayestopt_.mf;
pp     = size(Y,1);
mm     = size(T,1);
rr	   = size(Q,1);	
smpl   = size(Y,2);
T   = cat(1,cat(2,T,zeros(mm,pp)),zeros(pp,mm+pp));
R   = cat(1,cat(2,R,zeros(mm,pp)),cat(2,zeros(pp,rr),eye(pp)));
Q   = cat(1,cat(2,Q,zeros(rr,pp)),cat(2,zeros(pp,rr),H));
if size(Pinf,1) % Otherwise Pinf = 0 (no unit root) 
	Pinf   = cat(1,cat(2,Pinf,zeros(mm,pp)),zeros(pp,mm+pp));
end
Pstar  = cat(1,cat(2,Pstar,zeros(mm,pp)),cat(2,zeros(pp,mm),H));
a      = zeros(mm+pp,1);
QQ     = R*Q*transpose(R);
t      = 0;
lik    = zeros(smpl,1);
notsteady 	= 1;
crit      	= options_.kalman_tol;
newRank	  	= rank(Pinf,crit);

while rank(Pinf,crit) & t < smpl %% Matrix Finf is assumed to be zero
	t = t+1;
	for i=1:pp
		v(i) 	= Y(i,t)-a(mf(i))-a(mm+i)-trend(i,t);
		Fstar 	= Pstar(mf(i),mf(i))+Pstar(mm+i,mm+i);
		Finf	= Pinf(mf(i),mf(i));
		Kstar 	= Pstar(:,mf(i))+Pstar(:,mm+i);
		if Finf > crit
			Kinf	= Pinf(:,mf(i));
			a		= a + Kinf*v(i)/Finf;
			Pstar	= Pstar + Kinf*transpose(Kinf)*Fstar/(Finf*Finf) - ...
						(Kstar*transpose(Kinf)+Kinf*transpose(Kstar))/Finf;
			Pinf	= Pinf - Kinf*transpose(Kinf)/Finf;
			lik(t) 	= lik(t) + log(Finf);
		else %% Note that : (1) rank(Pinf)=0 implies that Finf = 0, (2) outside this loop (when for some i and t the condition
			 %% rank(Pinf)=0 is satisfied we have P = Pstar and F = Fstar and (3) Finf = 0 does not imply that
			 %% rank(Pinf)=0. [stéphane,11-03-2004].	  
			if rank(Pinf) == 0
				lik(t)	= lik(t) + log(Fstar) + v(i)*v(i)/Fstar;
			end
			a 		= a + Kstar*v(i)/Fstar;
			Pstar	= Pstar - Kstar*transpose(Kstar)/Fstar;					
		end
		oldRank = rank(Pinf,crit);
		a 		= T*a;
		Pstar 	= T*Pstar*transpose(T)+QQ;
		Pinf	= T*Pinf*transpose(T);
		newRank = rank(Pinf,crit);
		if oldRank ~= newRank
			disp('DiffuseLiklihoodH3 :: T does influence the rank of Pinf!')	
		end		 		
	end
end
if t == smpl                                                           
  error(['There isn''t enough information to estimate the initial' ... 
	 ' conditions of the nonstationary variables']);                   
end                                                                    
while notsteady & t < smpl
	t = t+1;
	for i=1:pp
		v(i) = Y(i,t) - a(mf(i)) - trend(i,t) -a(mm+i);
		Fi   = Pstar(mf(i),mf(i))+Pstar(mm+i,mm+i);
		if Fi > crit
			Ki		= Pstar(:,mf(i))+Pstar(:,mm+i);
            a		= a + Ki*v(i)/Fi;
            Pstar 	= Pstar - Ki*transpose(Ki)/Fi;
            lik(t) 	= lik(t) + log(Fi) + v(i)*v(i)/Fi;
		end
	end	
    oldP 		= Pstar;
    a 			= T*a;
    Pstar 		= T*Pstar*transpose(T) + QQ;
    notsteady 	= ~(max(max(abs(Pstar-oldP)))<crit);
end
while t < smpl
	t = t+1;
	for i=1:pp
		v(i) = Y(i,t) - a(mf(i)) - trend(i,t) - a(mm+i);
		Fi   = Pstar(mf(i),mf(i))+Pstar(mm+i,mm+i);
		if Fi > crit
			Ki 		= Pstar(:,mf(i))+Pstar(:,mm+i);
            a 		= a + Ki*v(i)/Fi;
            Pstar 	= Pstar - Ki*transpose(Ki)/Fi;
            lik(t) 	= lik(t) + log(Fi) + v(i)*v(i)/Fi;
		end
	end	
    a = T*a;
end
% adding log-likelihhod constants
lik = (lik + pp*log(2*pi))/2;

LIK = sum(lik(start:end)); % Minus the log-likelihood.
