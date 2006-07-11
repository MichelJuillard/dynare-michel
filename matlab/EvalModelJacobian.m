function J = EvalModelJacobian()
% Evaluates the Jacobian matrix associated to a model.
% 
% INPUTS 
%   None.
%  
% OUTPUTS 
%   o J a matlab structure: J.Eq1.<VARIABLE_NAME>.[<lag_><lead_>h]or[current]
%
% ALGORITHM 
%          
%
% SPECIAL REQUIREMENTS
%   None.
%  
%  
% part of DYNARE, copyright S. Adjemian, M. Juillard (2006)
% Gnu Public License.
global options_ M_ oo_   

klen = M_.maximum_lag + M_.maximum_lead + 1;
iyv  = M_.lead_lag_incidence';
% [(t-M_.maximum_lag);...;(t-1);t;(t+1);...;t+M_.maximum_lead] for each
% column (variable) of M_.lead_lag_incidence.
iyv  = iyv(:);
iyr0 = find(iyv) ;
if M_.exo_nbr == 0
  oo_.exo_steady_state = [] ;
end
tx = oo_.exo_simul;
dr = oo_.dr;

z = repmat(dr.ys,1,klen);
z = z(iyr0) ;
[junk,jacobia_] = feval([M_.fname '_dynamic'],z,tx);

for eq = 1:M_.endo_nbr
  for var = 1:M_.endo_nbr
    for h = -M_.maximum_lag:1:M_.maximum_lead
      linee = h+M_.maximum_lag+1;
      if h<0
        name = ['lag_' int2str(abs(h))];
      elseif h == 0
        name = ['current'];
      else
        name = ['lead_' int2str(h)];
      end
      if M_.lead_lag_incidence(linee,var) == 0
        eval(['J.Eq' int2str(eq) '.' deblank(M_.endo_names(var,:)) '.' name ' = 0;'])
      else
        idx = M_.lead_lag_incidence(linee,var);
        eval(['J.Eq' int2str(eq) '.' deblank(M_.endo_names(var,:)) '.' name ' = jacobia_(eq,idx);'])
      end
    end
  end
end