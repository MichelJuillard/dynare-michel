function y = solve_one_boundary(fname, y, x, y_index_eq, nze, periods, is_linear, Block_Num, y_kmin, maxit_, solve_tolf, lambda, cutoff, simulation_method, forward_backward)
% Computes the deterministic simulation of a block of equation containing
% lead or lag variables 
%
% INPUTS
%   fname               [string]        name of the file containing the block
%                                       to simulate
%   y                   [matrix]        All the endogenous variables of the model
%   x                   [matrix]        All the exogenous variables of the model
%   y_index_eq          [vector of int] The index of the endogenous variables of
%                                       the block
%   nze                 [integer]       number of non-zero elements in the
%                                       jacobian matrix
%   periods             [integer]       number of simulation periods
%   is_linear           [integer]       if is_linear=1 the block is linear
%                                       if is_linear=0 the block is not linear
%   Block_Num           [integer]       block number
%   y_kmin              [integer]       maximum number of lag in the model
%   maxit_              [integer]       maximum number of iteration in Newton
%   solve_tolf          [double]        convergence criteria
%   lambda              [double]        initial value of step size in
%   Newton
%   cutoff              [double]        cutoff to correct the direction in Newton in case
%                                       of singular jacobian matrix
%   simulation_method   [integer]       linear solver method used in the
%                                       Newton algorithm : 
%                                            - 0 sprse LU
%                                            - 2 GMRES
%                                            - 3 BicGStab
%   forward_backward    [integer]       The block has to be solve forward
%                                       (1) or backward (0)
%
% OUTPUTS
%   y                  [matrix]         All endogenous variables of the model      
%  
% ALGORITHM
%   Newton with LU or GMRES or BicGstab
%    
% SPECIAL REQUIREMENTS
%   none.
%  

% Copyright (C) 1996-2008 Dynare Team
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


    global oo_;
    Blck_size=size(y_index_eq,2);
    g2 = [];
    g3 = [];
    correcting_factor=0.01;
    luinc_tol=1e-10;
    max_resa=1e100;
    reduced = 0;
    if(forward_backward)
        incr = 1;
        start = y_kmin+1;
        finish = periods+y_kmin;
    else
        incr = -1;
        start = periods+y_kmin;
        finish = y_kmin+1;
    end
    
    for it_=start:incr:finish
       cvg=0;
       iter=0;
       g1=spalloc( Blck_size, Blck_size, nze);
       %Per_y_=it_*Blck_size;
       while ~(cvg==1 | iter>maxit_),
           [r, g1] = feval(fname, y, x, it_, 0, g1, g2, g3);
           %r
           if(~isreal(r))
              max_res=(-(max(max(abs(r))))^2)^0.5;
           else
              max_res=max(max(abs(r)));
           end;
           if(iter>0)
             if(~isreal(max_res) | isnan(max_res) | (max_resa<max_res && iter>1))
               if(isnan(max_res))
                 detJ=det(g1a);
                 if(abs(detJ)<1e-7)
                   max_factor=max(max(abs(g1a)));
                   ze_elem=sum(diag(g1a)<cutoff);
                   disp([num2str(full(ze_elem),'%d') ' elements on the Jacobian diagonal are below the cutoff (' num2str(cutoff,'%f') ')']);
                   if(correcting_factor<max_factor)
                     correcting_factor=correcting_factor*4;
                     disp(['The Jacobain matrix is singular, det(Jacobian)=' num2str(detJ,'%f') '.']);
                     disp(['    trying to correct the Jacobian matrix:']);
                     disp(['    correcting_factor=' num2str(correcting_factor,'%f') ' max(Jacobian)=' num2str(full(max_factor),'%f')]);
                     dx = -r/(g1+correcting_factor*speye(Blck_size));
                     y(it_,y_index_eq)=ya_save+lambda*dx;
                     continue;
                   else
                     disp('The singularity of the jacobian matrix could not be corrected');
                     return;
                   end;
                 end;
               elseif(lambda>1e-8)
                 lambda=lambda/2;
                 reduced = 1;
                 disp(['reducing the path length: lambda=' num2str(lambda,'%f')]);
                 y(it_,y_index_eq)=ya_save+lambda*dx;
                 continue;
               else
                 if(cutoff == 0)
                   fprintf('Error in simul: Convergence not achieved in block %d, at time %d, after %d iterations.\n Increase "options_.maxit_".\n',Block_Num, it_, iter);
                 else
                   fprintf('Error in simul: Convergence not achieved in block %d, at time %d, after %d iterations.\n Increase "options_.maxit_" or set "cutoff=0" in model options.\n',Block_Num, it_, iter);
                 end;
                 oo_.deterministic_simulation.status = 0;
                 oo_.deterministic_simulation.error = max_res;
                 oo_.deterministic_simulation.iterations = iter;
                 oo_.deterministic_simulation.block(Block_Num).status = 0;% Convergency failed.
                 oo_.deterministic_simulation.block(Block_Num).error = max_res;
                 oo_.deterministic_simulation.block(Block_Num).iterations = iter;
                 return;
               end;
             else
               if(lambda<1)
                 lambda=max(lambda*2, 1);
               end;
             end;
           end;
           ya = y(it_,y_index_eq)';
           ya_save=ya;
           g1a=g1;
           if(simulation_method==0),
              dx = (-r/g1)';
              ya = ya + lambda*dx;
              y(it_,y_index_eq)=ya';
           elseif(simulation_method==2),
              flag1=1;
              while(flag1>0)
                 [L1, U1]=luinc(g1,luinc_tol);
                 [za,flag1] = gmres(g1,-r',Blck_size,1e-6,Blck_size,L1,U1);
                 if (flag1>0 | reduced)
                    if(flag1==1)
                       disp(['Error in simul: No convergence inside GMRES after ' num2str(iter,'%6d') ' iterations, in block' num2str(Block_Num,'%3d')]);
                    elseif(flag1==2)
                       disp(['Error in simul: Preconditioner is ill-conditioned, in block' num2str(Block_Num,'%3d')]);
                    elseif(flag1==3)
                       disp(['Error in simul: GMRES stagnated (Two consecutive iterates were the same.), in block' num2str(Block_Num,'%3d')]);
                    end;
                    luinc_tol = luinc_tol/10;
                    reduced = 0;
                 else
                    dx = za - ya;
                    ya = ya + lambda*dx;
                    y(it_,y_index_eq)=ya';
                 end;
              end;
           elseif(simulation_method==3),
              flag1=1;
              while(flag1>0)
                 [L1, U1]=luinc(g1,luinc_tol);
                 [za,flag1] = bicgstab(g1,-r',1e-7,Blck_size,L1,U1);
                 if (flag1>0 | reduced)
                    if(flag1==1)
                       disp(['Error in simul: No convergence inside BICGSTAB after ' num2str(iter,'%6d') ' iterations, in block' num2str(Block_Num,'%3d')]);
                    elseif(flag1==2)
                       disp(['Error in simul: Preconditioner is ill-conditioned, in block' num2str(Block_Num,'%3d')]);
                    elseif(flag1==3)
                       disp(['Error in simul: GMRES stagnated (Two consecutive iterates were the same.), in block' num2str(Block_Num,'%3d')]);
                    end;
                    luinc_tol = luinc_tol/10;
                    reduced = 0;
                 else
                    dx = za - ya;
                    ya = ya + lambda*dx;
                    y(it_,y_index_eq)=ya';
                 end;
              end;
           end;
           if(is_linear)
              cvg = 1;
           else
              cvg=(max_res<solve_tolf);
           end;
           iter=iter+1;
       end
       if cvg==0
           if(cutoff == 0)
               fprintf('Error in simul: Convergence not achieved in block %d, at time %d, after %d iterations.\n Increase "options_.maxit_\".\n',Block_Num, it_,iter);
           else
               fprintf('Error in simul: Convergence not achieved in block %d, at time %d, after %d iterations.\n Increase "options_.maxit_" or set "cutoff=0" in model options.\n',Block_Num, it_,iter);
           end;
           oo_.deterministic_simulation.status = 0;
           oo_.deterministic_simulation.error = max_res;
           oo_.deterministic_simulation.iterations = iter;
           oo_.deterministic_simulation.block(Block_Num).status = 0;% Convergency failed.
           oo_.deterministic_simulation.block(Block_Num).error = max_res;
           oo_.deterministic_simulation.block(Block_Num).iterations = iter;
           return;
       end
    end
    oo_.deterministic_simulation.status = 1;
    oo_.deterministic_simulation.error = max_res;
    oo_.deterministic_simulation.iterations = iter;
    oo_.deterministic_simulation.block(Block_Num).status = 1;% Convergency failed.
    oo_.deterministic_simulation.block(Block_Num).error = max_res;
    oo_.deterministic_simulation.block(Block_Num).iterations = iter;
    return;