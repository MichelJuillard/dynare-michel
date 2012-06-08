function [info,tmp] = homotopic_steps(endo_simul0,exo_simul0,initial_weight,step_length,pfm)

% Copyright (C) 2012 Dynare Team
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

global options_ oo_

%Set bytecode flag
bytecode_flag = options_.ep.use_bytecode;

% Set increase and decrease factors.
increase_factor = 5.0;
decrease_factor = 0.2;

% Save current state of oo_.endo_simul and oo_.exo_simul.
endo_simul = endo_simul0;
exxo_simul = exo_simul0;

initial_step_length = step_length;
max_iter = 1000/step_length;
weight   = initial_weight;
verbose  = options_.ep.debug;

reduce_step_flag = 0;

if verbose
    format long
end

% (re)Set iter.
iter = 0;
% (re)Set iter.
jter = 0;
% (re)Set weight.
weight = initial_weight;
% (re)Set exo_simul to zero.
exo_simul0 = zeros(size(exo_simul0));
while weight<1
    iter = iter+1;
    exo_simul0(2,:) = weight*exxo_simul(2,:);
    if bytecode_flag
        oo_.endo_simul = endo_simul_1;
        oo_.exo_simul = exo_simul_1;
        [flag,tmp] = bytecode('dynamic');
    else
        flag = 1;
    end
    if flag
        [flag,tmp] = solve_perfect_foresight_model(endo_simul0,exo_simul0,pfm);
    end
    info.convergence = ~flag;% Equal to one if the perfect foresight solver converged for the current value of weight.
    if verbose
        if info.convergence
            disp(['Iteration n° ' int2str(iter) ', weight is ' num2str(weight,8) ', Ok!' ])
        else
            disp(['Iteration n° ' int2str(iter) ', weight is ' num2str(weight,8) ', Convergence problem!' ])
        end
    end
    if info.convergence
        %if d<stochastic_extended_path_depth
            endo_simul0 = tmp;
            %end
        jter = jter + 1;
        if jter>3
            if verbose
                disp('I am increasing the step length!')
            end
            step_length=step_length*increase_factor;
            jter = 0;
        end
        if abs(1-weight)<options_.dynatol.x;
            break
        end
        weight = weight+step_length;
    else% Perfect foresight solver failed for the current value of weight.
        if initial_weight>0 && abs(weight-initial_weight)<1e-12% First iteration, the initial weight is too high.
            if verbose
                disp('I am reducing the initial weight!')
            end
            initial_weight = initial_weight/2;
            weight = initial_weight;
            if weight<1e-12
                endo_simul0 = endo_simul;
                exo_simul0 = exxo_simul;
                info.convergence = 0;
                info.depth = d;
                tmp = [];
                return
            end
            continue
        else% Initial weight is OK, but the perfect foresight solver failed on some subsequent iteration.
            if verbose
                disp('I am reducing the step length!')
            end
            jter = 0;
            if weight>0
                weight = weight-step_length;
            end
            step_length=step_length*decrease_factor;
            weight = weight+step_length;
            if step_length<options_.dynatol.x
                break
            end
            continue
        end
    end
    if iter>max_iter
        info = NaN;
        return
    end
end
if weight<1
    exo_simul0 = exxo_simul;
    if bytecode_flag
        oo_.endo_simul = endo_simul_1;
        oo_.exo_simul = exo_simul_1;
        [flag,tmp] = bytecode('dynamic');
    else
        flag = 1;
    end
    if flag
        [flag,tmp] = solve_perfect_foresight_model(endo_simul0,exo_simul0,pfm);
    end
    info.convergence = ~flag;
    if info.convergence
        endo_simul0 = tmp;
        return
    else
        if step_length>options_.dynatol.x
            endo_simul0 = endo_simul;
            exo_simul0 = exxo_simul;
            info.convergence = 0;
            info.depth = d;
            tmp = [];
            return
        else
            error('extended_path::homotopy: Oups! I did my best, but I am not able to simulate this model...')
        end
    end
end
