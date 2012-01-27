function [info,tmp] = homotopic_steps(initial_weight,step_length)
global oo_ options_ M_

% Set increase and decrease factors.
increase_factor = 5.0;
decrease_factor = 0.2;

% Save current state of oo_.endo_simul and oo_.exo_simul.
endo_simul = oo_.endo_simul;
exxo_simul = oo_.exo_simul;

[idx,jdx] = find(abs(exxo_simul)>1e-12);
idx = unique(idx);

if idx(1)-2
    % The first scalar in idx must be equal to two.
    error('extended_path::homotopy: Something is wrong in oo_.exo_simul!')
end

stochastic_extended_path_depth = length(idx);

if stochastic_extended_path_depth>1
    for i=2:length(idx)
        if (idx(i)-idx(i-1)-1)
            error('extended_path::homotopy: Something is wrong in oo_.exo_simul!')
        end
    end
end

initial_step_length = step_length;
max_iter = 1000/step_length;
weight   = initial_weight;
verbose  = options_.ep.debug;

reduce_step_flag = 0;

if verbose
    format long
end

for d=1:stochastic_extended_path_depth
    % (re)Set iter.
    iter = 0;
    % (re)Set iter.
    jter = 0;
    % (re)Set weight.
    weight = initial_weight;
    % (re)Set exo_simul to zero.
    oo_.exo_simul = zeros(size(oo_.exo_simul));
    while weight<1
        iter = iter+1;
        oo_.exo_simul(idx(d),:) = weight*exxo_simul(idx(d),:);
        if d>1
            oo_.exo_simul(1:idx(d-1),:) = exxo_simul(1:idx(d-1),:);
        end
        [flag,tmp] = bytecode('dynamic');
        info.convergence = ~flag;% Equal to one if the perfect foresight solver converged for the current value of weight.
        if verbose
            if info.convergence
                disp(['Iteration n° ' int2str(iter) ', weight is ' num2str(weight,8) ', Ok!' ])
            else
                disp(['Iteration n° ' int2str(iter) ', weight is ' num2str(weight,8) ', Convergence problem!' ])
            end
        end
        if info.convergence
            if d<stochastic_extended_path_depth
                oo_.endo_simul = tmp;
            end
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
                    oo_.endo_simul = endo_simul;
                    oo_.exo_simul = exxo_simul;
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
        oo_.exo_simul = exxo_simul;
        [flag,tmp] = bytecode('dynamic');
        info.convergence = ~flag;
        if info.convergence
            oo_.endo_simul = tmp;
            return
        else
            if step_length>options_.dynatol.x
                oo_.endo_simul = endo_simul;
                oo_.exo_simul = exxo_simul;
                info.convergence = 0;
                info.depth = d;
                tmp = [];
                return
            else
                error('extended_path::homotopy: Oups! I did my best, but I am not able to simulate this model...')
            end
        end
    end
end