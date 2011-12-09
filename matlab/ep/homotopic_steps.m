function [info,tmp] = homotopic_steps(initial_weight,step_length,time)
global oo_ options_ M_

% Save current state of oo_.endo_simul and oo_.exo_simul.
endo_simul = oo_.endo_simul;
exxo_simul = oo_.exo_simul;

% Reset exo_simul to zero.
oo_.exo_simul = zeros(size(oo_.exo_simul));


initial_step_length = step_length;
max_iter = 1000/step_length;
weight   = initial_weight;
verbose  = 1;
iter     = 0;
ctime    = 0;

reduce_step_flag = 0;

if verbose
    format long
end

homotopy_1 = 1; % Only innovations are rescaled. Starting from weight equal to initial_weight.
homotopy_2 = 0; % Only innovations are rescaled. Starting from weight equal to zero.

disp(' ')

if homotopy_1
    while weight<1
        iter = iter+1;
        oo_.exo_simul(2,:) = weight*exxo_simul(2,:);
        t0 = tic;
        [flag,tmp] = bytecode('dynamic');
        TeaTime = toc(t0);
        ctime = ctime+TeaTime;
        %old_weight = weight;
        info.convergence = ~flag;
        if verbose
            if ~info.convergence
                disp(['Iteration n° ' int2str(iter) ', weight is ' num2str(weight,8) ', Convergence problem!' ])
            else
                disp(['Iteration n° ' int2str(iter) ', weight is ' num2str(weight,8) ', Ok!' ])
            end
        end
        if ~info.convergence
            if abs(weight-initial_weight)<1e-12% First iterations.
                if verbose
                    disp('I am reducing the initial weight!')
                end
                initial_weight = initial_weight/1.1;
                weight = initial_weight;
                if weight<1/4
                    homotopy_1 = 0;
                    homotopy_2 = 1;
                    break
                end
                continue
            else% A good initial weight has been obtained. In case of convergence problem we have to reduce the step length.
                if verbose
                    disp('I am reducing the step length!')
                end
                weight = weight-step_length;
                step_length=step_length/10;
                weight = weight+step_length;
                if 10*step_length<options_.dynatol
                    homotopy_1 = 0;
                    homotopy_2 = 0;
                    break
                end
                continue
            end
        else
            oo_.endo_simul = tmp;
            info.time = ctime;
            if abs(1-weight)<=1e-12;
                homotopy_1 = 0;
                homotopy_2 = 0;
                break
            end
            weight = weight+step_length;
            %step_length = initial_step_length;
        end
        if iter>max_iter
            info = NaN;
            return
        end
    end
    if weight<1 && homotopy_1
        oo_.exo_simul(2,:) = exxo_simul(2,:);
        t0 = tic; 
        [flag,tmp] = bytecode('dynamic');
        TeaTime = toc(t0);
        ctime = ctime+TeaTime;
        info.convergence = ~flag;
        info.time = ctime;
        if info.convergence
            oo_.endo_simul = tmp;
            homotopy_1 = 0;
            homotopy_2 = 0;
            return
        else
            if step_length>1e-12
                if verbose
                    disp('I am reducing step length!')
                end
                step_length=step_length/2;
            else
                weight = initial_weight;
                step_length = initial_step_length;
                info = NaN;
                homotopy_2 = 1;
                homotopy_1 = 0;
            end
        end
    end
end

iter   = 0;
weight = 0; 

if homotopy_2
    while weight<1
        iter = iter+1;
        oo_.exo_simul(2,:) = weight*exxo_simul(2,:);
        if time==1
            oo_.endo_simul = repmat(oo_.steady_state,1,size(oo_.endo_simul,2));
        else
            oo_.endo_simul = endo_simul;
        end
        t0 = tic;
        [flag,tmp] = bytecode('dynamic');
        TeaTime = toc(t0);
        ctime = ctime+TeaTime;
        old_weight = weight;
        info.convergence = ~flag;
        if verbose
            if ~info.convergence
                disp(['Iteration n° ' int2str(iter) ', weight is ' num2str(old_weight,8) ', Convergence problem!' ])
            else
                disp(['Iteration n° ' int2str(iter) ', weight is ' num2str(old_weight,8) ', Ok!' ])
            end
        end
        if ~info.convergence
            if iter==1
                disp('I am not able to simulate this model!')
                disp('There is something wrong with the initial condition of the homotopic')
                disp('approach...')
                error(' ')
            else
                if verbose
                    disp('I am reducing the step length!')
                end
                step_length=step_length/10;
                if 10*step_length<options_.dynatol
                    homotopy_1 = 0;
                    homotopy_2 = 0;
                    break
                end
                weight = old_weight+step_length;
            end
        else
            oo_.endo_simul = tmp;
            info.time = ctime;
            if abs(1-weight)<=1e-12;
                homotopy_2 = 1;
                break
            end
            weight = weight+step_length;
            step_length = initial_step_length;
        end
        if iter>max_iter
            info = NaN;
            return
        end
    end
    if weight<1 && homotopy_2
        oo_.exo_simul(2,:) = exxo_simul(2,:);
        t0 = tic; 
        [flag,tmp] = bytecode('dynamic');
        TeaTime = toc(t0);
        ctime = ctime+TeaTime;
        info.convergence = ~flag;
        info.time = ctime;
        if info.convergence
            oo_.endo_simul = tmp;
            homotopy_1 = 0;
        else
            if step_length>1e-12
                if verbose
                    disp('I am reducing step length!')
                end
                step_length=step_length/2;
            else
                weight = initial_weight;
                step_length = initial_step_length;
                info = NaN;
                homotopy_2 = 1;
                homotopy_1 = 0;
            end
        end
    end
end