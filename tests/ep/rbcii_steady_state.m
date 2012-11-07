function [ys_, params, info] = rbcii_steadystate2(ys_, exo_, params)
     
    % Flag initialization (equal to zero if the deterministic steady state exists) 
    info = 0;
    
    % efficiency
    ys_(13)=0;
    
    % Efficiency
    ys_(12)=params(8);
    
    % Steady state ratios 
    Output_per_unit_of_Capital=((1/params(1)-1+params(6))/params(4))^(1/(1-params(5)));
    Consumption_per_unit_of_Capital=Output_per_unit_of_Capital-params(6);
    Labour_per_unit_of_Capital=(((Output_per_unit_of_Capital/ys_(12))^params(5)-params(4))/(1-params(4)))^(1/params(5));
    Output_per_unit_of_Labour=Output_per_unit_of_Capital/Labour_per_unit_of_Capital;
    Consumption_per_unit_of_Labour=Consumption_per_unit_of_Capital/Labour_per_unit_of_Capital;

    % Steady state share of capital revenues in total revenues (calibration check) 
    ShareOfCapital=params(4)/(params(4)+(1-params(4))*Labour_per_unit_of_Capital^params(5));

    % Steady state level of labour
    ys_(3)=1/(1+Consumption_per_unit_of_Labour/((1-params(4))*params(2)/(1-params(2))*Output_per_unit_of_Labour^(1-params(5))));
    
    % Steady state level of consumption
    ys_(4)=Consumption_per_unit_of_Labour*ys_(3);
    
    % Steady state level of physical capital stock
    ys_(1)=ys_(3)/Labour_per_unit_of_Capital;
    
    % Steady state level of output
    ys_(2)=Output_per_unit_of_Capital*ys_(1);
    
    % Steady state level of investment
    ys_(5)=params(6)*ys_(1);
    
    % Steady state level of the expected term appearing in the Euler equation
    ys_(14)=params(1)*(ys_(4)^params(2)*(1-ys_(3))^(1-params(2)))^(1-params(3))/ys_(4)*(1+params(4)*(ys_(2)/ys_(1))^(1-params(5))-params(6));

    % Steady state level of output in the unconstrained regime (positive investment)
    ys_(6)=ys_(2);

    % Steady state level of labour in the unconstrained regime
    ys_(7)=ys_(3);
    
    % Steady state level of consumption in the unconstrained regime 
    ys_(8)=ys_(4);
        
    % Steady state level of labour in the constrained regime (noinvestment)
    [lss,info] = l_solver(ys_(3),params(4),params(5),params(2),params(8),ys_(1),100);
    if info, return, end
    ys_(10) = lss;

    % Steady state level of consumption in the constrained regime
    ys_(11)=params(8)*(params(4)*ys_(1)^params(5)+(1-params(4))*ys_(10)^params(5))^(1/params(5));
    
    % Steady state level of output in the constrained regime
    ys_(9)=ys_(11);

end




function r = p0(labour,alpha,psi,theta,effstar,kstar)
    r = labour * ( alpha*kstar^psi/labour^psi + 1-alpha + theta*(1-alpha)/(1-theta)/effstar^psi ) - theta*(1-alpha)/(1-theta)/effstar^psi;
end
    
function d = p1(labour,alpha,psi,theta,effstar,kstar)
    d = alpha*(1-psi)*kstar^psi/labour^psi + 1-alpha + theta*(1-alpha)/(1-alpha)/effstar^psi;
end

function [labour,info] = l_solver(labour,alpha,psi,theta,effstar,kstar,maxiter)
    iteration = 1; info = 0;
    r = p0(labour,alpha,psi,theta,effstar,kstar);
    condition = abs(r);
    while condition
        if iteration==maxiter
            info = 1;
            break
        end
        d = p1(labour,alpha,psi,theta,effstar,kstar);
        labour = labour - r/d;
        r = p0(labour,alpha,psi,theta,effstar,kstar);
        condition = abs(r)>1e-9;
        iteration = iteration + 1; 
    end
end