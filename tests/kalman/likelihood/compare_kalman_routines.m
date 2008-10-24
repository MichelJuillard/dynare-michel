function compare_kalman_routines(experience)
%% This function compares the kalman filter routines.
%%
%%     
%% stephane [DOT] adjemian [AT] ens [DOT] fr

pp = experience.Number0fObservedVariables;
mm = experience.SizeOfTheStateVector;
rr = experience.NumberOfStructuralShocks;
measurement_error_flag = experience.MeasurementErrors;
gend = experience.NumberOfPeriods;
PercentageOfMissingObservations = experience.PercentageOfMissingObservations;
PeriodsWithMissingObservations = experience.PeriodsWithMissingObservations;



%% SET VARIOUS PARAMETERS 
kalman_tol = 1e-12;
riccati_tol =1e-9;
start = 1;



%% SET THE STATE SPACE MODEL:

% I randomly choose the mm eigenvalues of the transition matrix.
TransitionEigenvalues  = rand(mm,1)*2-1;
% I randomly choose the mm eigenvectors of the transition matrix
tmp = rand(mm,mm*100);
TransitionEigenvectors = tmp*tmp'/(mm*100);
TransitionEigenvectors = rand(mm,mm);
% I build the transition matrix
T = TransitionEigenvectors*diag(TransitionEigenvalues)*inv(TransitionEigenvectors);
% I randomly choose matrix R
R = randn(mm,rr);
% I randomly choose the covariance matrix of the structurtal innovations
E = randn(rr,20*rr);
Q = E*transpose(E)/(20*rr);
% If needed I randomly choose the covariance matrix of teh measurement errors 
if measurement_error_flag == 0
    H = zeros(pp,1);
elseif measurement_error_flag == 1
    H = rand(pp,1);
elseif measurement_error_flag == 2
    E = randn(pp,20*pp);
    H = E*transpose(E)/(20*pp);
else
    disp('compare_kalman_routines:: Unknown option!')
end
% Set the selection vector (mf) 
MF = transpose(randperm(mm));
mf = MF(1:pp);

P = lyapunov_symm(T,R*Q*R',1.000001);


%% BUILD DATA SET (zero mean):

if measurement_error_flag == 0
    Y = simul_state_space_model(T,R,Q,mf,gend);
elseif measurement_error_flag == 1
    H = rand(pp,1);
    Y = simul_state_space_model(T,R,Q,mf,gend,diag(H));
elseif measurement_error_flag == 2
    E = randn(pp,20*pp);
    H = E*transpose(E)/(20*pp);
    Y = simul_state_space_model(T,R,Q,mf,gend,H);
else
    disp('compare_kalman_routines:: Unknown option!')
end

if PeriodsWithMissingObservations
    data_1 = Y(:,1:PeriodsWithMissingObservations)';
    if PeriodsWithMissingObservations<gend
        data_2 = Y(:,PeriodsWithMissingObservations+1:end);
    else
        data_2 = [];
    end
    tmp = randperm(1:PeriodsWithMissingObservations*pp);
    tmp = tmp(1:round(PercentageOfMissingObservations*length(tmp)));
    data_1(tmp) = NaN;
    Y = [ transpose(data_1) ; data_2 ];
end

[data_index,number_of_observations,no_more_missing_observations] = describe_missing_data(Y,gend,pp);

if PeriodsWithMissingObservations==0
    
    % kalman_filter.m
    if measurement_error_flag==0
        HH = 0;
    elseif measurement_error_flag==1
        HH = diag(H);
    elseif measurement_error_flag==2
        HH = H;
    end
    instant0 = clock;
    [LIK1,lik1] = kalman_filter(T,R,Q,HH,P,Y,start,mf,kalman_tol,riccati_tol);
    T1 = etime(clock, instant0);
    disp(['kalman_filter = ' num2str(T1)])
    
    % missing_observations_kalman_filter.m
    if measurement_error_flag==0
        HH = zeros(pp,pp);
    elseif measurement_error_flag==1
        HH = diag(H);
    elseif measurement_error_flag==2
        HH = H;
    end
    instant0 = clock;  
    [LIK2,lik2] = missing_observations_kalman_filter(T,R,Q,HH,P,Y,start,mf,kalman_tol,riccati_tol,data_index,number_of_observations,no_more_missing_observations);
    T2 = etime(clock, instant0);
    disp(['missing_observations_kalman_filter = ' num2str(T2)])
    
    
    % univariate_kalman_filter.m
    if measurement_error_flag==0
        HH = zeros(pp,pp);
    elseif measurement_error_flag==1
        HH = H;
    elseif measurement_error_flag==2
        HH = H;
    end
    if measurement_error_flag<2
        instant0 = clock;
        [LIK3,lik3] = univariate_kalman_filter(T,R,Q,HH,P,Y,start,mf,kalman_tol,riccati_tol,data_index,number_of_observations,no_more_missing_observations);
        T3 = etime(clock, instant0);
    else
        instant0 = clock;
        [LIK3,lik3] = univariate_kalman_filter_corr(T,R,Q,HH,P,Y,start,mf,kalman_tol,riccati_tol,data_index,number_of_observations,no_more_missing_observations);
        T3 = etime(clock, instant0);
    end
    
    disp(['univariate_kalman_filter = ' num2str(T3)])
    
    disp(' ')
    disp(' ')
    
    if abs(LIK2/LIK1-1)<1e-9
        disp('missing data version is Ok')
    else
        disp('missing data version is wrong')
        disp(['percentage dev. = ' num2str((LIK2/LIK1-1)*100)])
    end
    if abs(LIK3/LIK1-1)<1e-9
        disp('univariate version is Ok')
    else
        disp('univariate version is wrong')
        disp(['percentage dev. = ' num2str((LIK3/LIK1-1)*100)])
        %[lik1,lik3,lik1-lik3]        
        %[lik1-lik3]
    end
else
    % missing_observations_kalman_filter.m
    if measurement_error_flag==0
        HH = zeros(pp,pp);
    elseif measurement_error_flag==1
        HH = diag(H);
    elseif measurement_error_flag==2
        HH = H;
    end    
    instant0 = clock;  
    LIK2 = missing_observations_kalman_filter(T,R,Q,HH,P,Y,start,mf,kalman_tol,riccati_tol,data_index,number_of_observations,no_more_missing_observations);
    T2 = etime(clock, instant0);
    disp(['missing_observations_kalman_filter = ' num2str(T2)])
    
    % univariate_kalman_filter.m
    if measurement_error_flag==0
        HH = zeros(pp,pp);
    elseif measurement_error_flag==1
        HH = H;
    elseif measurement_error_flag==2
        error('The univariate approach cannot handle correlated measurement errors')
    end
    instant0 = clock;
    LIK3 = univariate_kalman_filter(T,R,Q,HH,P,Y,start,mf,kalman_tol,riccati_tol,data_index,number_of_observations,no_more_missing_observations);
    T3 = etime(clock, instant0);
    disp(['univariate_kalman_filter = ' num2str(T2)])
    
    disp(' ')
    disp(' ')
    
    if abs(LIK3/LIK2-1)<1e-9
        disp('univariate version is Ok')
    else
        disp('univariate version is wrong')
        disp(['percentage dev. = ' num2str((LIK3/LIK2-1)*100)])        
    end
end
    