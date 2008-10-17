function compare_kalman_routines(experience)
%% This function compare the kalman filter routines.
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


%% BUILD DATA SET (zero mean):
Y = randn(pp,gend);
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

[data_index,number_of_observations,no_more_missing_observations] = describe_missing_data(data,gend,nvarobs);


%% SET THE STATE SPACE MODEL:

% I randomly choose the mm eigenvalues of the transition matrix.
TransitionEigenvalues  = rand(mm,1)*2-1;
% I randomly choose the mm eigenvectors of the transition matrix
TransitionEigenvectors = rand(mm,A);
% I build the transition matrix
T = TransitionEigenvectors*diag(TransitionEigenvectors)*inv(TransitionEigenvectors);
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
% Set the selec tion vector (mf) 
MF = transpose(randperm(mm));
mf = tmp(1:pp);



if PeriodsWithMissingObservations==0
    instant0 = clock;  
    LIK1 = kalman_filter(T,R,Q,H,P,Y,start,mf,kalman_tol,riccati_tol);
    T1 = etime(clock, instant0);
    disp(['kalman_filter = ' num2str(T1)])
    instant0 = clock;  
    LIK2 = missing_observations_kalman_filter(T,R,Q,H,P,Y,start,mf,kalman_tol,riccati_tol,data_index,number_of_observations,no_more_missing_observations);
    T2 = etime(clock, instant0);
    disp(['missing_observations_kalman_filter = ' num2str(T2)])
    instant0 = clock;
    LIK3 = univariate_kalman_filter(T,R,Q,H,P,Y,start,mf,kalman_tol,riccati_tol,data_index,number_of_observations,no_more_missing_observations);
    T3 = etime(clock, instant0);
    disp(['univariate_kalman_filter = ' num2str(T2)])
    if T1==T2
        disp('missing data version is Ok')
    else
        LIK1-LIK2
    end
    if T1==T3
        disp('univariate version is Ok')
    else
        LIK1-LIK3
    end
else
    instant0 = clock;  
    LIK2 = missing_observations_kalman_filter(T,R,Q,H,P,Y,start,mf,kalman_tol,riccati_tol,data_index,number_of_observations,no_more_missing_observations);
    T2 = etime(clock, instant0);
    disp(['missing_observations_kalman_filter = ' num2str(T2)])
    instant0 = clock;
    LIK3 = univariate_kalman_filter(T,R,Q,H,P,Y,start,mf,kalman_tol,riccati_tol,data_index,number_of_observations,no_more_missing_observations);
    T3 = etime(clock, instant0);
    disp(['univariate_kalman_filter = ' num2str(T2)])
    if T1==T3
        disp('univariate version is Ok')
    else
        LIK1-LIK3
    end
else
    