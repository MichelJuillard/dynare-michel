n = 10;
p = 3;
q = 4;

% With measurment errors.
for i=1:100
    tmp = randn(q,q+100);
    Q = tmp*transpose(tmp)/1000;
    R = randn(n,q);
    QQ = R*Q*transpose(R);
    tmp = randn(p,p+100);
    H = tmp*transpose(tmp)/2000;
    eig_val_T = rand(n)*2-1;
    eig_vec_T = randn(n,n);
    T = eig_vec_T*eig_val_T*inv(eig_vec_T);
    Z = zeros(p,n);
    Z(1,1) = 1;
    Z(2,3) = 1;
    Z(3,6) = 1;
    [err P] = kalman_steady_state(transpose(T),QQ,transpose(Z),H);
    mexErrCheck('kalman_steady_state',err);
end

% Without measurment errors.
for i=1:100
    tmp = randn(q,q+100);
    Q = tmp*transpose(tmp)/1000;
    R = randn(n,q);
    QQ = R*Q*transpose(R);
    H = zeros(p,p);
    eig_val_T = rand(n)*2-1;
    eig_vec_T = randn(n,n);
    T = eig_vec_T*eig_val_T*inv(eig_vec_T);
    Z = zeros(p,n);
    Z(1,1) = 1;
    Z(2,3) = 1;
    Z(3,6) = 1;
    [err P] = kalman_steady_state(transpose(T),QQ,transpose(Z),H);
    mexErrCheck('kalman_steady_state',err);
end
