load test
% $$$ Y = Y(1:2,:);
% $$$ mf = mf(1:2);
% $$$ H=H(1:2,1:2);
% $$$ pp = pp-1;
% $$$ trend =trend(1:2,:);
Pinf1(1,1) = 1;
Pstar1(1,1) = 0;
Pstar1(4,1) = 0;
Pstar1(1,4) = 0;
[alphahat1,epsilonhat1,etahat1,a11, aK1] = DiffuseKalmanSmootherH1(T,R,Q,H, ...
						  Pinf1,Pstar1,Y,trend,pp,mm,smpl,mf);
[alphahat2,epsilonhat2,etahat2,a12, aK2] = DiffuseKalmanSmootherH3(T,R,Q,H, ...
						  Pinf1,Pstar1,Y,trend, ...
						  pp,mm,smpl,mf);
max(max(abs(alphahat1-alphahat2)))
max(max(abs(epsilonhat1-epsilonhat2)))
max(max(abs(etahat1-etahat2)))
max(max(abs(a11-a12)))
max(max(abs(aK1-aK2)))

return
[alphahat1,etahat1,a11, aK1] = DiffuseKalmanSmoother1(T,R,Q, ...
						  Pinf1,Pstar1,Y,trend,pp,mm,smpl,mf);
[alphahat2,etahat2,a12, aK2] = DiffuseKalmanSmoother3(T,R,Q, ...
						  Pinf1,Pstar1,Y,trend, ...
						  pp,mm,smpl,mf);


max(max(abs(alphahat1-alphahat2)))
max(max(abs(etahat1-etahat2)))
max(max(abs(a11-a12)))
%max(max(abs(aK1-aK2)))


H = zeros(size(H));
[alphahat1,etahat1,a11, aK1] = DiffuseKalmanSmoother1(T,R,Q, ...
						  Pinf1,Pstar1,Y,trend,pp,mm,smpl,mf);
[alphahat2,epsilonhat2,etahat2,a12, aK2] = DiffuseKalmanSmootherH1(T,R,Q,H, ...
						  Pinf1,Pstar1,Y,trend, ...
						  pp,mm,smpl,mf);
max(max(abs(alphahat1-alphahat2)))
max(max(abs(etahat1-etahat2)))
max(max(abs(a11-a12)))
%max(max(abs(aK1-aK2)))


[alphahat1,etahat1,a11, aK1] = DiffuseKalmanSmoother3(T,R,Q, ...
						  Pinf1,Pstar1,Y,trend,pp,mm,smpl,mf);
[alphahat2,epsilonhat2,etahat2,a12, aK2] = DiffuseKalmanSmootherH3(T,R,Q, H, ...
						  Pinf1,Pstar1,Y,trend,pp,mm,smpl,mf);

max(max(abs(alphahat1-alphahat2)))
max(max(abs(etahat1-etahat2)))
max(max(abs(a11-a12)))
%max(max(abs(aK1-aK2)))
