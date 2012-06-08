//options_.usePartInfo=1;

var MC E EF R_KF QF CF IF YF LF PIEF WF RF R_K Q C I Y L PIE W R EE_A PIE_BAR EE_B EE_G EE_L EE_I KF K one BIGTHETA;    

varexo E_A E_B E_G E_L E_I ETA_R E_PIE_BAR ETA_Q ETA_P ETA_W  ;  
 
parameters 
xi_e lambda_w alpha czcap beta phi_i tau sig_c hab ccs cinvs phi_y gamma_w xi_w gamma_p xi_p sig_l r_dpi 
r_pie r_dy r_y rho rho_a rho_pb rho_b rho_g rho_l rho_i LMP  ;



alpha=.30;
beta=.99;
tau=0.025;
ccs=0.6;
cinvs=.22;  //% alpha*(tau+ctrend)/R_K   R_K=ctrend/beta-1+tau  
lambda_w = 0.5;
phi_i= 6.771;
sig_c=   1.353; 
hab=    0.573;    
xi_w=   0.737;
sig_l=    2.400;
xi_p=   0.908;
xi_e= 0.599;
gamma_w=    0.763;
gamma_p=    0.469;
czcap=    0.169;
phi_y=    1.408;
r_pie=     1.684;
r_dpi=    0.14;
rho=      0.961;
r_y=      0.099;
r_dy=     0.159;
rho_a=    0.823;
rho_b=    0.855;
rho_g=    0.949;
rho_l=   0.889;
rho_i=   0.927;
rho_pb=  0.924;
LMP = 0.0 ; //NEW.

model(linear, use_dll); 
          CF = (1/(1+hab))*(CF(1)+hab*CF(-1))-((1-hab)/((1+hab)*sig_c))*(RF-PIEF(1)-EE_B) ;
	      0 =  alpha*R_KF+(1-alpha)*WF -EE_A ;
          PIEF = 0*one;
	      IF = (1/(1+beta))* ((  IF(-1) + beta*(IF(1)))+(1/phi_i)*QF)+0*ETA_Q+EE_I ;
	      QF = -(RF-PIEF(1))+(1-beta*(1-tau))*((0+czcap)/czcap)*R_KF(1)+beta*(1-tau)*QF(1) +0*EE_I ;
          KF =  (1-tau)*KF(-1)+tau*IF(-1) ;
	      YF = (ccs*CF+cinvs*IF)+EE_G  ;

	      YF = 1*phi_y*( alpha*KF+alpha*(1/czcap)*R_KF+(1-alpha)*LF+EE_A ) ;
	      WF = (sig_c/(1-hab))*(CF-hab*CF(-1)) + sig_l*LF - EE_L ;
	      LF = R_KF*((1+czcap)/czcap)-WF+KF ;
          EF = EF(-1)+EF(1)-EF+(LF-EF)*((1-xi_e)*(1-xi_e*beta)/(xi_e));
         
	      C = (hab/(1+hab))*C(-1)+(1/(1+hab))*C(1)-((1-hab)/((1+hab)*sig_c))*(R-PIE(1)-EE_B) ;
	      I = (1/(1+beta))* ((  I(-1) + beta*(I(1)))+(1/phi_i)*Q )+1*ETA_Q+1*EE_I ;
	      Q = -(R-PIE(1))+(1-beta*(1-tau))*((0+czcap)/czcap)*R_K(1)+beta*(1-tau)*Q(1) +EE_I*0+0*ETA_Q ;
	      K =  (1-tau)*K(-1)+tau*I(-1) ;
	      Y = (ccs*C+cinvs*I)+ EE_G   ;
	      Y = phi_y*( alpha*K+alpha*(1/czcap)*R_K+(1-alpha)*L ) +phi_y*EE_A  ;
	      PIE = (1/(1+beta*gamma_p))*
	            ( 
	            (beta)*(PIE(1)) +(gamma_p)*(PIE(-1)) 
	            +((1-xi_p)*(1-beta*xi_p)/(xi_p))*(MC)
	            )  + ETA_P ; 
	            
	      MC = alpha*R_K+(1-alpha)*W -EE_A;
	      W =  (1/(1+beta))*(beta*W(+1)+W(-1))
                +(beta/(1+beta))*(PIE(+1))
                -((1+beta*gamma_w)/(1+beta))*(PIE)
                +(gamma_w/(1+beta))*(PIE(-1))
                -(1/(1+beta))*(((1-beta*xi_w)*(1-xi_w))/(((1+(((1+lambda_w)*sig_l)/(lambda_w))))*xi_w))*(W-sig_l*L-(sig_c/(1-hab))*(C-hab*C(-1))+EE_L)
                +ETA_W;
	      L = R_K*((1+czcap)/czcap)-W+K ;

//	      R = r_dpi*(PIE-PIE(-1))
//              +(1-rho)*(r_pie*(PIE(-1)-PIE_BAR)+r_y*(Y-YF))
//              +r_dy*(Y-YF-(Y(-1)-YF(-1)))
//              +rho*(R(-1)-PIE_BAR)
//              +PIE_BAR
//              +ETA_R;


	      R = 

r_dpi*(PIE-PIE(-1))

              +(1-rho)*(r_pie*(BIGTHETA)+r_y*(Y-YF))
              +r_dy*(Y-YF-(Y(-1)-YF(-1)))
              +rho*(R(-1)-PIE_BAR)
              +PIE_BAR
              +ETA_R;


          E = E(-1)+E(1)-E+(L-E)*((1-xi_e)*(1-xi_e*beta)/(xi_e));
          
          
          EE_A = (rho_a)*EE_A(-1)  + E_A;
	      PIE_BAR = rho_pb*PIE_BAR(-1)+ E_PIE_BAR ;
	      EE_B = rho_b*EE_B(-1) + E_B ;
	      EE_G = rho_g*EE_G(-1) + E_G ;
	      EE_L = rho_l*EE_L(-1) + E_L ;
	      EE_I = rho_i*EE_I(-1) + E_I ;
	      one = 0*one(-1) ;

		LMP*BIGTHETA(1) = BIGTHETA - (PIE(-1) - PIE_BAR) ; 

end; 

 
shocks;
var E_A; stderr 0.598;
var E_B; stderr 0.336;
var E_G; stderr 0.325;
var E_I; stderr 0.085;
var E_L; stderr 3.520;
var ETA_P; stderr 0.160;
var ETA_W; stderr 0.289;
var ETA_R; stderr 0.081;
var ETA_Q; stderr 0.604;
var E_PIE_BAR; stderr 0.017;
end;

//stoch_simul(irf=20) Y C PIE R W R_K L Q I K ;

// stoch_simul generates what kind of standard errors for the shocks ?

//steady;
//check;
//stoch_simul(periods=200,irf=20,simul_seed=3) Y C PIE MC R W R_K E L I ;

//datatomfile('ddd',[]);

// new syntax 

estimated_params;
// PARAM NAME, INITVAL, LB, UB, PRIOR_SHAPE, PRIOR_P1, PRIOR_P2, PRIOR_P3, PRIOR_P4, JSCALE
// PRIOR_SHAPE: BETA_PDF, GAMMA_PDF, NORMAL_PDF, INV_GAMMA_PDF
stderr E_A,0.543,0.01,4,INV_GAMMA_PDF,0.4,2;
stderr E_PIE_BAR,0.072,0.001,4,INV_GAMMA_PDF,0.02,10;
stderr E_B,0.2694,0.01,4,INV_GAMMA_PDF,0.2,2;
stderr E_G,0.3052,0.01,4,INV_GAMMA_PDF,0.3,2;
stderr E_L,1.4575,0.1,6,INV_GAMMA_PDF,1,2;
stderr E_I,0.1318,0.01,4,INV_GAMMA_PDF,0.1,2;
stderr ETA_R,0.1363,0.01,4,INV_GAMMA_PDF,0.1,2;
stderr ETA_Q,0.4842,0.01,4,INV_GAMMA_PDF,0.4,2;
stderr ETA_P,0.1731,0.01,4,INV_GAMMA_PDF,0.15,2;
stderr ETA_W,0.2462,0.1,4,INV_GAMMA_PDF,0.25,2;
rho_a,.9722,.1,.9999,BETA_PDF,0.85,0.1;
rho_pb,.85,.1,.999,BETA_PDF,0.85,0.1;
rho_b,.7647,.1,.99,BETA_PDF,0.85,0.1;
rho_g,.9502,.1,.9999,BETA_PDF,0.85,0.1;
rho_l,.9542,.1,.9999,BETA_PDF,0.85,0.1;
rho_i,.6705,.1,.99,BETA_PDF,0.85,0.1;
phi_i,5.2083,1,15,NORMAL_PDF,4,1.5;
sig_c,0.9817,0.25,3,NORMAL_PDF,1,0.375;
hab,0.5612,0.3,0.95,BETA_PDF,0.7,0.1;
xi_w,0.7661,0.3,0.9,BETA_PDF,0.75,0.05;
sig_l,1.7526,0.5,5,NORMAL_PDF,2,0.75;
xi_p,0.8684,0.3,0.95,BETA_PDF,0.75,0.05;
xi_e,0.5724,0.1,0.95,BETA_PDF,0.5,0.15;
gamma_w,0.6202,0.1,0.99,BETA_PDF,0.75,0.15;
gamma_p,0.6638,0.1,0.99,BETA_PDF,0.75,0.15;
czcap,0.2516,0.01,2,NORMAL_PDF,0.2,0.075;
phi_y,1.3011,1.001,2,NORMAL_PDF,1.45,0.125;
r_pie,1.4616,1.2,2,NORMAL_PDF,1.7,0.1;
r_dpi,0.1144,0.01,0.5,NORMAL_PDF,0.3,0.1;
rho,0.8865,0.5,0.99,BETA_PDF,0.8,0.10;
r_y,0.0571,0.01,0.2,NORMAL_PDF,0.125,0.05;
r_dy,0.2228,0.05,0.5,NORMAL_PDF,0.0625,0.05;
end;

varobs Y C I E PIE W R;

//estimation(datafile=rawdata_euromodel_1,presample=40, first_obs=1, nobs=118, lik_init=2, mode_compute=1,mh_replic=0);
estimation(datafile=rawdata_euromodel_1,presample=40, first_obs=1, nobs=118,mh_jscale=0.2,mh_replic=150000, mode_check); //


//stoch_simul(periods=200,irf=20,simul_seed=3) Y C PIE R W R_K L Q I K ;

