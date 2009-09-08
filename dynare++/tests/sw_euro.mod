var MC EH EF R_KF QF CF IF YF LF PIEF WF RF R_KH QH CH IH YH LH PIEH WH RH EE_A PIE_BAR EE_B EE_G EE_L EE_I KF KH ONE;    

varexo E_A E_B E_G E_L E_I ETA_R E_PIE_BAR ETA_Q ETA_P ETA_W  ;  
 
parameters xi_e lambda_w alpha czcap beta phi_i tau sig_c hab ccs cinvs phi_y gamma_w xi_w gamma_p xi_p sig_l r_dpi r_pie r_dy r_y rho rho_a rho_pb rho_b rho_g rho_l rho_i  ;
alpha=.30;
beta=0.99;
tau=0.025;
ccs=0.6;
cinvs=.22;
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


model; 
          CF = (1/(1+hab))*(CF(1)+hab*CF(-1))-((1-hab)/((1+hab)*sig_c))*(RF-PIEF(1)-EE_B) ;
	      0 =  alpha*R_KF+(1-alpha)*WF -EE_A ;
          PIEF = 0*ONE;
	      IF = (1/(1+beta))* ((  IF(-1) + beta*(IF(1)))+(1/phi_i)*QF)+0*ETA_Q+EE_I ;
	      QF = -(RF-PIEF(1))+(1-beta*(1-tau))*((1+czcap)/czcap)*R_KF(1)+beta*(1-tau)*QF(1) +0*EE_I ;
          KF =  (1-tau)*KF(-1)+tau*IF(-1) ;
	      YF = (ccs*CF+cinvs*IF)+EE_G  ;
	      YF = 1*phi_y*( alpha*KF+alpha*(1/czcap)*R_KF+(1-alpha)*LF+EE_A ) ;
	      WF = (sig_c/(1-hab))*(CF-hab*CF(-1)) + sig_l*LF - EE_L ;
	      LF = R_KF*((1+czcap)/czcap)-WF+KF ;
          EF = EF(-1)+EF(1)-EF+(LF-EF)*((1-xi_e)*(1-xi_e*beta)/(xi_e));
         
	      CH = (hab/(1+hab))*CH(-1)+(1/(1+hab))*CH(1)-((1-hab)/((1+hab)*sig_c))*(RH-PIEH(1)-EE_B) ;
	      IH = (1/(1+beta))* ((  IH(-1) + beta*(IH(1)))+(1/phi_i)*QH )+1*ETA_Q+1*EE_I ;
	      QH = -(RH-PIEH(1))+(1-beta*(1-tau))*((1+czcap)/czcap)*R_KH(1)+beta*(1-tau)*QH(1) +EE_I*0+0*ETA_Q ;
	      KH =  (1-tau)*KH(-1)+tau*IH(-1) ;
	      YH = (ccs*CH+cinvs*IH)+ EE_G   ;
	      YH = phi_y*( alpha*KH+alpha*(1/czcap)*R_KH+(1-alpha)*LH ) +phi_y*EE_A  ;
	      PIEH = (1/(1+beta*gamma_p))*
	            ( 
	            (beta)*(PIEH(1)) +(gamma_p)*(PIEH(-1)) 
	            +((1-xi_p)*(1-beta*xi_p)/(xi_p))*(MC)
	            )  + ETA_P ; 
	            
	      MC = alpha*R_KH+(1-alpha)*WH -EE_A;
	      WH =  (1/(1+beta))*(beta*WH(+1)+WH(-1))
                +(beta/(1+beta))*(PIEH(+1))
                -((1+beta*gamma_w)/(1+beta))*(PIEH)
                +(gamma_w/(1+beta))*(PIEH(-1))
                -(1/(1+beta))*(((1-beta*xi_w)*(1-xi_w))/(((1+(((1+lambda_w)*sig_l)/(lambda_w))))*xi_w))*(WH-sig_l*LH-(sig_c/(1-hab))*(CH-hab*CH(-1))+EE_L)
                +ETA_W;
	      LH = R_KH*((1+czcap)/czcap)-WH+KH ;
	      RH = r_dpi*(PIEH-PIEH(-1))
              +(1-rho)*(r_pie*(PIEH(-1)-PIE_BAR)+r_y*(YH-YF))
              +r_dy*(YH-YF-(YH(-1)-YF(-1)))
              +rho*(RH(-1)-PIE_BAR)
              +PIE_BAR
              +ETA_R;
          EH = EH(-1)+EH(1)-EH+(LH-EH)*((1-xi_e)*(1-xi_e*beta)/(xi_e));
          
          
          EE_A = (rho_a)*EE_A(-1)  + E_A;
	      PIE_BAR = rho_pb*PIE_BAR(-1)+ E_PIE_BAR ;
	      EE_B = rho_b*EE_B(-1) + E_B ;
	      EE_G = rho_g*EE_G(-1) + E_G ;
	      EE_L = rho_l*EE_L(-1) + E_L ;
	      EE_I = rho_i*EE_I(-1) + E_I ;
	      ONE = 0*ONE(-1) ;
end; 

vcov = [0.357604 0 0 0 0 0 0 0 0 0;
        0 0.112896 0 0 0 0 0 0 0 0;
        0 0 0.105625 0 0 0 0 0 0 0;
        0 0 0 12.39040 0 0 0 0 0 0;
        0 0 0 0 0.722500 0 0 0 0 0;
        0 0 0 0 0 0.656100 0 0 0 0;
        0 0 0 0 0 0 0.000289 0 0 0;
        0 0 0 0 0 0 0 0.364816 0 0;
        0 0 0 0 0 0 0 0 0.025600 0;
        0 0 0 0 0 0 0 0 0 0.083521];

order = 1;
