var A, Disp, G, Int, L,
		   LStar, pi, Welf, WelfStar,  x0,
		   Y, YGap, YStar, z1, z2, Cbar, Cequiv; 
varexo eps1 eps2 eps3;

parameters alpha  beta  gamma  rhoa  rhog  rho phi  chi  chi0  theta  xi  
ABar  GBar  KBar  ZBar  piBar Istar;
alpha = 0.3;
  beta = 0.99;
  gamma = 15;
  rhoa = 0.8;
  rhog = 0.7;
  phi = 1.5;
  chi = 1.5;
  chi0 = 1;
  theta = 0.333333333333;
  xi = 0.75;
  ABar = 4.0266;
  GBar = 0.3163;
  KBar = 9.489;
  ZBar = .03;
  piBar = 1;
  rho=.8;
  Istar=1.01010101010101;


model;
z1 - ((Y-G)^(1-phi) + beta *xi *piBar *pi(+1)^(1/theta) *z1(+1));
z2 - (Y *chi0 *(1-L-ZBar)^(-chi) / ((1-alpha) *A *KBar^alpha
	*L^(-alpha)) + beta *xi *pi(+1)^((1+theta)/theta) *z2(+1));
x0 - (1+theta)*z2 /z1;
pi^(-1/theta) - ((1-xi) *(x0*pi)^(-1/theta) + xi *piBar^(-1/theta));
Y - (Disp^(-1) *A *KBar^alpha *L^(1-alpha));
Disp - ((1-xi) *x0^(-(1+theta)/theta)
	+ xi *(piBar/pi)^(-(1+theta)/theta) *Disp(-1));
log(A/ABar) - (rhoa *log(A(-1)/ABar) + eps1);
log(G/GBar) - (rhog *log(G(-1)/GBar) + eps2);
(Y-G)^(-phi) - (beta *(Int/pi(+1)) *(Y(+1)-G(+1))^(-phi));
Welf - ((Y-G)^(1-phi) /(1-phi)
	+ chi0*(1-L-ZBar)^(1-chi) /(1-chi) + beta *Welf(+1));
Cequiv = (((1-beta)*Welf-chi0*(1-LStar-ZBar)^(1-chi) /(1-chi))*(1-phi))^(1/(1-phi));
(1-alpha) *A *KBar^alpha *LStar^(-alpha)
	- (1+theta) *YStar *(YStar-G)^(phi-1) *chi0
	*(1-LStar-ZBar)^(-chi);
YStar - A *KBar^alpha *LStar^(1-alpha);
YGap - (log(Y/YStar));
WelfStar - ((YStar-G)^(1-phi) /(1-phi)
	+ chi0*(1-LStar-ZBar)^(1-chi) /(1-chi) + beta *WelfStar(+1));
Int = (Int(-1)^rho)*((Istar*(pi/piBar)^gamma)^(1-rho))*exp(eps3);
Cbar=(1/100)*((1-phi)*((1-beta)*WelfStar-chi0*(1-LStar-ZBar)^(1-chi)/(1-chi)))^(1/(1-phi));
end;

initval;
A=            4.022;	       
Disp=	      1;	       
G=	      0.3;	       
Int=	      1.0101;	       
L=	      0.22;	       
LStar=	       0.22;	       
pi=	      1;	       
Welf=	      -359;	       
WelfStar=     -359;	 
x0=	      1;	       
Y=	       2.8;	       
YGap=	      0;	       
YStar=	      2.8;	       
z1=	      2.5;	       
z2=	      1.8;             
Cbar= 0.024;
Cequiv = 0.024;
end;

vcov = [0.001 0 0 ; 0 0.001 0; 0 0 0.001];

order=4;


