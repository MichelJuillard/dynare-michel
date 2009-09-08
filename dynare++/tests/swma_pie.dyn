// this model has sticky wages and adjustment costs in
// investment, consumer goods sector is perfectly competitive, thus MC=1
// with money and transaction costs based on money velocity
// and it has a financial accelerator
// wage is indexed to past consumer price inflation

// LAMBDA Lagrange multiplier on household's budget constraint (divided by price level)
// PIE inflation of CPI
// PIETILDE to what inflation new wage setters index (here PIE(-1) but could be PIEW(-1))
// INT nominal interest rate
// C real consumption
// I real investment
// K real capital
// R real rental rate of capital
// W real wage
// L labour
// Y real output
// PIEW nominal wage inflation
// VW wage front loading term for newly set wages
// BBD, BBE, BBF, BBG terms in nominator and denominator in wage FOC
// G government
// SL process for labor shock
// SC process for consumption shock
// SY process for technology shock
// RM real money balances hold
// Q real price of capital
// Q_M1 lagged Q
// RK nominal return of capital for enterpreneurs
// OMEGABAR threshold value for idiosyncratic shock
// N real net worth of borrowers
// WF lifetime utility

var LAMBDA PIE PIETILDE INT C I K R W L Y PIEW VW BBD BBE BBF BBG G SL SC SY RM
    Q Q_M1 RK OMEGABAR N ACAL ACALPRIME BCAL BCALPRIME WF;

varexo E_C E_L E_Y E_GOV E_INT;

parameters dep beta gamma eta biga alpha sigmaw phiw deltaw sg pietar h psi nu osigma mu tc1 tc2 ksi1 ksi2 c_weight rho_g rho_l rho_c rho_y;
dep = 0.025;
beta = 0.99;
gamma = 1;
eta = 2;
alpha = 0.30;
biga = alpha^(-alpha)*(1-alpha)^(alpha-1);
sigmaw = 11;
phiw = 2;
deltaw = 0.75;
sg = 0.18;
pietar = 1.03^0.25;
h = 0.8;
// investment adjustment costs
psi = 12;
// enterpreneur saving rate
nu = 0.94;
// stderr of enterpreneur's idiosyncratic shocks
osigma = 0.5;
// monitoring cost for lender
mu = 0.2;
// consumption transaction costs
tc1 = 0.05;
tc2 = 0.5;
// Taylor rule
ksi1 = 0.106;
ksi2 = 3;
rho_g = 0.90;
rho_l = 0.90;
rho_c = 0.90;
rho_y = 0.90;
// weight of consumption utility 
c_weight = 1;

model;
// capital accumulation
K = (1 - dep - psi/2*(I(-1)/K(-1)-dep)^2)*K(-1) + I(-1);
// FOC bonds
LAMBDA = beta*INT*LAMBDA(+1)/PIE(+1);
// FOC consumption (right hand side is equal to LAMBDA*(1+TC+TCPRIME*C/RM))
SC*c_weight*(C-h*C(-1))^(-eta) = LAMBDA*(1+2*tc1*C/RM-2*sqrt(tc1*tc2));
// FOC money (right hand side is equal to 1 - TCPRIME*C*C/RM/RM)
beta*LAMBDA(+1)/LAMBDA/PIE(+1) = 1 - tc1*C*C/RM/RM + tc2;
// FOC investment removed
// FOC capital(+1) removed
// real price of capital
Q = (1-psi*(I/K-dep))^(-1);
// nominal return on capital
RK = PIE*(R + Q*(1 - dep + psi*(I/K-dep)*I/K -psi/2*(I/K-dep)^2))/Q(-1);
// FOC in optimal contract for K(+1)
RK(+1)*(BCAL(+1)*ACALPRIME(+1)/BCALPRIME(+1)-ACAL(+1)) = INT(+1)*ACALPRIME(+1)/BCALPRIME(+1);
// Participation constraint
//RK(+1)*BCAL(+1) = INT(+1)*(1-N(+1)*PIE(+1)/Q/K(+1));
RK*BCAL = INT*(1-N*PIE/Q(-1)/K);
// evolution of net worth (real)
N*PIE*PIE(-1) = nu*(ACAL(-1)+BCAL(-1))*RK(-1)*Q_M1(-1)*K(-1) - nu*INT(-1)*(Q_M1(-1)*K(-1)-N(-1)*PIE);
// marginal cost is 1
1 = biga*(W/SY)^(1-alpha)*R^alpha;
// labor attaining minimal MC
L = (1-alpha)/W*Y;
// capital attaining minimal MC
K = alpha/R*Y;
// FOC for newly set wages
W*VW = sigmaw/(sigmaw-1)*(BBD*VW^(-sigmaw*gamma) + phiw*BBE*VW^(-sigmaw) - phiw*BBF)/BBG;
// definition of BBD
BBD = SL*L^(1+gamma) + deltaw*beta*(PIETILDE(+1)/PIEW(+1))^(-sigmaw*(1+gamma))*BBD(+1);
// definition of BBE
BBE = LAMBDA*L*W + deltaw*beta*(PIETILDE(+1)/PIEW(+1))^(-2*sigmaw)*BBE(+1);
// definition of BBF
BBF = LAMBDA*L*W + deltaw*beta*(PIETILDE(+1)/PIEW(+1))^(-sigmaw)*BBF(+1);
// definition of BBG
BBG = LAMBDA*L + deltaw*beta*(PIETILDE(+1)/PIEW(+1))^(-sigmaw)*PIETILDE(+1)/PIE(+1)*BBG(+1);
// price index
1 = (1-deltaw)*VW^(1-sigmaw) + deltaw*(PIETILDE/PIEW)^(1-sigmaw);
// definition of ACAL
ACAL = 0.5*erfc((log(OMEGABAR) - 0.5*osigma^2)/osigma/sqrt(2.0)) - OMEGABAR/2*erfc((log(OMEGABAR) + 0.5*osigma^2)/osigma/sqrt(2.0));
// definition of BCAL
BCAL = OMEGABAR/2*erfc((log(OMEGABAR) + 0.5*osigma^2)/osigma/sqrt(2.0)) + (1-mu)/2*(1+erf((log(OMEGABAR) - 0.5*osigma^2)/osigma/sqrt(2.0)));
// definition of ACALPRIME
ACALPRIME = -0.5*erfc((log(OMEGABAR) + 0.5*osigma^2)/osigma/sqrt(2.0));
// definition of BCALPRIME
BCALPRIME = -ACALPRIME - mu/osigma/2.506628274631*exp(-((log(OMEGABAR) + 0.5*osigma)^2)/2/osigma/osigma);
// identity for PIEW
PIEW = PIE*W/W(-1);
// welfare identity
WF = SC*c_weight*(C-h*C(-1))^(1-eta)/(1-eta) - SL*L^(1+gamma)/(1+gamma) + beta*WF(+1);
// interest rate rule
INT = INT(-1)^ksi1*((PIE/beta)*(PIE/pietar)^ksi2)^(1-ksi1)*exp(E_INT);
// aggregate constraint
Y = C + I + G + (1-ACAL-BCAL)*RK*Q(-1)*K;
//Y = C + I + G;
// process for government
G/Y = (G(-1)/Y(-1))^rho_g*sg^(1-rho_g)*exp(E_GOV/sg);
// to what do they index (pietar, past inflation, past indexed inflation)
PIETILDE = PIE(-1);
//PIETILDE = pietar;
// exo processes
SL = SL(-1)^rho_l*exp(E_L);
SC = SC(-1)^rho_c*exp(E_C);
SY = SY(-1)^rho_y*exp(E_Y);
// lagged Q
Q_M1 = Q(-1);
end;

initval;
RM = 0.1;
INT = pietar/beta;
PIE = pietar;
PIEW = pietar;
PIETILDE = pietar;
//R = dep/beta;
R = 0.1;
W = (1/biga/(R)^alpha)^(1/(1-alpha));
LAMBDA = ((1-dep*alpha/R-sg)*(1-h)*c_weight/(1-alpha)*W^(1/gamma+1)*((sigmaw-1)/sigmaw)^(1/gamma))^(-1/(1/eta+1/gamma));
L = (W*LAMBDA*(sigmaw-1)/sigmaw)^(1/gamma);
Y = W*L/(1-alpha);
K = alpha/R*Y;
I = dep*K;
G = sg*Y;
VW = 1;
BBD = L^(1+gamma)/(1-deltaw*beta);
BBE = LAMBDA*L*W/(1-deltaw*beta);
BBF = LAMBDA*L*W/(1-deltaw*beta);
BBG = LAMBDA*L/(1-deltaw*beta);
Q = 1;
Q_M1 = Q;
RK = 1/Q*PIE*(R+(1-dep)*Q);
OMEGABAR = 0.5;
ACAL = 0.5*erfc((log(OMEGABAR) - 0.5*osigma^2)/osigma/sqrt(2.0)) - OMEGABAR/2*erfc((log(OMEGABAR) + 0.5*osigma^2)/osigma/sqrt(2.0));
BCAL = OMEGABAR/2*erfc((log(OMEGABAR) + 0.5*osigma^2)/osigma/sqrt(2.0)) + (1-mu)/2*(1+erf((log(OMEGABAR) - 0.5*osigma^2)/osigma/sqrt(2.0)));
ACALPRIME = -0.5*erfc((log(OMEGABAR) + 0.5*osigma^2)/osigma/sqrt(2.0));
BCALPRIME = -ACALPRIME - mu/osigma/2.506628274631*exp(-((log(OMEGABAR) + 0.5*osigma)^2)/2/osigma/osigma);
N = (nu*(ACAL+BCAL)*RK*Q*K-nu*INT*Q*K)/(PIE*PIE-nu*INT*PIE);
C = Y - I - G - (1-ACAL-BCAL)*RK*Q*K;
SL = 1;
SC = 1;
SY = 1;
WF = 1/(1-beta)*(SC*c_weight*((1-h)*C)^(1-eta)/(1-eta) - SL*L^(1+gamma)/(1+gamma));
end;

vcov = [
0.0001 0 0 0 0;
0 0.0001 0 0 0;
0 0 0.0001 0 0;
0 0 0 0.0001 0;
0 0 0 0 0.0001
];

order = 4;

