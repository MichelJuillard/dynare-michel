#include "dynareR.cpp"

int main(void) {
	const char *parameters[] = {"beta","gamma","rho","alpha","delta"};
	const char *varendo[] = {"k","c","a"};
	const char *varexo[] = {"eps"};
	const int numpar = 5;
	const int numendo = 3;
	const int numexo = 1;
	const int ord = 2;
	const int numsteps = 0;
	const double parval[] = {.99,2,.9,.3,.025};
	const double vcov[] = {0.001};
	const double initval[] = {0.066, 0.43, 0.01};

	int e;
	double tensorbuffer[100];
	int num_state;
	int ordering_state[] = {0,0,0};
	int ordering_endo[] = {0,0,0};
	int ordering_exo[] = {0};
	double newinitval[] = {0,0,0};
	
	const char *modeleq[] = {"(c/c(1))^gamma*beta*(alpha*exp(a(1))*k^(alpha-1)+1-delta)=1; a=rho*a(-1)+eps; k+c=exp(a)*k(-1)^alpha+(1-delta)*k(-1);"};

	dynareR(varendo, &numendo, varexo, &numexo, parameters, &numpar, modeleq,
			&ord, "journal", parval, vcov, initval,
			&numsteps, tensorbuffer,
			&num_state, ordering_state, ordering_endo, ordering_exo,
			newinitval,&e);
	printf("error code: %d\n", e);
}
