    @#define N = 4
    var A B;
    varexo eB;
    parameters Bss rho;

    Bss=1;
    rho=0.9;

    model;
    A = B(@{N});
    B = (1-rho)*Bss +rho*B(-1) +eB;
    end;

    steady_state_model;
    A = Bss;
    B = Bss;
    end;

    resid;
    steady;
    check;

    shocks;
    var eB ; stderr 1 ;
    end;

    stoch_simul(order=1,irf=0) ;
