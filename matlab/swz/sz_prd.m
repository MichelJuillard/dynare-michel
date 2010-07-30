function sz_prd(M,options_)
%==========================================================================
%== Directory structure
%==========================================================================

%generation of mhm file
%generateMHM_input(options_);

swz_root = strrep(which('swz_sbvar'),'/swz_sbvar.m','');

% path for C executables
%c_path='./c-executables';
c_path = [swz_root '/bin'];

% path for Markov specification
m_spec_path=[swz_root '/switching_specification'];

% path for MHM specification
mhm_spec_path=[swz_root '/mhm_specification'];

%==========================================================================
%== Processing control
%==========================================================================
% 1 to use standalone, 0 to use mex
options_.ms.standalone = 0;

%==========================================================================
%== Output control
%==========================================================================
% tag for output files %create an option
%options_.ms.output_file_tag='test_2v';

% 1 to create init_<options_.ms.output_file_tag>.dat, 0 otherwise
%options_.ms.create_initialization_file = 1; %0 originally

% 1 to perform estimation, 0 otherwise
%options_.ms.estimate_msmodel = 1;

% 1 to perform estimation, 0 otherwise
%options_.ms.compute_mdd = 1;

% 1 to compute probabilites, 0 otherwise
%options_.ms.compute_probabilities = 1;%1 in the original

% 1 to Prints draws of the posterior
%options_.ms.print_draws = 1;
%options_.ms.n_draws=1000;
%options_.ms.thinning_factor=1;

%==========================================================================
%== Markov Process Specification File
%==========================================================================
%options_.ms.markov_file = 'specification_2v2c.dat'; %create an option
markov_file = [options_.ms.markov_file '.dat'];

%==========================================================================
%== Markov Process Specification File
%==========================================================================
%options_.ms.mhm_file = 'MHM_input.dat';


mhm_file = '/MHM_input.dat';
%options_.ms.proposal_draws = 100000;

%==========================================================================
%== Var Specification
%==========================================================================
% Number of options_.ms.nlags
%options_.ms.nlags = 4;

% Var restriction function
%options_.ms.restriction_fname = 'ftd_upperchol3v'; %create an option


%==========================================================================
%== BVAR prior
%==========================================================================

%=== The following mu is effective only if indxPrior==1.
mu = zeros(6,1);   % hyperparameters

% mu(1): overall tightness for A0 and Aplus
mu(1) = 1.0;

% mu(2): relative tightness for Aplus
mu(2) = 1.0;

% mu(3): relative tightness for the constant term
mu(3) = 0.1;

% mu(4): tightness on lag decay.  (1.2 - 1.5 faster decay produces better
% inflation forecasts
mu(4) = 1.2;  

% mu(5): weight on nvar sums of coeffs dummy observations (unit roots).
mu(5) = 1; 

% mu(6): weight on single dummy initial observation including constant
% (cointegration, unit roots, and stationarity).
mu(6) = 1; 

% Alpha on p. 66 for squared time-varying structural shock lambda.
galp = 1.0; 

% Beta on p. 66 for squared time-varying structural shock lambda.
gbeta = 1.0;   

% Case 3 (no state change across options_.ms.nlags (l) but allows all variables for a
% given lag to switch states). Normal prior variance for glamda
% (nvar-by-nvar for each state) for different variables in lagged D+.  See
% p.71v.  
gsig2_lmdm = 50^2;  


%==========================================================================
%== Data
%==========================================================================
% Read in data to produce rectangular array named xdd.  Each column is one
% data series.
%load ./data/artificial_data
xdd=options_.data;

% Information about timing of the data for consistancy checks
% quarters (4) or months (12)
%q_m = 4;  
%options_.ms.freq = 4;
q_m = options_.ms.freq;
% beginning year in data set
%yrBin=1978;   
%options_.ms.initial_year = 1959;
yrBin=options_.ms.initial_year;
% beginning quarter or month in data set
%qmBin=1;  
%options_.ms.initial_subperiod = 1;
qmBin=options_.ms.initial_subperiod;
% final year in data set
%yrFin=2007;  
%options_.ms.final_year = 2005;
yrFin=options_.ms.final_year;
% final month or quarter in data set
%qmFin=4;  
%options_.ms.final_subperiod = 4;
qmFin=options_.ms.final_subperiod;
% first year to use in estimation
%yrStart=yrBin;
yrStart=options_.ms.initial_year;
% first quarter or month to use in estimation
%qmStart=qmBin;
qmStart=options_.ms.initial_subperiod;
% last year to use in estimation
%yrEnd=yrFin;
yrEnd=options_.ms.final_year;
% last quater or month to use in estimation
%qmEnd=qmFin;
qmEnd=options_.ms.final_subperiod;
% Log variables in xdd
logindx = []; 

% Convert percent to decimal in xdd
pctindx = [];  

% Select the variable to use and rearrange columns if desired
%vlist = [3 1 2];
%options_.ms.vlist = [1 2 3];
options_.ms.vlist = [1:size(options_.varobs,1)];
vlist1=options_.ms.vlist;

%==========================================================================
%== Linux or Windows system - differences in some naming conventions
%==========================================================================
use_linux = 1;


%==========================================================================
%== Random number seed for selecting starting point in constant parameter 
%== optimization.
%==========================================================================
% Set to zero to set from clock
seednumber = 7910;
if seednumber
    randn('state',seednumber);
    rand('state',seednumber);
else
    rand('state',fix(100*sum(clock)));
    randn('state',1000000000*rand);
    rand('state',1000000000*rand);
end

%==========================================================================
%==========================================================================
%==========================================================================
%== Beginning of code.  Modify below at own risk.
%==========================================================================

% options that may at some point become user specified
%indxC0Pres = 0;    % 1: cross-A0-and-A+ restrictions; 0: idfile_const is all we have
indxC0Pres =options_.ms.cross_restrictions;
                   % Example for indxOres==1: restrictions of the form P(t) = P(t-1).
%Rform = 0;         % 1: contemporaneous recursive reduced form; 0: restricted (non-recursive) form
Rform =options_.ms.contemp_reduced_form;
% % % Pseudo = 0;  % 1: Pseudo forecasts; 0: real time forecasts
%indxPrior = 1;     % 1: Bayesian prior; 0: no prior
indxPrior =options_.ms.bayesian_prior;
%indxDummy = indxPrior;  % 1: add dummy observations to the data; 0: no dummy added.
indxDummy = options_.ms.bayesian_prior;
%ndobs = 0;         % No dummy observations for xtx, phi, fss, xdatae, etc.  Dummy observations are used as an explicit prior in fn_rnrprior_covres_dobs.m.
ndobs =options_.ms.dummy_obs;
%if indxDummy
%   ndobs=nvar+1;         % number of dummy observations
%else
%   ndobs=0;    % no dummy observations
%end

% 
hpmsmd = [0.0; 0.0];
indxmsmdeqn = [0; 0; 0; 0];  %This option disenable using this in fn_rnrprior_covres_dobs.m

nStates = -1;


%==========================================================================
%== Create initialization file
%==========================================================================
if options_.ms.create_initialization_file == 1
    %======================================================================
    %== Check and setup data
    %======================================================================
    % log data
    xdd(:,logindx) = log(xdd(:,logindx));

    % convert percentage to decimal
    xdd(:,pctindx)=.01*xdd(:,pctindx);

    if (q_m ~= 12) && (q_m ~= 4)
        disp('Warning: data must be monthly or quarterly!')
        return
    end

    % number of data points
    nData=(yrFin-yrBin)*q_m + (qmFin-qmBin+1);
    % number of data points in estimation sample
    nSample=(yrEnd-yrStart)*q_m + (qmEnd-qmEnd+1);
    % number of periods not used at beginning of data (non-negative number)
    nStart=(yrStart-yrBin)*q_m + (qmStart-qmBin);
    % number of periods not used at end of data (non-positive number)
    nEnd=(yrEnd-yrFin)*q_m + (qmEnd-qmFin);

    if (nEnd > 0) || (nStart < 0)
        disp('Warning: desired estimation period not in data set!')
        return
    end
    if (nSample <= 0)
        disp('Warning: no data points in estimation period!')
        return
    end

    % reorder variables and create estimation data set
    xdgel=xdd(nStart+1:nData+nEnd,vlist1);

    % bad data points
    baddata = find(isnan(xdgel));
    if ~isempty(baddata)
        disp('Warning: some data for estimation period are unavailable.')
        return
    end

    % set nvar and nexo
    nvar=size(xdgel,2);
    nexo=1;

    % Arranged data information, WITHOUT dummy obs when 0 after mu is used.
    % See fn_rnrprior_covres_dobs.m for using the dummy observations as part of
    % an explicit prior.
    [xtx,xty,yty,fss,phi,y,ncoef,xr,Bh] = fn_dataxy(nvar,options_.ms.nlags,xdgel,mu,0,nexo);


    %======================================================================
    %== Linear Restrictions
    %======================================================================
    if Rform
        Ui=cell(nvar,1); Vi=cell(ncoef,1);
        for kj=1:nvar
            Ui{kj} = eye(nvar);  Vi{kj} = eye(ncoef);
        end
    else
        [Ui,Vi,n0,np,ixmC0Pres] = feval(options_.ms.restriction_fname,nvar,nexo,options_.ms);
        if min(n0)==0
            disp('A0: restrictions give no free parameters in one of equations')
            return
        elseif min(np)==0
            disp('Aplus: Restrictions in give no free parameters in one of equations')
            return
        end
    end


    %======================================================================
    %== Estimation
    %======================================================================
    if indxPrior
        %*** Obtains asymmetric prior (with no linear restrictions) with dummy observations as part of an explicit prior (i.e,
        %      reflected in Hpmulti and Hpinvmulti).  See Forecast II, pp.69a-69b for details.
        if 1  % Liquidity effect prior on both MS and MD equations.
            [Pi,H0multi,Hpmulti,H0invmulti,Hpinvmulti] = fn_rnrprior_covres_dobs(nvar,q_m,options_.ms.nlags,xdgel,mu,indxDummy,hpmsmd,indxmsmdeqn);
        else
            [Pi,H0multi,Hpmulti,H0invmulti,Hpinvmulti] = fn_rnrprior(nvar,q_m,options_.ms.nlags,xdgel,mu);
        end

        %*** Combines asymmetric prior with linear restrictions
        [Ptld,H0invtld,Hpinvtld] = fn_rlrprior(Ui,Vi,Pi,H0multi,Hpmulti,nvar);

        %*** Obtains the posterior matrices for estimation and inference
        [Pmat,H0inv,Hpinv] = fn_rlrpostr(xtx,xty,yty,Ptld,H0invtld,Hpinvtld,Ui,Vi);
    else
        %*** Obtain the posterior matrices for estimation and inference
        [Pmat,H0inv,Hpinv] = fn_dlrpostr(xtx,xty,yty,Ui,Vi);
    end

    if Rform
        %*** Obtain the ML estimate
        A0hatinv = chol(H0inv{1}/fss);   % upper triangular but lower triangular choleski
        A0hat=inv(A0hatinv);

        Aphat = Pmat{1}*A0hat;
    else
        %*** Obtain the ML estimate
        %   load idenml
        x = 10*rand(sum(n0),1);
        H0 = eye(sum(n0));
        crit = 1.0e-9;
        nit = 10000;
        %
        tic
        [fhat,xhat,grad,Hhat,itct,fcount,retcodehat] = csminwel('fn_a0freefun',x,H0,'fn_a0freegrad',crit,nit,Ui,nvar,n0,fss,H0inv);
        endtime = toc;

        A0hat = fn_tran_b2a(xhat,Ui,nvar,n0);

        xhat = fn_tran_a2b(A0hat,Ui,nvar,n0);
        [Aphat,ghat] = fn_gfmean(xhat,Pmat,Vi,nvar,ncoef,n0,np);
        if indxC0Pres
            Fhatur0P = Fhat;  % ur: unrestriced across A0 and A+
            for ki = 1:size(ixmC0Pres,1)   % loop through the number of equations in which
                % cross-A0-A+ restrictions occur. See St. Louis Note p.5.
                ixeq = ixmC0Pres{ki}(1,1);   % index for the jth equation in consideration.
                Lit = Vi{ixeq}(ixmC0Pres{ki}(:,2),:);  % transposed restriction matrix Li
                % V_j(i,:) in f_j(i) = V_j(i,:)*g_j
                ci = ixmC0Pres{ki}(:,4) .* A0hat(ixmC0Pres{ki}(:,3),ixeq);
                % s * a_j(h) in the restriction f_j(i) = s * a_j(h).
                LtH = Lit/Hpinv{ixeq};
                HLV = LtH'/(LtH*Lit');
                gihat = Vi{ixeq}'*Fhatur0P(:,ixeq);
                Aphat(:,ixeq) = Vi{ixeq}*(gihat + HLV*(ci-Lit*gihat));
            end
        end
    end


    %======================================================================
    %== Create matlab initialization file
    %======================================================================
    matlab_filename = ['matlab_',options_.ms.output_file_tag,'.prn'];
    fidForC = fopen(matlab_filename,'w');

    fprintf(fidForC,'\n%s\n','//== gxia: alpha parameter for gamma prior of xi ==//');
    fprintf(fidForC,' %20.15f ', galp);
    fprintf(fidForC, '\n\n');

    fprintf(fidForC,'\n%s\n','//== gxib: beta parameter for gamma prior of xi ==//');
    fprintf(fidForC,' %20.15f ', gbeta);
    fprintf(fidForC, '\n\n');

    fprintf(fidForC,'\n%s\n','//== glamdasig: sigma parameter for normal prior of lamda ==//');
    fprintf(fidForC,' %20.15f ', sqrt(gsig2_lmdm));
    fprintf(fidForC, '\n\n');

    %=== lags, nvar, nStates, sample size (excluding options_.ms.nlags where, with dummyies, fss=nSample-options_.ms.nlags+ndobs).
    fprintf(fidForC,'\n%s\n','//== lags, nvar, nStates, T ==//');
    fprintf(fidForC,' %d  %d  %d  %d\n\n\n',options_.ms.nlags, nvar, nStates, fss);

    %=== A0hat nvar-by-nvar from the constant VAR.
    fprintf(fidForC,'\n%s\n','//== A0hat: nvar-by-nvar ==//');
    indxFloat = 1;
    xM = A0hat;
    nrows = nvar;
    ncols = nvar;
    fn_fprintmatrix(fidForC, xM, nrows, ncols, indxFloat)

    %=== Aphat ncoef-by-nvar from the constant VAR.
    %=== Each column of Aphat is in the order of [nvar variables for 1st lag, ..., nvar variables for last lag, constant term].
    fprintf(fidForC,'\n%s\n','//== Aphat: ncoef(lags*nvar+1)-by-nvar ==//');
    indxFloat = 1;
    xM = Aphat;
    nrows = ncoef;
    ncols = nvar;
    fn_fprintmatrix(fidForC, xM, nrows, ncols, indxFloat)

    %=== n0const: nvar-by-1, whose ith element represents the number of free A0 parameters in ith equation for the case of constant parameters.
    fprintf(fidForC,'\n%s\n','//== n0const: nvar-by-1 ==//');
    indxFloat = 0;
    xM = n0;
    nrows = 1;
    ncols = nvar;
    fn_fprintmatrix(fidForC, xM', nrows, ncols, indxFloat)

    %=== npconst: nvar-by-1, whose ith element represents the number of free A+ parameters in ith equation for the case of constant parameters.
    fprintf(fidForC,'\n%s\n','//== npconst: nvar-by-1 ==//');
    indxFloat = 0;
    xM = np;
    nrows = 1;
    ncols = nvar;
    fn_fprintmatrix(fidForC, xM', nrows, ncols, indxFloat)

    %=== Uiconst: nvar-by-1 cell.  In each cell, nvar-by-qi orthonormal basis for the null of the ith
    %           equation contemporaneous restriction matrix where qi is the number of free parameters.
    %           With this transformation, we have ai = Ui*bi or Ui'*ai = bi where ai is a vector
    %           of total original parameters and bi is a vector of free parameters. When no
    %           restrictions are imposed, we have Ui = I.  There must be at least one free
    %           parameter left for the ith equation.
    fprintf(fidForC,'\n%s\n','//== Uiconst: cell(nvar,1) and nvar-by-n0const(i) for the ith cell (equation) ==//');
    for i_=1:nvar
        fn_fprintmatrix(fidForC, Ui{i_}, nvar, n0(i_), 1);
    end

    %=== Viconst: nvar-by-1 cell.  In each cell, k-by-ri orthonormal basis for the null of the ith
    %           equation lagged restriction matrix where k is a total of exogenous variables and
    %           ri is the number of free parameters. With this transformation, we have fi = Vi*gi
    %           or Vi'*fi = gi where fi is a vector of total original parameters and gi is a
    %           vector of free parameters. There must be at least one free parameter left for
    %           the ith equation.
    fprintf(fidForC,'\n%s\n','//== Viconst: cell(nvar,1) and ncoef-by-n0const(i) for the ith cell (equation) ==//');
    for i_=1:nvar
        fn_fprintmatrix(fidForC, Vi{i_}, ncoef, np(i_), 1);
    end

    %=== H0barconstcell: cell(nvar,1) (equations) and n-by-n for each cell (equaiton).
    %=== H0barconst:  prior covariance matrix for each column of A0 under asymmetric prior (including SZ dummy obs.) with NO linear restrictions imposed yet.
    fprintf(fidForC,'\n%s\n','//== H0barconstcell: cell(nvar,1) and n-by-n for the ith cell (equation) ==//');
    for i_=1:nvar
        fn_fprintmatrix(fidForC, H0multi(:,:,i_), nvar, nvar, 1);
    end

    %=== Hpbarconstcell: cell(nvar,1) (equations) and ncoef-by-ncoef for each cell (equaiton).
    %=== Hpbarconst:  prior covariance matrix for each column of A+ under asymmetric prior (including SZ dummy obs.) with NO linear restrictions imposed yet.
    fprintf(fidForC,'\n%s\n','//== Hpbarconstcell: cell(nvar,1) and ncoef-by-ncoef for the ith cell (equation) ==//');
    for i_=1:nvar
        fn_fprintmatrix(fidForC, Hpmulti(:,:,i_), ncoef, ncoef, 1);
    end

    %=== phi:  X; T-by-k; column: [nvar for 1st lag, ..., nvar for last lag, other exogenous terms, const term]
    fprintf(fidForC,'\n%s\n','//== Xright -- X: T-by-ncoef ==//');
    xM = phi;
    nrows = fss;
    ncols = ncoef;
    for ki=1:nrows
        for kj=1:ncols
            fprintf(fidForC,' %20.15f ',xM((kj-1)*nrows+ki));
            if (kj==ncols)
                fprintf(fidForC,'\n');
            end
        end
        if (ki==nrows)
            fprintf(fidForC,'\n\n');
        end
    end

    %=== y:    Y: T-by-nvar where T=fss
    fprintf(fidForC,'\n%s\n','//== Yleft -- Y: T-by-nvar ==//');
    xM = y;
    nrows = fss;
    ncols = nvar;
    for ki=1:nrows
        for kj=1:ncols
            fprintf(fidForC,' %20.15f ',xM((kj-1)*nrows+ki));
            if (kj==ncols)
                fprintf(fidForC,'\n');
            end
        end
        if (ki==nrows)
            fprintf(fidForC,'\n\n');
        end
    end

    fclose(fidForC);

    %======================================================================
    %== Create C initialization filename
    %======================================================================
    swz_write_markov_file(markov_file,M,options_)
    if options_.ms.standalone == 1
        if use_linux == 1
            create_init_file=[c_path,'/sbvar_init_file ',matlab_filename,' ',markov_file,' ',options_.ms.output_file_tag];
            system(create_init_file); %Run operating system command and return result
        else
            create_init_file=[c_path,'\sbvar_init.exe ',matlab_filename,' ',markov_file,' ',options_.ms.output_file_tag];
            dos(create_init_file)            
        end
    else           
        create_init_file=[matlab_filename,' ',markov_file,' ',options_.ms.output_file_tag];
        mex_sbvar_init_file(create_init_file);
    end
end

%==========================================================================
%== Perform estimation
%==========================================================================
if options_.ms.estimate_msmodel == 1
    if options_.ms.standalone == 1
        if use_linux == 1
            perform_estimation=[c_path,'/sbvar_estimation -cseed 5 -ft ',options_.ms.output_file_tag];
            system(perform_estimation);
        else
            perform_estimation=[c_path,'\sbvar_estimation.exe -cseed 5 -ft ',options_.ms.output_file_tag];
            dos(perform_estimation)            
        end
    else
        perform_estimation=['-cseed 5 -ft ',options_.ms.output_file_tag];
        mex_sbvar_estimation(perform_estimation);
    end
end


%==========================================================================
%== Compute marginal data density
%==========================================================================
if options_.ms.compute_mdd == 1
    mhm_file = ['mhm_input_' M.fname '.dat'];
    swz_write_mhm_input(mhm_file,options_.ms);
    if options_.ms.standalone == 1
        if use_linux == 1
            compute_mdd1=[c_path,'/sbvar_mhm_1 -cseed 5 -ft ',options_.ms.output_file_tag,' -fi ',mhm_file];
            system(compute_mdd1);
            compute_mdd2=[c_path,'/sbvar_mhm_2 -cseed 5 -ft ',options_.ms.output_file_tag,' -d ',int2str(options_.ms.proposal_draws),' -t 3'];
            system(compute_mdd2);
        else
            compute_mdd1=[c_path,'\sbvar_mhm_1.exe -cseed 5 -ft ',options_.ms.output_file_tag,' -fi ',mhm_file];
            system(compute_mdd1);
            compute_mdd2=[c_path,'\sbvar_mhm_2.exe -cseed 5 -ft ',options_.ms.output_file_tag,' -d ',int2str(options_.ms.proposal_draws),' -t 3'];
            system(compute_mdd2);            
        end
    else
        compute_mdd1=['-cseed 5 -ft ',options_.ms.output_file_tag,' -fi ',mhm_file];
        mex_sbvar_mhm_1(compute_mdd1);
        compute_mdd2=['-cseed 5 -ft ',options_.ms.output_file_tag,' -d ',int2str(options_.ms.proposal_draws),' -t 3'];
        mex_sbvar_mhm_2(compute_mdd2);
    end
end


%==========================================================================
%== Compute posterior mode regime probabilities
%==========================================================================
if options_.ms.compute_probabilities == 1 %error registers here
    if options_.ms.standalone == 1
        if use_linux == 1
            compute_prob=[c_path,'/sbvar_probabilities -ft ',options_.ms.output_file_tag];
            system(compute_prob);
        else
            compute_prob=[c_path,'\sbvar_probabilities -ft ',options_.ms.output_file_tag];
            system(compute_prob);            
        end
    else
        compute_prob=['-ft ',options_.ms.output_file_tag];
        mex_sbvar_probabilities(compute_prob);
    end
end

%==========================================================================
%== Print Draws
%==========================================================================
if options_.ms.print_draws == 1 %error here as well
    if options_.ms.standalone == 1
        if use_linux == 1
            print_draws=[c_path,'/sbvar_draws -cseed 5 -ft ',options_.ms.output_file_tag,' -i ',int2str(options_.ms.n_draws),' -t ',int2str(options_.ms.thinning_factor)];
            system(print_draws);
        else
            print_draws=[c_path,'\sbvar_draws -cseed 5 -ft ',options_.ms.output_file_tag,' -i ',int2str(options_.ms.n_draws),' -t ',int2str(options_.ms.thinning_factor)];
            system(print_draws);            
        end
    else
        print_draws=['-cseed 5 -ft ',options_.ms.output_file_tag,' -i ',int2str(options_.ms.n_draws),' -t ',int2str(options_.ms.thinning_factor)];
        mex_sbvar_draws(print_draws);
    end
end