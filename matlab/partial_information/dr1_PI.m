function [dr,info,M_,options_,oo_] = dr1_PI(dr,task,M_,options_,oo_)
% function [dr,info,M_,options_,oo_] = dr1_PI(dr,task,M_,options_,oo_)
% Computes the reduced form solution of a rational expectation model first
% order
% approximation of the Partial Information stochastic model solver around the deterministic steady state).
% Prepares System as
%        A0*E_t[y(t+1])+A1*y(t)=A2*y(t-1)+c+psi*eps(t)
% with z an exogenous variable process.
% and calls PI_Gensys.m solver 
% based on Pearlman et al 1986 paper and derived from
% C.Sims' gensys linear solver.
% to return solution in format
%       [s(t)' x(t)' E_t x(t+1)']'=G1pi [s(t-1)' x(t-1)' x(t)]'+C+impact*eps(t),
%
% INPUTS
%   dr         [matlab structure] Decision rules for stochastic simulations.
%   task       [integer]          if task = 0 then dr1 computes decision rules.
%                                 if task = 1 then dr1 computes eigenvalues.
%   M_         [matlab structure] Definition of the model.           
%   options_   [matlab structure] Global options.
%   oo_        [matlab structure] Results 
%    
% OUTPUTS
%   dr         [matlab structure] Decision rules for stochastic simulations.
%   info       [integer]          info=1: the model doesn't define current variables uniquely
%                                 info=2: problem in mjdgges.dll info(2) contains error code. 
%                                 info=3: BK order condition not satisfied info(2) contains "distance"
%                                         absence of stable trajectory.
%                                 info=4: BK order condition not satisfied info(2) contains "distance"
%                                         indeterminacy.
%                                 info=5: BK rank condition not satisfied.
%   M_         [matlab structure]            
%   options_   [matlab structure]
%   oo_        [matlab structure]
%  
% ALGORITHM
%   ...
%    
% SPECIAL REQUIREMENTS
%   none.
%  

% Copyright (C) 1996-2010 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

    info = 0;
  
    options_ = set_default_option(options_,'loglinear',0);
    options_ = set_default_option(options_,'noprint',0);
    options_ = set_default_option(options_,'olr',0);
    options_ = set_default_option(options_,'olr_beta',1);
    options_ = set_default_option(options_,'qz_criterium',1.000001);
    
    xlen = M_.maximum_endo_lead + M_.maximum_endo_lag + 1;

    if (options_.aim_solver == 1)
        options_.aim_solver == 0;
       warning('You can not use AIM with Part Info solver. AIM ignored'); 
    end
    if (options_.order > 1) 
        warning('You can not use order higher than 1 with Part Info solver. Order 1 assumed');
        options_.order =1;
    end

    % expanding system for Optimal Linear Regulator
    if options_.ramsey_policy
        if isfield(M_,'orig_model')
            orig_model = M_.orig_model;
            M_.endo_nbr = orig_model.endo_nbr;
            M_.endo_names = orig_model.endo_names;
            M_.lead_lag_incidence = orig_model.lead_lag_incidence;
            M_.maximum_lead = orig_model.maximum_lead;
            M_.maximum_endo_lead = orig_model.maximum_endo_lead;
            M_.maximum_lag = orig_model.maximum_lag;
            M_.maximum_endo_lag = orig_model.maximum_endo_lag;
        end
        old_solve_algo = options_.solve_algo;
        %  options_.solve_algo = 1;
        oo_.steady_state = dynare_solve('ramsey_static',oo_.steady_state,0,M_,options_,oo_,it_);
        options_.solve_algo = old_solve_algo;
        [junk,junk,multbar] = ramsey_static(oo_.steady_state,M_,options_,oo_,it_);
        [jacobia_,M_] = ramsey_dynamic(oo_.steady_state,multbar,M_,options_,oo_,it_);
        klen = M_.maximum_lag + M_.maximum_lead + 1;
        dr.ys = [oo_.steady_state;zeros(M_.exo_nbr,1);multbar];
    else
        klen = M_.maximum_lag + M_.maximum_lead + 1;
        iyv = M_.lead_lag_incidence';
        iyv = iyv(:);
        iyr0 = find(iyv) ;
        it_ = M_.maximum_lag + 1 ;
        
        if M_.exo_nbr == 0
            oo_.exo_steady_state = [] ;
        end
        
        z = repmat(dr.ys,1,klen);
        z = z(iyr0) ;
        [junk,jacobia_] = feval([M_.fname '_dynamic'],z,[oo_.exo_simul ...
              oo_.exo_det_simul], M_.params, it_);
    end
    
    if options_.debug
        save([M_.fname '_debug.mat'],'jacobia_')
    end
    
    dr=set_state_space(dr,M_);
    kstate = dr.kstate;
    kad = dr.kad;
    kae = dr.kae;
    nstatic = dr.nstatic;
    nfwrd = dr.nfwrd;
    npred = dr.npred;
    nboth = dr.nboth;
    order_var = dr.order_var;
    nd = size(kstate,1);
    nz = nnz(M_.lead_lag_incidence);
    
    sdyn = M_.endo_nbr - nstatic;
    
    k0 = M_.lead_lag_incidence(M_.maximum_endo_lag+1,order_var);
    k1 = M_.lead_lag_incidence(find([1:klen] ~= M_.maximum_endo_lag+1),:);
    b = jacobia_(:,k0);

    if (options_.aim_solver == 1) 
      error('Anderson and Moore AIM solver is not compatible with Partial Information models');
    end % end if useAIM and...

    %forward--looking models
    if nstatic > 0
        [Q,R] = qr(b(:,1:nstatic));
        aa = Q'*jacobia_;
    else
        aa = jacobia_;
    end

    % If required, try PCL86 solver, that is, if not the check being
    % performed only and if it is 1st order 
    % create sparse, extended jacobia AA:
    nendo=M_.endo_nbr; % = size(aa,1)
    %%% OLD: aax=zeros(size(aa,1),size(aa,1)*klen);
     % partition jacobian:
    jlen=dr.nspred+dr.nsfwrd+M_.endo_nbr+M_.exo_nbr; % length of jacobian


    % find size xlen of the state vector Y and of A0, A1 and A2 transition matrices:
    % it is the sum the all i variables's lag/lead representations,
    % for each variable i representation being defined as:
    % Max (i_lags-1,0)+ Max (i_leads-1,0)+1
    % so that if variable x appears with 2 lags and 1 lead, and z
    % with 2 lags and 3 leads, the size of the state space is:
    % 1+0+1   +   1+2+1   =6
    % e.g. E_t Y(t+1)=
    %     E_t x(t)
    %     E_t x(t+1)
    %     E_t z(t)
    %     E_t z(t+1)
    %     E_t z(t+2)
    %     E_t z(t+3) 

      
    % first transpose M_.lead_lag_incidence';
    lead_lag=M_.lead_lag_incidence';
    max_lead_lag=zeros(nendo,2); % lead/lag representation in Y for each endogenous variable i
    if ( M_.maximum_lag <= 1) && (M_.maximum_lead <= 1)
      xlen=nendo; %%=0;
      AA0=zeros(xlen,xlen);  % empty A0
      AA2=AA0; % empty A2 and A3
      AA3=AA0;
      AA1=jacobia_(:,npred+1:npred+nendo); 
      fnd = find(lead_lag(:,3));
      AA0(:, fnd)= jacobia_(:,nonzeros(lead_lag(:,3))); %forwd jacobian
      fnd = find(lead_lag(:,1));
      AA2(:, fnd)= jacobia_(:,nonzeros(lead_lag(:,1))); %backward

      if M_.orig_endo_nbr<nendo
        exp_0= strmatch('AUX_EXPECT_LEAD_0_', M_.endo_names);
        num_exp_0=size(exp_0,1);
        if num_exp_0>0
          AA3(:,exp_0)=AA1(:,exp_0);
          XX0=zeros(nendo,num_exp_0);
          AA1(:,exp_0)=XX0(:,[1:num_exp_0])
        end
      end

    else
      xlen=0;
      for i=1:nendo
          llmask=find(lead_lag(i,:)); % mask of leads and lags for var i
          nlag = max((M_.maximum_lag+1-min(llmask)), 0); % reduced no of lags and
          nlead = max((max(llmask)-(M_.maximum_lag+1)), 0); % reduced no of leads for var i
          max_lead_lag(i,:)=[nlag nlead]; % store for future reference
          %xlen=xlen+(nlag+nlead+1); % size as the sum over all the i variables
          xlen=xlen+(max(nlag-1,0)+max(nlead-1,0)+1); % size as the sum over all the i variables
      end
      AA0=zeros(xlen,xlen);  % empty A0
      AA2=AA0; % empty A2 and A3
      AA3=AA0;
    end
    if (xlen>nendo )||( M_.maximum_lag >1) || (M_.maximum_lead >1) % we could not use shortcut above        
        start=xlen -nendo+1;
        offset=0;
        for i=1:nendo
            llmask=find(lead_lag(i,:)); % mask of leads and lags for var i
            nlag=max_lead_lag(i,1); % size for the i'th variable lags
            nlead=max_lead_lag(i,2); % size for the i'th variable lead
            ilen=max(nlag-1,0)+max(nlead-1,0)+1; % size for the i'th variable
            if  lead_lag(i,M_.maximum_lag-nlag+1) && nlag>0 %(j0==1 )&& lead_lag(i,j0) % !=0 %(ilen - iLagLen) %<=max(nlag,1)
                AA2( start:end, offset+1)=aa(:,lead_lag(i,M_.maximum_lag-nlag+1));
            else                        
                AA2( start+i-1, offset+1)=Inf;
            end
            if  lead_lag(i,klen-M_.maximum_lead+nlead) && nlead>0 %(j0==1 )&& lead_lag(i,j0) % !=0 %(ilen - iLagLen) %<=max(nlag,1)
                AA0( start:end, offset+ilen)=aa(:,lead_lag(i,klen-M_.maximum_lead+nlead));
            else                        
                AA0( start+i-1, offset+ilen)=Inf;
            end
                for j0= 1:ilen % M_.maximum_lag
                    if (j0<max(nlag-1,0)+1) && (ilen>nlag)&&  lead_lag(i,j0+1) % !=0 %(ilen - iLagLen) %
                        AA1( start:end, offset+j0)=aa(:,lead_lag(i,(j0+1)));
                    elseif (j0==max(nlag-1,0)+1) && lead_lag(i,(M_.maximum_lag+1)) %&&  (j0<=ilen-max(nlead-1,0) ) ...
                        %AA1( start+i-1, offset+j0)=Inf;  // jacobian at t
                        AA1( start:end, offset+j0)=aa(:,lead_lag(i,(M_.maximum_lag+1)));
                    elseif (j0>max(nlag-1,0)+1)&& (ilen>nlead)&& lead_lag(i,M_.maximum_lag+j0+1)
                        AA1( start:end, offset+j0)=aa(:,lead_lag(i,(M_.maximum_lag+1+j0)));
                    end
                end

            offset=offset+ilen;
            if offset>xlen 
                error(' dr1_PI: offset exceeds max xlen!');
            end
        end
    end
        
        PSI=-[[zeros(xlen-nendo,M_.exo_nbr)];[jacobia_(:, jlen-M_.exo_nbr+1:end)]]; % exog
        cc=0;
        NX=M_.exo_nbr; % no of exogenous varexo shock variables.
        NETA=nfwrd+nboth; % total no of exp. errors  set to no of forward looking equations
        FL_RANK=rank(AA0); % nfwrd+nboth; % min total no of forward looking equations and vars
        
        try

        % call [G1pi,C,impact,nmat,TT1,TT2,gev,eu]=PI_gensys(a0,a1,a2,c,PSI,NX,NETA,NO_FL_EQS)
        % System given as
        %        a0*E_t[y(t+1])+a1*y(t)=a2*y(t-1)+c+psi*eps(t)
        % with eps an exogenous variable process.
        % Returned system is
        %       [s(t)' x(t)' E_t x(t+1)']'=G1pi [s(t-1)' x(t-1)' x(t)]'+C+impact*eps(t),
        %  and (a) the matrix nmat satisfying   nmat*E_t z(t)+ E_t x(t+1)=0
        %      (b) matrices TT1, TT2  that relate y(t) to these states:
        %      y(t)=[TT1 TT2][s(t)' x(t)']'.
            if(options_.ACES_solver==1)
                SSbar= diag(dr.ys);%(oo_.steady_state);
                AA0=AA0*SSbar;
                AA1=AA1*SSbar;
                AA2=AA2*SSbar;
                AA3=AA3*SSbar;
            end
            %%if any(AA3)  % for expectational models when complete
            %%  [G1pi,CC,impact,nmat,TT1,TT2,gev,eu, DD, E3,E5]=PI_gensysEXP(AA0,AA1,-AA2,AA3,cc,PSI,NX,NETA,FL_RANK, M_, options_);
            %%else
              [G1pi,CC,impact,nmat,TT1,TT2,gev,eu, DD, E3,E5]=PI_gensys(AA0,AA1,-AA2,AA3,cc,PSI,NX,NETA,FL_RANK, M_, options_);
            %%end

            % reuse some of the bypassed code and tests that may be needed 
            if eu ~=[1; 1]
                info(1) = abs(eu(1)+eu(2));
                info(2) = 1.0e+8;
           %     return
            end
            
            dr.PI_ghx=G1pi;
            dr.PI_ghu=impact;
            dr.PI_TT1=TT1;
            dr.PI_TT2=TT2;
            dr.PI_nmat=nmat;
            dr.PI_CC=CC;
            dr.PI_gev=gev;
            dr.PI_eu=eu;
            dr.PI_FL_RANK=FL_RANK;
            %dr.ys=zeros(nendo); % zero steady state
            dr.ghx=G1pi;
            dr.ghu=impact;
            dr.eigval = eig(G1pi);
            dr.rank=FL_RANK;
            
            if options_.ACES_solver==1
                ACES.A=G1pi;
                ACES.C=impact; % (:,1);
                ACES.D=DD; %=impact (:,20);
                ACES.E3=E3;
                ACES.E5=E5;
                save ACES ;%ACES_A ACES_C ACES_D ACES_E2 ACES_E5
                %save([M_.fname '_ACES_IN'], 'ACES')
                save ([M_.fname '_ACES_A.txt'], 'G1pi', '-ascii', '-double', '-tabs');
                save ([M_.fname '_ACES_C.txt'], 'impact','-ascii', '-double', '-tabs');
                save ([M_.fname '_ACES_D.txt'], 'DD', '-ascii', '-double', '-tabs');
                save ([M_.fname '_ACES_E3.txt'], 'E3', '-ascii', '-double', '-tabs');
                save ([M_.fname '_ACES_E5.txt'], 'E5', '-ascii', '-double', '-tabs');
            end

        catch
            disp('Problem with using Part Info solver - Using Dynare solver instead');
            lerror=lasterror;
            disp (lerror);
            options_.PartInfo = 0; % and then try mjdgges instead
        end

    % TODO: 
    % if options_.loglinear == 1
    % if exogenous deterministic variables
        
return;
