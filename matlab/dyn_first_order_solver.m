function [dr,info] = dyn_first_order_solver(jacobia,DynareModel,dr,DynareOptions,task)

%@info:
%! @deftypefn {Function File} {[@var{dr},@var{info}] =} dyn_first_order_solver (@var{jacobia},@var{DynareModel},@var{dr},@var{DynareOptions},@var{task})
%! @anchor{dyn_first_order_solver}
%! @sp 1
%! Computes the first order reduced form of the DSGE model
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item jacobia
%! Matrix containing the Jacobian of the model
%! @item DynareModel
%! Matlab's structure describing the model (initialized by @code{dynare}).
%! @item dr
%! Matlab's structure describing the reduced form solution of the model.
%! @item qz_criterium
%! Double containing the criterium to separate explosive from stable eigenvalues
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item dr
%! Matlab's structure describing the reduced form solution of the model.
%! @item info
%! Integer scalar, error code.
%! @sp 1
%! @table @ @code
%! @item info==0
%! No error.
%! @item info==1
%! The model doesn't determine the current variables uniquely.
%! @item info==2
%! MJDGGES returned an error code.
%! @item info==3
%! Blanchard & Kahn conditions are not satisfied: no stable equilibrium.
%! @item info==4
%! Blanchard & Kahn conditions are not satisfied: indeterminacy.
%! @item info==5
%! Blanchard & Kahn conditions are not satisfied: indeterminacy due to rank failure.
%! @item info==7
%! One of the generalized eigenvalues is close to 0/0
%! @end table
%! @end table
%! @end deftypefn
%@eod:

% Copyright (C) 2001-2012 Dynare Team
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

persistent reorder_jacobian_columns innovations_idx index_s index_m index_c index_p row_indx index_0m index_0p k1 k2 j3 j4
persistent ndynamic nstatic nfwrd npred nboth nd nyf n


if ~nargin
    reorder_jacobian_columns = [];
    dr = [];
    info = [];
    return
end

if isempty(reorder_jacobian_columns)

    kstate   = dr.kstate;
    nfwrd    = dr.nfwrd;
    nboth    = dr.nboth;
    npred    = dr.npred-nboth;
    nstatic  = dr.nstatic;
    ndynamic = npred+nfwrd+nboth;
    nyf      = nfwrd + nboth;
    n        = ndynamic+nstatic;

    k1 = 1:(npred+nboth);
    k2 = 1:(nfwrd+nboth);

    order_var = dr.order_var;
    nd = size(kstate,1);
    lead_lag_incidence = DynareModel.lead_lag_incidence;
    nz = nnz(lead_lag_incidence);

    %lead variables actually present in the model
    j4 = nstatic+npred+1:nstatic+npred+nboth+nfwrd;   % Index on the forward and both variables
    j3 = nonzeros(lead_lag_incidence(2,j4)) - nstatic - 2 * npred - nboth;  % Index on the non-zeros forward and both variables
    j4 = find(lead_lag_incidence(2,j4));

    no_lead_id = find(lead_lag_incidence(3,:)==0);
    no_lag_id = find(lead_lag_incidence(1,:)==0);

    static_id = intersect(no_lead_id,no_lag_id);
    lag_id = setdiff(no_lead_id,static_id);
    lead_id = setdiff(no_lag_id,static_id);
    both_id = intersect(setdiff(1:n,no_lead_id),setdiff(1:n,no_lag_id));

    lead_idx = lead_lag_incidence(3,lead_id);
    lag_idx = lead_lag_incidence(1,lag_id);
    both_lagged_idx = lead_lag_incidence(1,both_id);
    both_leaded_idx = lead_lag_incidence(3,both_id);
    innovations_idx = (size(jacobia,2)-DynareModel.exo_nbr+1):size(jacobia,2);
    dr.state_var  = [lag_idx, both_lagged_idx];

    indexi_0 = 0;
    if DynareModel.maximum_endo_lag > 0 && (npred > 0  || nboth > 0)
        indexi_0 = min(lead_lag_incidence(2,:));
    end

    index_c  = lead_lag_incidence(2,:);             % Index of all endogenous variables present at time=t
    index_s  = lead_lag_incidence(2,1:nstatic);     % Index of all static endogenous variables present at time=t
    index_0m = (nstatic+1:nstatic+npred)+indexi_0-1;
    index_0p = (nstatic+npred+1:n)+indexi_0-1;
    index_m  = 1:(npred+nboth);
    index_p  = lead_lag_incidence(3,find(lead_lag_incidence(3,:)));
    row_indx = nstatic+1:n;

    reorder_jacobian_columns = [lag_idx, both_lagged_idx, npred+nboth+[static_id lag_id both_id lead_id], both_leaded_idx, lead_idx, innovations_idx ];

end

info = 0;

dr.ghx = [];
dr.ghu = [];

jacobia = jacobia(:,reorder_jacobian_columns);

if nstatic > 0
    [Q, junk] = qr(jacobia(:,index_s));
    aa = Q'*jacobia;
else
    aa = jacobia;
end

A = aa(:,index_m);  % Jacobain matrix for lagged endogeneous variables
B = aa(:,index_c);  % Jacobian matrix for contemporaneous endogeneous variables
C = aa(:,index_p);  % Jacobain matrix for led endogeneous variables

if task ~= 1 && DynareOptions.dr_cycle_reduction == 1
    A1 = [aa(row_indx,index_m ) zeros(ndynamic,nfwrd)];
    B1 = [aa(row_indx,index_0m) aa(row_indx,index_0p) ];
    C1 = [zeros(ndynamic,npred) aa(row_indx,index_p)];
    [ghx, info] = cycle_reduction(A1, B1, C1, DynareOptions.dr_cycle_reduction_tol);
    ghx = ghx(:,index_m);
    hx = ghx(1:npred+nboth,:);
    gx = ghx(1+npred:end,:);
end

if (task ~= 1 && ((DynareOptions.dr_cycle_reduction == 1 && info ==1) || DynareOptions.dr_cycle_reduction == 0)) || task == 1
    D = [[aa(row_indx,index_0m) zeros(ndynamic,nboth) aa(row_indx,index_p)] ; [zeros(nboth, npred) eye(nboth) zeros(nboth, nboth + nfwrd)]];
    E = [-aa(row_indx,[index_m index_0p])  ; [zeros(nboth,nboth+npred) eye(nboth,nboth+nfwrd) ] ];

    [err, ss, tt, w, sdim, dr.eigval, info1] = mjdgges(E,D,DynareOptions.qz_criterium);
    mexErrCheck('mjdgges', err);

    if info1
        if info1 == -30
            % one eigenvalue is close to 0/0
            info(1) = 7;
        else
            info(1) = 2;
            info(2) = info1;
            info(3) = size(E,2);
        end
        return
    end

    nba = nd-sdim;

    if task == 1
        dr.rank = rank(w(1:nyf,nd-nyf+1:end));
        % Under Octave, eig(A,B) doesn't exist, and
        % lambda = qz(A,B) won't return infinite eigenvalues
        if ~exist('OCTAVE_VERSION')
            dr.eigval = eig(E,D);
        end
        return
    end

    if nba ~= nyf
        temp = sort(abs(dr.eigval));
        if nba > nyf
            temp = temp(nd-nba+1:nd-nyf)-1-DynareOptions.qz_criterium;
            info(1) = 3;
        elseif nba < nyf;
            temp = temp(nd-nyf+1:nd-nba)-1-DynareOptions.qz_criterium;
            info(1) = 4;
        end
        info(2) = temp'*temp;
        return
    end

    %First order approximation
    if task ~= 1
        indx_stable_root = 1: (nd - nyf);     %=> index of stable roots
        indx_explosive_root = npred + nboth + 1:nd;  %=> index of explosive roots
                                                     % derivatives with respect to dynamic state variables
                                                     % forward variables
        Z = w';
        Z11t = Z(indx_stable_root,    indx_stable_root)';
        Z21  = Z(indx_explosive_root, indx_stable_root);
        Z22  = Z(indx_explosive_root, indx_explosive_root);
        if ~isfloat(Z21) && (condest(Z21) > 1e9)
            info(1) = 5;
            info(2) = condest(Z21);
            return;
        else
            gx = - Z22 \ Z21;
        end
        % predetermined variables
        hx =  Z11t * inv(tt(indx_stable_root, indx_stable_root)) * ss(indx_stable_root, indx_stable_root) * inv(Z11t);
        ghx = [hx(k1,:); gx(k2(nboth+1:end),:)];
    end;
end

if  task~= 1

    if nstatic > 0
        B_static = B(:,1:nstatic);  % submatrix containing the derivatives w.r. to static variables
    else
        B_static = [];
    end;
    %static variables, backward variable, mixed variables and forward variables
    B_pred = B(:,nstatic+1:nstatic+npred+nboth);
    B_fyd = B(:,nstatic+npred+nboth+1:end);

    % static variables
    if nstatic > 0
        temp = - C(1:nstatic,j3)*gx(j4,:)*hx;
        b = aa(:,index_c);
        b10 = b(1:nstatic, 1:nstatic);
        b11 = b(1:nstatic, nstatic+1:n);
        temp(:,index_m) = temp(:,index_m)-A(1:nstatic,:);
        temp = b10\(temp-b11*ghx);
        ghx = [temp; ghx];
        temp = [];
    end

    A_ = real([B_static C(:,j3)*gx+B_pred B_fyd]); % The state_variable of the block are located at [B_pred B_both]

    if DynareModel.exo_nbr
        if nstatic > 0
            fu = Q' * jacobia(:,innovations_idx);
        else
            fu = jacobia(:,innovations_idx);
        end;

        ghu = - A_ \ fu;
    else
        ghu = [];
    end;
end


dr.ghx = ghx;
dr.ghu = ghu;

if DynareOptions.aim_solver ~= 1 && DynareOptions.use_qzdiv
    % Necessary when using Sims' routines for QZ
    dr.ghx = real(ghx);
    dr.ghu = real(ghu);
    hx = real(hx);
end

dr.Gy = hx;