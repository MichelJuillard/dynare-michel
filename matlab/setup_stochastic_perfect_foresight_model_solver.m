function pfm = setup_stochastic_perfect_foresight_model_solver(DynareModel,DynareOptions,DynareOutput,IntegrationMethod)

% Copyright (C) 2013 Dynare Team
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
    
pfm.lead_lag_incidence = DynareModel.lead_lag_incidence;
pfm.ny = DynareModel.endo_nbr;
pfm.Sigma = DynareModel.Sigma_e;
pfm.Omega = chol(pfm.Sigma,'upper'); % Sigma = Omega'*Omega
pfm.number_of_shocks = length(pfm.Sigma);
pfm.stochastic_order = DynareOptions.ep.stochastic.order;
pfm.max_lag = DynareModel.maximum_endo_lag;
if pfm.max_lag > 0
    pfm.nyp = nnz(pfm.lead_lag_incidence(1,:));
    pfm.iyp = find(pfm.lead_lag_incidence(1,:)>0);
else
    pfm.nyp = 0;
    pfm.iyp = [];
end
pfm.ny0 = nnz(pfm.lead_lag_incidence(pfm.max_lag+1,:));
pfm.iy0 = find(pfm.lead_lag_incidence(pfm.max_lag+1,:)>0);
if DynareModel.maximum_endo_lead
    pfm.nyf = nnz(pfm.lead_lag_incidence(pfm.max_lag+2,:));
    pfm.iyf = find(pfm.lead_lag_incidence(pfm.max_lag+2,:)>0);
else
    pfm.nyf = 0;
    pfm.iyf = [];
end
pfm.nd = pfm.nyp+pfm.ny0+pfm.nyf;
pfm.nrc = pfm.nyf+1;
pfm.isp = [1:pfm.nyp];
pfm.is = [pfm.nyp+1:pfm.ny+pfm.nyp];
pfm.isf = pfm.iyf+pfm.nyp;
pfm.isf1 = [pfm.nyp+pfm.ny+1:pfm.nyf+pfm.nyp+pfm.ny+1];
pfm.iz = [1:pfm.ny+pfm.nyp+pfm.nyf];
pfm.periods = DynareOptions.ep.periods;
pfm.steady_state = DynareOutput.steady_state;
pfm.params = DynareModel.params;
if DynareModel.maximum_endo_lead
    pfm.i_cols_1 = nonzeros(pfm.lead_lag_incidence(pfm.max_lag+(1:2),:)');
    pfm.i_cols_A1 = find(pfm.lead_lag_incidence(pfm.max_lag+(1:2),:)');
else
    pfm.i_cols_1 = nonzeros(pfm.lead_lag_incidence(pfm.max_lag+1,:)');
    pfm.i_cols_A1 = find(pfm.lead_lag_incidence(pfm.max_lag+1,:)');
end
if pfm.max_lag > 0
    pfm.i_cols_T = nonzeros(pfm.lead_lag_incidence(1:2,:)');
else
    pfm.i_cols_T = nonzeros(pfm.lead_lag_incidence(1,:)');
end
pfm.i_cols_j = 1:pfm.nd;
pfm.i_upd = pfm.ny+(1:pfm.periods*pfm.ny);
pfm.dynamic_model = str2func([DynareModel.fname,'_dynamic']);
pfm.verbose = DynareOptions.ep.verbosity;
pfm.maxit_ = DynareOptions.simul.maxit;
pfm.tolerance = DynareOptions.dynatol.f;

if nargin>3 && DynareOptions.ep.stochastic.order
    % Compute weights and nodes for the stochastic version of the extended path.
    switch IntegrationMethod
      case 'Tensor-Gaussian-Quadrature'
        % Get the nodes and weights from a univariate Gauss-Hermite quadrature.
        [nodes,weights] = gauss_hermite_weights_and_nodes(DynareOptions.ep.stochastic.quadrature.nodes);
        % Replicate the univariate nodes for each innovation and dates, and, if needed, correlate them. 
        nodes = repmat(nodes,1,pfm.number_of_shocks*pfm.stochastic_order)*kron(eye(pfm.stochastic_order),pfm.Omega);
        % Put the nodes and weights in cells
        for i=1:pfm.number_of_shocks
            rr(i) = {nodes(:,i)};
            ww(i) = {weights};
        end
        % Build the tensorial grid
        pfm.nodes = cartesian_product_of_sets(rr{:});
        pfm.weights = prod(cartesian_product_of_sets(ww{:}),2);
        pfm.nnodes = length(pfm.weights);
      case 'Stroud-Cubature-3'
        [nodes,weights] = cubature_with_gaussian_weight(pfm.number_of_shocks*pfm.stochastic_order,3,'Stroud')
        pfm.nodes = kron(eye(pfm.stochastic_order),transpose(Omega))*nodes;
        pfm.weights = weights;
        pfm.nnodes = length(pfm.weights);
      case 'Stroud-Cubature-5'
        [nodes,weights] = cubature_with_gaussian_weight(pfm.number_of_shocks*pfm.stochastic_order,5,'Stroud')
        pfm.nodes = kron(eye(pfm.stochastic_order),transpose(Omega))*nodes;
        pfm.weights = weights;
        pfm.nnodes = length(weights);
      otherwise
        error('setup_stochastic_perfect_foresight_model_solver:: Unknown integration algorithm!')
    end
end
