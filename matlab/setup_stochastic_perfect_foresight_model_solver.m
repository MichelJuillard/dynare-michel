function pfm = setup_stochastic_perfect_foresight_model_solver(DynareModel,DynareOptions,DynareOutput,Algorithm,IntegrationMethod)
    
pfm.lead_lag_incidence = DynareModel.lead_lag_incidence;
pfm.ny = DynareModel.endo_nbr;
pfm.Sigma_e = DynareModel.Sigma_e;
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
pfm.maxit_ = DynareOptions.maxit_;
pfm.tolerance = DynareOptions.dynatol.f;
