function simulate_debug(steady_state)
global M_ oo_ options_;
fid = fopen([M_.fname '_options.txt'],'wt');
if steady_state~=1
  fprintf(fid,'%d\n',options_.periods);
end;
fprintf(fid,'%d\n',options_.simul.maxit);
fprintf(fid,'%6.20f\n',options_.slowc);
fprintf(fid,'%6.20f\n',options_.markowitz);
fprintf(fid,'%6.20f\n',options_.dynatol.f);
fprintf(fid,'%d\n',options_.minimal_solving_periods);
fclose(fid);

fid = fopen([M_.fname '_M.txt'],'wt');
fprintf(fid,'%d\n',M_.maximum_lag);
fprintf(fid,'%d\n',M_.maximum_lead);
fprintf(fid,'%d\n',M_.maximum_endo_lag);
fprintf(fid,'%d\n',M_.param_nbr);
if steady_state==1
  fprintf(fid,'%d\n',size(oo_.exo_steady_state, 1));
  fprintf(fid,'%d\n',size(oo_.exo_steady_state, 2));
else
  fprintf(fid,'%d\n',size(oo_.exo_simul, 1));
  fprintf(fid,'%d\n',size(oo_.exo_simul, 2));
end;
fprintf(fid,'%d\n',M_.endo_nbr);
if steady_state==1
    fprintf(fid,'%d\n',size(oo_.steady_state, 2));
else
    fprintf(fid,'%d\n',size(oo_.endo_simul, 2));
end;
fprintf(fid,'%d\n',M_.exo_det_nbr);

fprintf(fid,'%d\n',size(oo_.steady_state,1));
fprintf(fid,'%d\n',size(oo_.steady_state,2));
fprintf(fid,'%d\n',size(oo_.exo_steady_state,1));
fprintf(fid,'%d\n',size(oo_.exo_steady_state,2));

fprintf(fid,'%6.20f\n',M_.params);

fclose(fid);
fid = fopen([M_.fname '_oo.txt'],'wt');
if steady_state==1
  fprintf(fid,'%6.20f\n',oo_.steady_state);
  fprintf(fid,'%6.20f\n',oo_.exo_steady_state);
else
  fprintf(fid,'%6.20f\n',oo_.endo_simul);
  fprintf(fid,'%6.20f\n',oo_.exo_simul);
end;
fprintf(fid,'%6.20f\n',oo_.steady_state);
fprintf(fid,'%6.20f\n',oo_.exo_steady_state);
fclose(fid);
