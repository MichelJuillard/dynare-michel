function [y_]=irf(dr, varexo, ssize, long_, drop_, replic, iorder)

if length(varexo) then
  shock_var = grep_exact(lgx_,varexo);
else
  shock_var = 1;
end

y_	= 0;
for j = 1: replic
  rand('seed',j,'n');
  ex1_ = rand(long_+drop_+xkmin_,exo_nbr,'n')*chol(Sigma_e_);
  ex2_ = ex1_;
  ex2_(drop_+1,shock_var) = ex2_(drop_+1,shock_var)+ssize;   
  y1_ = simult_(dr('ys'),dr,ex1_,iorder,long_+drop_);
  y2_ = simult_(dr('ys'),dr,ex2_,iorder,long_+drop_);
  y_ = y_+(y2_(:,ykmin_+drop_:$)-y1_(:,ykmin_+drop_:$));
end
y_ = y_/replic;
dsample(1,long_-drop_);

dyn2vec();




