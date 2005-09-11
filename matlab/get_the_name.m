function [nam,texnam] = get_the_name(k,TeX)
% stephane.adjemian@cepremap.cnrs.fr [07-13-2004]
global M_ estim_params_ options_

nam = [];
texnam = [];

if k <= estim_params_.nvx
    vname = deblank(M_.exo_names(estim_params_.var_exo(k,1),:));
    nam=['SE_',vname];
    if TeX
        tname  = deblank(M_.exo_names_tex(estim_params_.var_exo(k,1),:));
        texnam = ['$ SE_{' tname '} $'];
    end
elseif  k <= (estim_params_.nvx+estim_params_.nvn)
    vname = deblank(options_.varobs(estim_params_.var_endo(k-estim_params_.nvx,1),:));
    nam=['SE_EOBS_',vname];
    if TeX
        tname  = deblank(options_.TeX_varobs(estim_params_.var_endo(k-estim_params_.nvx,1),:));
        texnam = ['$ SE_{' tname '} $'];
    end
elseif  k <= (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx)
    jj = k - (estim_params_.nvx+estim_params_.nvn);
    k1 = estim_params_.corrx(jj,1);
    k2 = estim_params_.corrx(jj,2);
    vname = [deblank(M_.exo_names(k1,:)) '_' deblank(M_.exo_names(k2,:))];
    nam=['CC_',vname];
    if TeX
        tname  = [deblank(M_.exo_names_tex(k1,:)) ',' deblank(M_.exo_names_tex(k2,:))];
        texnam = ['$ CC_{' tname '} $'];
    end
elseif  k <= (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+ ...
		            estim_params_.ncn)
    jj = k - (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx);
    k1 = estim_params_.corrn(jj,1);
    k2 = estim_params_.corrn(jj,2);
    vname = [deblank(M_.endo_names(k1,:)) '_' deblank(M_.endo_names(k2,:))];
    nam=['CC_EOBS_' vname];
    if TeX
        tname  = [deblank(M_.endo_names_tex(k1,:)) ',' deblank(M_.endo_names_tex(k2,:))];
        texnam =['$ CC_{' tname '} $'];
    end
else
    jj = k - (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+ ...
	      estim_params_.ncn); 
    jj1 = estim_params_.param_vals(jj,1);
    nam = deblank(M_.param_names(jj1,:));
    if TeX
        texnam = ['$ '  deblank(M_.param_names_tex(jj1,:))  ' $'];
    end    
end


% SA 07-15-2004 Added TeX names.
% SA 12-02-2004 Changed non-TeX names format.
% SA 01-11-2005 v3TOv4