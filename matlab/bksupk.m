function d1 = bksupk(ny,fid,jcf,icc1)

global M_ options_

icf = [1:jcf-1] ;
ir = [(options_.periods-1)*ny+1:ny*options_.periods] ;
irf = icc1+(options_.periods-1)*ny ;
d1 = zeros(options_.periods*ny,1) ;

ofs = (((options_.periods-1)*ny+1)-1)*jcf*8 ;
junk = fseek(fid,ofs,-1) ;
c = fread(fid,[jcf,ny],'float64') ;

d1(ir) = c(:,jcf) ;
ir = ir-ny ;

i = 2 ;

while i <= M_.maximum_lead | i <= options_.periods
	irf1 = selif(irf,irf<=options_.periods*ny) ;

	ofs = (((options_.periods-i)*ny+1)-1)*jcf*8 ;
	junk = fseek(fid,ofs,-1) ;
	c = fread(fid,[jcf,ny],'float64') ;

	d1(ir) = c(:,jcf) - c(:,1:size(irf1,1))*d1(irf1) ;
	ir = ir - ny ;
	irf = irf - ny ;
	i = i + 1 ;
end

while i <= options_.periods

	ofs = (((options_.periods-i)*ny+1)-1)*jcf*8 ;
	junk = fseek(fid,ofs,-1) ;
	c = fread(fid,[jcf,ny],'float64') ;

	d1(ir) = c(:,jcf)-c(:,icf)*d1(irf) ;
	ir = ir-ny ;			
	irf = irf-ny ;
	i = i+1;
end

return ;
