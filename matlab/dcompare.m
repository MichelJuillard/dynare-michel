% Copyright (C) 2001 Michel Juillard
%
function dcompare(s1)

global options_ nvx nvy x y lag1

ftest(s1,0) ;

i = [lag1(1):size(x,2)-lag1(2)+1]' ;

if size(options_.smpl,1) == 1
	error(['DSAMPLE not specified.']) ;
end

if options_.smpl(3) > 0
	if options_.smpl(3) == 2
		if options_.smpl(1)<0 | options_.smpl(2)>size(x,2)-lag1(2)
			error ('Wrong sample.') ;
		end
		i = [options_.smpl(1)+lag1(1):options_.smpl(2)+lag1(1)]' ;
	elseif options_.smpl(3) == 1
		if options_.smpl(1)>size(x,2)-lag1(2)
			error ('Wrong sample.') ;
		end
		i = [lag1(1):options_.smpl(1)+lag1(1)]' ;
	end
end

j = bseastr(nvx,nvy) ;

if stop
	return ;
end

z = mean(mean(abs(x(j,i)-y(j,i)))) ;

disp (['The mean absolute difference between set ' s1(1,:) 'and set ' s1(2,:)]) ;
disp (['is : ' num2str(z)]) ;
return ;


