
function rplot(s1)
		
% function rplot(s1)
% This optionnal command creates the plot of the variable trajectory. By default, the entire simulation period is
% ploted. The instruction DSAMPLE permits to reduce the number of periods in the plot.
% RPLOT :		RPLOT ( ['var1'; 'var2';  ...] )	
%
% INPUTS
%    s1:           list of vector names
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2001-2007)
% Gnu Public License.


global M_ oo_ options_

rplottype = options_.rplottype;

col = ['y','c','r','g','b','w','m'] ;
ix = [1 - M_.maximum_lag:size(oo_.endo_simul,2)-M_.maximum_lag]' ;

y = [];
for k=1:size(s1,1)
  if isempty(strmatch(s1(k,:),M_.endo_names,'exact'))
    error (['One of the variable specified does not exist']) ;
  end

  y = [y; oo_.endo_simul(strmatch(s1(k,:),M_.endo_names,'exact'),:)] ;
end

if options_.smpl == 0
        i = [M_.maximum_lag:size(oo_.endo_simul,2)]' ;
else
	i = [options_.smpl(1)+M_.maximum_lag:options_.smpl(2)+M_.maximum_lag]' ;
end

t = ['Plot of '] ;
if rplottype == 0
	for j = 1:size(y,1)
		t = [t s1(j,:) ' '] ;
	end
        figure ;
        plot(ix(i),y(:,i)) ;
	title (t,'Interpreter','none') ;
	xlabel('Periods') ;
	if size(s1,1) > 1
	  legend(s1,0);
	end
elseif rplottype == 1
	for j = 1:size(y,1)
		figure ;
		plot(ix(i),y(j,i)) ;
		title(['Plot of ' s1(:,j)]) ;
		xlabel('Periods') ;
	end
elseif rplottype == 2
	figure ;
	nl = max(1,fix(size(y,1)/4)) ;
	nc = ceil(size(y,1)/nl) ;
	for j = 1:size(y,1)
		subplot(nl,nc,j) ;
		plot(ix(i),y(j,i)) ;
		hold on ;
		plot(ix(i),oo_.steady_state(j)*ones(1,size(i,1)),'w:') ;
		xlabel('Periods') ;
		ylabel([s1(:,j)]) ;
		title(['Plot of ' s1(:,j)]) ;
	end
end

% 02/28/01 MJ replaced bseastr by MATLAB's strmatch
% 06/19/01 MJ added 'exact' to strmatch calls
% 06/25/03 MJ correction when options_.smpl ~= 0









