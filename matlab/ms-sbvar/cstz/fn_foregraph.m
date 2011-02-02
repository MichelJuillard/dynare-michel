function fn_foregraph(yfore,yacte,keyindx,rnum,cnum,q_m,ylab,forelabel,conlab)
%
%   Graph annual (or calendar year) forecasts vs actual (all data from "msstart.m")
%
% yfore:  actual and forecast annual growth data with dates.
% yacte:  actual annual growth data with dates.
% keyindx:  index for the variables to be graphed
% rnum:  number of rows in subplot
% cnum:  number of columns in subplot
% q_m:  if 4 or 12, quarterly or monthly data
% ylab:  string array for the length(keyindx)-by-1 variables
% forelabel:  title label for as of time of forecast
% conlab:  label for what conditions imposed; e.g., conlab = 'All bar MS shocks inspl'
%-------------
% No output argument for this graph file
%  See fn_seriesgraph.m, fn_forerrgraph.m.
%
% Tao Zha, March 2000


vyrs = yfore(:,1);
hornum = cell(length(vyrs),1);    % horizontal year (number)
count=0;
for k=vyrs'
   count=count+1;
   jnk=num2str(k);
   hornum{count}=jnk(3:4);   % e.g., with '1990', we have '90'
end

count=0;
for i = keyindx
   count = count+1;
   subplot(rnum,cnum,count)
   plot(yacte(:,1)+yacte(:,2)/q_m,yacte(:,2+i),yfore(:,1)+yfore(:,2)/q_m,yfore(:,2+i),'--')

   if (yfore(1,2)==0)   % only for annual growth rates (not for, say, monthly annualized rates)
      set(gca,'XLim',[vyrs(1) vyrs(end)])
      set(gca,'XTick',vyrs)
      set(gca,'XTickLabel',char(hornum))
   end

   if i==keyindx(1)
      title(forelabel)
   elseif i>=length(keyindx)  %i>=length(keyindx)-1
      xlabel(conlab)
   end
   %
   grid
   ylabel(char(ylab(i)))
end
