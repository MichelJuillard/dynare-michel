function formdata(fname,date)

% function formdata(fname,date)
% store endogenous and exogenous variables in a "FRM" TROLL text format file
% INPUT
%   fname: name of the FRM file
%   date:  the date of first observation (i.e. 2007A for an annual dataset)
% OUTPUT
%   none
% ALGORITHM
%   none
% SPECIAL REQUIREMENT
%   none
%    
% part of DYNARE, copyright (2007-2008)
% Gnu Public License.

  global M_ oo_
  fid = fopen([fname '_endo.frm'],'w');
  n=size(oo_.endo_simul,1);
  t=size(oo_.endo_simul,2);
  SN=upper(cellstr(M_.endo_names));
  for i=1:n
    str=strvcat(SN(i));
    fprintf(fid,'USER: x x DATAFILE: x %s\n',str);
    fprintf(fid,'PER: 1    YEAR: %s   FRAC: 1   NOBS: %d   CLINES: 0   DLINES: ???\n',date,t);
    fprintf(fid,'%10.5f %10.5f %10.5f %10.5f\n',reshape(oo_.endo_simul(i,1:floor(t/4)*4),floor(t/4),4));
    if(floor(t/4)*4<t)
        switch(t-floor(t/4)*4)
            case 1
              fprintf(fid,'%10.5f\n',oo_.endo_simul(i,floor(t/4)*4+1:t));
            case 2
              fprintf(fid,'%10.5f %10.5f\n',oo_.endo_simul(i,floor(t/4)*4+1:t));
            case 3
              fprintf(fid,'%10.5f %10.5f %10.5f\n',oo_.endo_simul(i,floor(t/4)*4+1:t));
        end;
    %else
    %    fprintf(fid,'\n');
    end;
  end;
  fclose(fid);
  
  fid = fopen([fname '_exo.frm'],'w');
  n=size(oo_.exo_simul,2);
  t=size(oo_.exo_simul,1);
  SN=upper(cellstr(M_.exo_names));
  for i=1:n
    str=strvcat(SN(i));
    fprintf(fid,'USER: x x DATAFILE: x %s\n',str);
    fprintf(fid,'PER: 1    YEAR: %s   FRAC: 1   NOBS: %d   CLINES: 0   DLINES: ???\n',date,t);
    fprintf(fid,'%10.5f %10.5f %10.5f %10.5f\n',reshape(oo_.exo_simul(1:floor(t/4)*4,i),floor(t/4),4));
    if(floor(t/4)*4<t)
        switch(t-floor(t/4)*4)
            case 1
              fprintf(fid,'%10.5f\n',oo_.exo_simul(floor(t/4)*4+1:t,i)');
            case 2
              fprintf(fid,'%10.5f %10.5f\n',oo_.exo_simul(floor(t/4)*4+1:t,i)');
            case 3
              fprintf(fid,'%10.5f %10.5f %10.5f\n',oo_.exo_simul(floor(t/4)*4+1:t,i)');
        end;
    %else
    %    fprintf(fid,'\n');
    end;
  end;
  fclose(fid);
return;