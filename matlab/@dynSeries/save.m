function save(A,basename,format)
% Saves a dynSeries object on disk.
    
if nargin<3 || isempty(format)
    format = 'csv';
end

if nargin<2 || isempty(basename)
    basename = inputname(1);
end

switch format
  case 'm'
    fid = fopen([basename, '.m'],'w');
    fprintf(fid,'%% File created on %s.\n',datestr(now));
    fprintf(fid,'\n');
    fprintf(fid,'FREQ__ = %s;\n',num2str(A.freq));
    fprintf(fid,'INIT__ = '' %s'';\n',A.init.format);
    fprintf(fid,'\n');
    fprintf(fid,'NAMES__ = {');
    for i=1:A.vobs
        fprintf(fid,[ '''' A.name{i}  '''']);
        if i<A.vobs
            fprintf(fid,'; ');
        end
    end
    fprintf(fid,'};\n');
    fprintf(fid,'TEX__ = {');
    for i=1:A.vobs
        fprintf(fid,['''' A.tex{i}  '''']);
        if i<A.vobs
            fprintf(fid,'; ');
        end
    end
    fprintf(fid,'};\n');
    for v=1:A.vobs
        fprintf(fid,'%s = [\n', A.name{v});
        fprintf(fid,'%15.8g\n',A.data(1:end-1,v));
        fprintf(fid,'%15.8g];\n\n',A.data(end,v));
    end
    fclose(fid);
  case 'mat'
    FREQ__ = A.freq;
    INIT__ = A.init.format;
    NAMES__ = A.name;
    TEX__ = A.tex;
    str = [];
    for v=1:A.vobs
        str = [str, A.name{v} ' = A.data(:,' num2str(v) ');' ];
    end
    eval(str);
    save([basename '.mat'],'INIT__','FREQ__','NAMES__','TEX__',A.name{:});
  case 'csv'
    fid = fopen([basename, '.csv'],'w');
    fprintf(fid,', %s', A.name{:});
    fprintf(fid,'\n');
    for t=1:A.nobs
        date = A.init+(t-1);
        str = sprintf(', %15.8g',A.data(t,:));
        fprintf(fid, '%s%s\n',date.format,str);
    end
    fclose(fid);
end

%@test:1
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dynSeries(A,[],A_name,[]);
%$    save(ts1,[],'csv');
%$    t = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dynSeries(A,[],A_name,[]);
%$    save(ts1,[],'m');
%$    t = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dynSeries(A,[],A_name,[]);
%$    save(ts1,[],'mat');
%$    t = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ T = all(t);
%@eof:3

%@test:4
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$
%$ % Instantiate a time series object.
%$ try
%$    ts1 = dynSeries(A,[],A_name,[]);
%$    ts1.save;
%$    t = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ T = all(t);
%@eof:4