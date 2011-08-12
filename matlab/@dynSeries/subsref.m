function B = subsref(ts, S)
    if isequal(S.type,'.')
        switch S.subs
          case {'data','nobs','vobs','name','tex','freq','time','init','last'} % Public members.
            B = builtin('subsref', ts, S);
          case {'log','exp'}                                                   % Give "dot access" to public methods.
            B = feval(S.subs,ts);
          otherwise                                                            % Extract a sub-object by selecting one variable.
            ndx = strmatch(S.subs,ts.name);
            if ~isempty(ndx)
                B = dynSeries();
                B.data = ts.data(:,ndx);
                B.name = deblank(ts.name(ndx,:));
                B.tex  = deblank(ts.tex(ndx,:));
                B.nobs = ts.nobs;
                B.vobs = 1;
                B.freq = ts.freq;
                B.time = ts.time;
                B.init = ts.init;
                B.last = ts.last;
                return
            else
                error('dynSeries::subsref: Unknown public method, public member or variable!')
            end
        end
        return
    end
    if isequal(S.type,'()')                                                    % Extract a sub-object by selecting a sub-sample.
        B = dynSeries();
        if size(ts.data,2)>1
            S.subs = [S.subs, ':'];
        end
        B.data = builtin('subsref', ts.data, S);
        B.nobs = size(B.data,1);
        B.vobs = ts.vobs;
        B.freq = ts.freq;
        B.time = builtin('subsref', ts.time, S);
        B.init = B.time(1,:);
        B.last = B.time(end,:);
        B.name = ts.name;
        B.tex  = ts.tex;
    end

%@test:1
%$ addpath ../matlab
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = char('A1','A2');
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ a = ts1(2:9);
%$
%$ % Expected results.
%$ e.data = [transpose(2:9),2*transpose(2:9)];
%$ e.nobs = 8;
%$ e.vobs = 2;
%$ e.name = char('A1','A2');
%$ e.freq = 1;
%$ tmp = ts1.time; e.time = tmp(2:9,:);
%$ e.init = e.time(1,:);
%$ e.last = e.time(end,:);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(a.data,e.data);
%$ t(2) = dyn_assert(a.time,e.time);
%$ t(3) = dyn_assert(a.nobs,e.nobs);
%$ t(4) = dyn_assert(a.vobs,e.vobs);
%$ t(5) = dyn_assert(a.freq,e.freq);
%$ t(6) = dyn_assert(a.init,e.init);
%$ t(7) = dyn_assert(a.last,e.last);
%$ T = all(t);
%@eof:1

%@test:2
%$ addpath ../matlab
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = char('A1','A2');
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ a = ts1.A1;
%$
%$ % Expected results.
%$ e.data = transpose(1:10);
%$ e.nobs = 10;
%$ e.vobs = 1;
%$ e.name = char('A1');
%$ e.freq = 1;
%$ e.time = [transpose(1:10),ones(10,1)];
%$ e.init = e.time(1,:);
%$ e.last = e.time(end,:);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(a.data,e.data);
%$ t(2) = dyn_assert(a.time,e.time);
%$ t(3) = dyn_assert(a.nobs,e.nobs);
%$ t(4) = dyn_assert(a.vobs,e.vobs);
%$ t(5) = dyn_assert(a.freq,e.freq);
%$ t(6) = dyn_assert(a.init,e.init);
%$ t(7) = dyn_assert(a.last,e.last);
%$ T = all(t);
%@eof:2

%@test:3
%$ addpath ../matlab
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = char('A1','A2');
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ a = ts1.log;
%$
%$ % Expected results.
%$ e.data = log(A);
%$ e.nobs = 10;
%$ e.vobs = 2;
%$ e.name = char('A1','A2');
%$ e.freq = 1;
%$ tmp = ts1.time; e.time = tmp(1:10,:);
%$ e.init = e.time(1,:);
%$ e.last = e.time(end,:);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(a.data,e.data);
%$ t(2) = dyn_assert(a.time,e.time);
%$ t(3) = dyn_assert(a.nobs,e.nobs);
%$ t(4) = dyn_assert(a.vobs,e.vobs);
%$ t(5) = dyn_assert(a.freq,e.freq);
%$ t(6) = dyn_assert(a.init,e.init);
%$ t(7) = dyn_assert(a.last,e.last);
%$ T = all(t);
%@eof:3
