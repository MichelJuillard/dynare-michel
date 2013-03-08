function A = subsasgn(A,S,B)

%@info:
%! @deftypefn {Function File} {@var{A} =} subsasgn (@var{A}, @var{S}, @var{B})
%! @anchor{@dynSeries/subsasgn}
%! @sp 1
%! Overloads the subsasgn method for the Dynare time series class (@ref{dynSeries}).
%! @end deftypefn
%@eod:

% Copyright (C) 2012 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT fr

if isa(A,'dynSeries') && isa(B,'dynSeries') 
    if length(S)==1 && isequal(S.type,'{}')
        if isequal(A.nobs,B.nobs) && isequal(A.init,B.init)
            id = NaN(length(S.subs),1);
            for i=1:length(S.subs)
                tmp = strmatch(S.subs{i},A.name,'exact');
                if isempty(tmp)
                    error(['dynSeries::subsasgn: variable ' S.subs{i} ' is not a member of ' inputname(1) ' dynSeries object!'])
                else
                    id(i) = tmp;
                end 
            end
            if isequal(B.vobs,length(S.subs))
                A.name(id) = B.name;
                A.data(:,id) = B.data;
                return
            end
        end
        return
    elseif length(S)==1 && isequal(S.type,'.')
        A = merge(A,B);
        return
    end
end

error('dynSeries::subsasgn: Wrong calling sequence!')
    
%@test:1
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(10,1);
%$
%$ % Instantiate two dynSeries object.
%$ ts1 = dynSeries(A,[],{'A1';'A2';'A3'},[]);
%$ ts2 = dynSeries(B,[],{'B1'},[]);
%$
%$ % modify first object.
%$ ts1{'A2'} = ts2;
%$ t(1) = 1;
%$ % Instantiate a time series object.
%$
%$    t(2) = dyn_assert(ts1.vobs,3);
%$    t(3) = dyn_assert(ts1.nobs,10);
%$    t(4) = dyn_assert(ts1.name{2},'B1');
%$    t(5) = dyn_assert(ts1.name{1},'A1');
%$    t(6) = dyn_assert(ts1.name{3},'A3');
%$    t(7) = dyn_assert(ts1.data,[A(:,1), B, A(:,3)],1e-15);
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define a datasets.
%$ A = rand(10,3);
%$
%$ % Instantiate two dynSeries object.
%$ ts1 = dynSeries(A,[],{'A1';'A2';'A3'},[]);
%$
%$ % Apply the exponential function to the second variable.
%$ ts1{'A2'} = ts1{'A2'}.exp;
%$
%$ % Instantiate a time series object.
%$
%$    t(1) = dyn_assert(ts1.vobs,3);
%$    t(2) = dyn_assert(ts1.nobs,10);
%$    t(3) = dyn_assert(ts1.name{2},'A2');
%$    t(4) = dyn_assert(ts1.name{1},'A1');
%$    t(5) = dyn_assert(ts1.name{3},'A3');
%$    t(6) = dyn_assert(ts1.data,[A(:,1), exp(A(:,2)), A(:,3)],1e-15);
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define a datasets.
%$ A = rand(10,3);
%$
%$ % Instantiate two dynSeries object.
%$ ts1 = dynSeries(A,[],{'A1';'A2';'A3'},[]);
%$
%$ % Apply the logarithm function to the first and third variables.
%$ ts1{'A1'} = ts1{'A1'}.log;
%$ ts1{'A3'} = ts1{'A3'}.log;
%$
%$ % Instantiate a time series object.
%$
%$    t(1) = dyn_assert(ts1.vobs,3);
%$    t(2) = dyn_assert(ts1.nobs,10);
%$    t(3) = dyn_assert(ts1.name{2},'A2');
%$    t(4) = dyn_assert(ts1.name{1},'A1');
%$    t(5) = dyn_assert(ts1.name{3},'A3');
%$    t(6) = dyn_assert(ts1.data,[log(A(:,1)), A(:,2), log(A(:,3))],1e-15);
%$ T = all(t);
%@eof:3

%@test:4
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(10,3);
%$
%$ % Instantiate two dynSeries object.
%$ ts1 = dynSeries(A,[],{'A1';'A2';'A3'},[]);
%$ ts2 = dynSeries(B,[],{'B1';'B2';'B3'},[]);
%$
%$ % Apply the logarithm function to the first and third variables.
%$ ts1.A1 = ts2.B1;
%$
%$ % Instantiate a time series object.
%$
%$ t(1) = dyn_assert(ts1.vobs,4);
%$ t(2) = dyn_assert(ts1.nobs,10);
%$ t(3) = dyn_assert(ts1.name{1},'A1');
%$ t(3) = dyn_assert(ts1.name{4},'B1');
%$ T = all(t);
%@eof:4

%@test:5
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(10,3);
%$
%$ % Instantiate two dynSeries object.
%$ ts1 = dynSeries(A,[],{'A1';'A2';'A3'},[]);
%$ ts2 = dynSeries(B,[],{'A1';'B2';'B3'},[]);
%$
%$ % Apply the logarithm function to the first and third variables.
%$ ts1.A1 = ts2.A1;
%$
%$ % Instantiate a time series object.
%$
%$ t(1) = dyn_assert(ts1.vobs,3);
%$ t(2) = dyn_assert(ts1.nobs,10);
%$ t(3) = dyn_assert(ts1.name{1},'A1');
%$ t(4) = dyn_assert(ts1.data(:,1),B(:,1), 1e-15);
%$ t(5) = dyn_assert(ts1.data(:,2:3),A(:,2:3), 1e-15);
%$ T = all(t);
%@eof:5

