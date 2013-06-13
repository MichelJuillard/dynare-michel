function A = subsasgn(A,S,B)

%@info:
%! @deftypefn {Function File} {@var{A} =} subsasgn (@var{A}, @var{S}, @var{B})
%! @anchor{@dynSeries/subsasgn}
%! @sp 1
%! Overloads the subsasgn method for the Dynare time series class (@ref{dynSeries}).
%! @end deftypefn
%@eod:

% Copyright (C) 2012-2013 Dynare Team
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

if length(S)>1
    error('dynSeries::subsasgn: Wrong syntax!')
end

switch S.type
  case '{}'
    if ~isequal(numel(S.subs),numel(unique(S.subs)))
        error('dynSeries::subsasgn: Wrong syntax!')
    end
    for i=1:numel(S.subs)
        element = S.subs{i};
        idArobase = strfind(element,'@');
        if ~isempty(idArobase)
            switch length(idArobase)
              case 2
                idComma = strfind(element(idArobase(1)+1:idArobase(2)-1),',');
                if ~isempty(idComma)
                    elements = cell(1,numel(idComma)+1); j = 1;
                    expression = element(idArobase(1)+1:idArobase(2)-1);
                    while ~isempty(expression)
                        [token, expression] = strtok(expression,',');
                        elements(j) = {[element(1:idArobase(1)-1), token, element(idArobase(2)+1:end)]};
                        j = j + 1;
                    end
                    S.subs = replace_object_in_a_one_dimensional_cell_array(S.subs, elements(:), i);
                else
                    error('dynSeries::subsasgn: Wrong syntax, matlab''s regular expressions cannot be used here!')
                end
              case 4
                idComma_1 = strfind(element(idArobase(1)+1:idArobase(2)-1),',');
                idComma_2 = strfind(element(idArobase(3)+1:idArobase(4)-1),',');
                if ~isempty(idComma_1)
                    elements = cell(1,(numel(idComma_1)+1)*(numel(idComma_2)+1)); j = 1;
                    expression_1 = element(idArobase(1)+1:idArobase(2)-1);
                    while ~isempty(expression_1)
                        [token_1, expression_1] = strtok(expression_1,',');
                        expression_2 = element(idArobase(3)+1:idArobase(4)-1);
                        while ~isempty(expression_2)
                            [token_2, expression_2] = strtok(expression_2,',');
                            elements(j) = {[element(1:idArobase(1)-1), token_1, element(idArobase(2)+1:idArobase(3)-1), token_2, element(idArobase(4)+1:end)]};
                            j = j+1;
                        end
                    end
                    S.subs = replace_object_in_a_one_dimensional_cell_array(S.subs, elements(:), i);
                else
                    error('dynSeries::subsasgn: Wrong syntax, matlab''s regular expressions cannot be used here!')
                end
              otherwise
                error('dynSeries::subsasgn: Wrong syntax!')
            end
        end
    end
    if ~isequal(length(S.subs),B.vobs)
        error('dynSeries::subsasgn: Wrong syntax!')
    end
    if ~isequal(S.subs(:),B.name)
        for i = 1:B.vobs
            if ~isequal(S.subs{i},B.name{i})
                % Rename a variable.
                id = strmatch(S.subs{i},A.name);
                if isempty(id)
                    % Add a new variable a change its name.
                    B.name(i) = {S.subs{i}};
                    B.tex(i) = {name2tex(S.subs{i})};
                else
                    % Rename variable and change its content.
                    B.name(i) = A.name(id);
                    B.tex(i) = A.tex(id);
                end
            end
        end
    end
  case '.'
    if ~isequal(S.subs,B.name)
        if ~isequal(S.subs,B.name{1})
                % Rename a variable.
                id = strmatch(S.subs,A.name);
                if isempty(id)
                    % Add a new variable a change its name.
                    B.name(1) = {S.subs};
                    B.tex(1) = {name2tex(S.subs)};
                else
                    % Rename variable and change its content.
                    B.name(1) = A.name(id);
                    B.tex(1) = A.tex(id);
                end
            end
        end
  otherwise
    error('dynSeries::subsasgn: Wrong syntax!')
end
  
A = merge(A,B);

%@test:1
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(10,1);
%$
%$ % Instantiate two dynSeries object.
%$ ts1 = dynSeries(A,[],{'A1';'A2';'A3'},[]);
%$ ts2 = dynSeries(B,[],{'B1'},[]);
%$
%$ % modify first object.
%$ try
%$     ts1{'A2'} = ts2;
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end 
%$ 
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dyn_assert(ts1.vobs,3);
%$    t(3) = dyn_assert(ts1.nobs,10);
%$    t(4) = dyn_assert(ts1.name{2},'A2');
%$    t(5) = dyn_assert(ts1.name{1},'A1');
%$    t(6) = dyn_assert(ts1.name{3},'A3');
%$    t(7) = dyn_assert(ts1.data,[A(:,1), B, A(:,3)],1e-15);
%$ end
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
%$ A = rand(10,3);
%$
%$ % Instantiate two dynSeries object.
%$ ts1 = dynSeries(A,[],{'A1';'A2';'A3'},[]);
%$
%$ % Apply the logarithm function to the first and third variables of ts1.
%$ ts1{'A1','A3'} = ts1{'A1','A3'}.log;
%$
%$ t(1) = dyn_assert(ts1.vobs,3);
%$ t(2) = dyn_assert(ts1.nobs,10);
%$ t(3) = dyn_assert(ts1.name{1},'A1') && dyn_assert(ts1.name{2},'A2') && dyn_assert(ts1.name{3},'A3');
%$ t(4) = dyn_assert(ts1.data,[log(A(:,1)), A(:,2), log(A(:,3))],1e-15);
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

%@test:6
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(10,2);
%$
%$ % Instantiate two dynSeries object.
%$ ts1 = dynSeries(A,[],{'A1';'A2';'A3'},[]);
%$ ts2 = dynSeries(B,[],{'B1';'B2'},[]);
%$
%$ % Call tested routine.
%$ try
%$     ts1.B2 = ts2.B2.log;
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dyn_assert(ts1.vobs,4);
%$    t(3) = dyn_assert(ts1.nobs,10);
%$    t(4) = dyn_assert(ts1.name{1},'A1');
%$    t(5) = dyn_assert(ts1.name{2},'A2');
%$    t(6) = dyn_assert(ts1.name{3},'A3');
%$    t(7) = dyn_assert(ts1.name{4},'B2');
%$    t(8) = dyn_assert(ts1.data,[A(:,1), A(:,2), A(:,3), log(B(:,2))],1e-15);
%$ end
%$ T = all(t);
%@eof:6

%@test:7
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(10,2);
%$
%$ % Instantiate two dynSeries object.
%$ ts1 = dynSeries(A,[],{'A1';'A2';'A3'},[]);
%$ ts2 = dynSeries(B,[],{'B1';'B2'},[]);
%$
%$ % Append B2 to the first object.
%$ ts1{'B2'} = ts2{'B2'};
%$ t(1) = 1;
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dyn_assert(ts1.vobs,4);
%$    t(3) = dyn_assert(ts1.nobs,10);
%$    t(4) = dyn_assert(ts1.name{1},'A1');
%$    t(5) = dyn_assert(ts1.name{2},'A2');
%$    t(6) = dyn_assert(ts1.name{3},'A3');
%$    t(6) = dyn_assert(ts1.name{4},'B2');
%$    t(7) = dyn_assert(ts1.data,[A(:,1), A(:,2), A(:,3), B(:,2)],1e-15);
%$ end
%$ T = all(t);
%@eof:7

%@test:8
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(10,1);
%$
%$ % Instantiate two dynSeries object.
%$ ts1 = dynSeries(A,[],{'A1';'A2';'A3'},[]);
%$ ts2 = dynSeries(B,[],{'B1'},[]);
%$
%$ % modify first object.
%$ try
%$     ts1{'A4'} = ts2;
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end 
%$ 
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dyn_assert(ts1.vobs,4);
%$    t(3) = dyn_assert(ts1.nobs,10);
%$    t(4) = dyn_assert(ts1.name{2},'A2');
%$    t(5) = dyn_assert(ts1.name{1},'A1');
%$    t(6) = dyn_assert(ts1.name{3},'A3');
%$    t(7) = dyn_assert(ts1.name{4},'A4');
%$    t(8) = dyn_assert(ts1.data,[A, B],1e-15);
%$ end
%$ T = all(t);
%@eof:8

%@test:9
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(10,2);
%$
%$ % Instantiate two dynSeries object.
%$ ts1 = dynSeries(A,[],{'A1';'A2';'A3'},[]);
%$ ts2 = dynSeries(B,[],{'A1';'B1'},[]);
%$
%$ % modify first object.
%$ try
%$     ts1{'A1','A4'} = ts2;
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end 
%$ 
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dyn_assert(ts1.vobs,4);
%$    t(3) = dyn_assert(ts1.nobs,10);
%$    t(4) = dyn_assert(ts1.name{2},'A2');
%$    t(5) = dyn_assert(ts1.name{1},'A1');
%$    t(6) = dyn_assert(ts1.name{3},'A3');
%$    t(7) = dyn_assert(ts1.name{4},'A4');
%$    t(8) = dyn_assert(ts1.data,[B(:,1), A(:,2:3), B(:,2)],1e-15);
%$ end
%$ T = all(t);
%@eof:9


%@test:10
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(10,3);
%$
%$ % Instantiate two dynSeries object.
%$ ts1 = dynSeries(A,[],{'A1';'A2';'A3'},[]);
%$ ts2 = dynSeries(B,[],{'A1';'B1';'B2'},[]);
%$
%$ % modify first object.
%$ try
%$     ts1{'A@1,2@','A4'} = ts2;
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end 
%$ 
%$ % Instantiate a time series object.
%$ if t(1)
%$    t(2) = dyn_assert(ts1.vobs,4);
%$    t(3) = dyn_assert(ts1.nobs,10);
%$    t(4) = dyn_assert(ts1.name{1},'A1');
%$    t(5) = dyn_assert(ts1.name{2},'A2');
%$    t(6) = dyn_assert(ts1.name{3},'A3');
%$    t(7) = dyn_assert(ts1.name{4},'A4');
%$    t(8) = dyn_assert(ts1.data,[B(:,1:2), A(:,3), B(:,3)],1e-15);
%$ end
%$ T = all(t);
%@eof:10

%@test:10
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(10,5);
%$
%$ % Instantiate two dynSeries object.
%$ ts1 = dynSeries(A,[],{'A_1';'A_2';'A_3'},[]);
%$ ts2 = dynSeries(B,[],{'A_1';'A_2';'B_1';'B_2';'B_3'},[]);
%$
%$ % modify first object.
%$ try
%$     ts1{'@A,B@_@1,2@'} = ts2{'@A,B@_@1,2@'};
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$
%$ % Instantiate a time series object.
%$ if t(1)
%$    %t(2) = dyn_assert(ts1.vobs,4);
%$    %t(3) = dyn_assert(ts1.nobs,10);
%$    %t(4) = dyn_assert(ts1.name,{'A1','A2';'A3';'B1';'B2'});
%$    %t(5) = dyn_assert(ts1.data,[B(:,1:2), A(:,3), B(:,3:4)],1e-15);
%$ end
%$ T = all(t);
%@eof:10
