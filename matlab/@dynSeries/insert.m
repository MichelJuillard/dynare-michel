function ts = insert(ts,us,id) % --*-- Unitary tests --*--

% Add a variable in a dynSeries object.

%@info:
%! @deftypefn {Function File} {@var{ts} =} insert (@var{ts}, @var{us}, @var{id})
%! @anchor{dynSeries/insert}
%! @sp 1
%! Insert method for the dynSeries class. Insert new variables (@var{us}) in a dynSeries object (@var{ts}) at
%! positions @var{id}.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item ts
%! Object instantiated by @ref{dynSeries}.
%! @item us
%! Object instantiated by @ref{dynSeries}.
%! @item id
%! vector of integers.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item ts
%! Object instantiated by @ref{dynSeries}, without variable (@var{a}).
%! @end table
%! @end deftypefn
%@eod:

% Copyright (C) 2013 Dynare Team
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

[n,message] = common_strings_in_cell_arrays(ts.name,us.name);

if n
    error(['dynSeries::insert: Variable(s) ' message ' already exist in ''' inputname(1) '''!'])
end

if ~isequal(ts.freq,us.freq)
    error(['dynSeries::insert: ''' inputname(1) ''' and ''' inputname(2) ''' dynSeries objects must have common frequencies!'])
end

[ts,us] = align(ts, us);

n = length(id);

for i=1:n
    ts.data = insert_column_vector_in_a_matrix(ts.data,us.data(:,i),id(i));
    ts.name = insert_object_in_a_one_dimensional_cell_array(ts.name,us.name{i},id(i));
    ts.tex = insert_object_in_a_one_dimensional_cell_array(ts.tex,us.tex{i},id(i));
    id = id+1;
end

%@test:1
%$ % Define a datasets.
%$ A = rand(10,3); B = rand(5,2);
%$
%$ % Define names.
%$ A_name = {'A1'; 'A2';'A3'};
%$ B_name = {'B1'; 'B2'};
%$
%$ % Define initial dates.
%$ A_init = '1950Q1';
%$ B_init = '1950Q3';
%$
%$ % Instantiate two dynSeries objects.
%$ ts1 = dynSeries(A, A_init, A_name,[]);
%$ ts2 = dynSeries(B, B_init, B_name,[]);
%$
%$ try
%$    ts1 = insert(ts1,ts2,[1,2]);
%$    t(1) = 1;
%$ catch
%$    t = 0;
%$ end
%$
%$ if length(t)>1
%$    t(2) = dyn_assert(ts1.vobs,{'B1';'A1';'B2';'A3'});
%$    t(3) = dyn_assert(ts1.nobs,10);
%$    eB = [NaN(2,2); B; NaN(3,2)];
%$    t(4) = dyn_assert(ts1.data,[eB(:,1), A(:,1), eB(:,2), A(:,2:3)], 1e-15);
%$ end
%$ T = all(t);
%@eof:1