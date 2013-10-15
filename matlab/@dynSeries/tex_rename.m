function ts = tex_rename(ts, varargin) % --*-- Unitary tests --*--

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

assert(nargin <= 3, 'dynSeries::tex_rename: accepts at most two args');
if nargin == 2
    newtex = varargin{1};
    assert(ts.vobs == 1, ['dynSeries::tex_rename: with one argument, the ' ...
                        'dynSeries contain only one variable.']);
else
    newtex = varargin{2};
    name = varargin{1};
    assert(ischar(name), 'dynSeries::tex_rename: name must be a string');
end
assert(ischar(newtex), 'dynSeries::tex_rename: the newtex name must be a string');

if nargin == 2
    idname = 1;
else
    idname = find(strcmp(name, ts.name));
    if isempty(idname)
        error(['dynSeries::tex_rename: Variable ' name ' is unknown in dynSeries object ' inputname(1)  '!'])
    end
end
ts.tex(idname) = {newtex};

%@test:1
%$ t = zeros(8,1);
%$ ts = dynSeries([transpose(1:5), transpose(6:10)],'1950q1',{'Output'; 'Consumption'}, {'Y_t'; 'C_t'});
%$ try
%$     ts = tex_rename(ts,'Output','\\Delta Y_t');
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dyn_assert(ts.freq,4);
%$     t(3) = dyn_assert(ts.init.freq,4);
%$     t(4) = dyn_assert(ts.init.time,[1950, 1]);
%$     t(5) = dyn_assert(ts.vobs,2);
%$     t(6) = dyn_assert(ts.nobs,5);
%$     t(7) = dyn_assert(ts.name,{'Output'; 'Consumption'});
%$     t(8) = dyn_assert(ts.tex,{'\\Delta Y_t'; 'C_t'});
%$ end
%$
%$ T = all(t);
%@eof:1

%@test:2
%$ t = zeros(8,1);
%$ ts = dynSeries([transpose(1:5), transpose(6:10)],'1950q1',{'Output'; 'Consumption'}, {'Y_t'; 'C_t'});
%$ try
%$     ts = ts.tex_rename('Output','\\Delta Y_t');
%$     t(1) = 1;
%$ catch
%$     t = 0;
%$ end
%$
%$ if length(t)>1
%$     t(2) = dyn_assert(ts.freq,4);
%$     t(3) = dyn_assert(ts.init.freq,4);
%$     t(4) = dyn_assert(ts.init.time,[1950, 1]);
%$     t(5) = dyn_assert(ts.vobs,2);
%$     t(6) = dyn_assert(ts.nobs,5);
%$     t(7) = dyn_assert(ts.name,{'Output'; 'Consumption'});
%$     t(8) = dyn_assert(ts.tex,{'\\Delta Y_t'; 'C_t'});
%$ end
%$
%$ T = all(t);
%@eof:2
