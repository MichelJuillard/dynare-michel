function [freq,init,data,varlist,tex] = load_m_file_data(file)

%@info:
%! @deftypefn {Function File} {@var{freq}, @var{init}, @var{data}, @var{varlist} =} load_m_file_data (@var{file})
%! @anchor{load_m_file_data}
%! @sp 1
%! Loads data in a matlab/octave m-file.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item file
%! string, name of the m file (matlab/octave script).
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item freq
%! Scalar integer (1, 4, 12, 52).
%! @item init
%! dynDate object, initial date.
%! @item data
%! Matrix of doubles, data.
%! @item varlist
%! Cell of strings (names of the variables in the database).
%! @end table
%! @sp 2
%! @strong{Remarks}
%! @sp 1
%! The frequency and initial date can be specified with variables FREQ__ and INIT__ in the matlab/octave script. FREQ__ must be a scalar integer and INIT__ a string like '1938M11', '1945Q3', '1973W3' or '2009'. If these variables are not specified default values for freq and init are 1 and dynDate(1).
%! @end deftypefn
%@eod:


% Copyright (C) 2012, 2013 Dynare Team
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

[basename, ext] = strtok(file,'.');
    
run(basename);

if exist('INIT__','var')
    init = dynDate(INIT__);
    clear('INIT__')
else
    init = dynDate(1);
end

if exist('FREQ__','var')
    freq = FREQ__;
    clear('FREQ__');
else
    freq = 1;
end

if exist('NAMES__','var')
    varlist0 = NAMES__;
    clear('NAMES__');
else
    varlist0 = [];
    list_of_variables = [];
end

if exist('TEX__','var')
    tex = TEX__;
    clear('TEX__');
else
    tex = [];
end


if isempty(varlist0)
    list_of_variables = whos();
end

data = [];
varlist = {};

if isempty(varlist0)
    for i=1:length(list_of_variables)
        if isequal(list_of_variables(i).name,'freq') || isequal(list_of_variables(i).name,'time') || isequal(list_of_variables(i).name,'data') ...
                || isequal(list_of_variables(i).name,'varlist') ...
                || isequal(list_of_variables(i).name,'varlist0') ...
                || isequal(list_of_variables(i).name,'list_of_variables') ...
                || isequal(list_of_variables(i).name,'tex') ...                
            continue
        end
        if list_of_variables(i).global || list_of_variables(i).persistent
            continue
        end
        if list_of_variables(i).complex || ~strcmp(list_of_variables(i).class,'double')
            continue
        end
        try
            eval(['data = [data, ' list_of_variables(i).name '];'])
            eval(['varlist = {varlist{:}, ''' list_of_variables(i).name '''};']) 
        catch
            error(['load_m_file:: All the vectors (variables) in ' inputname(1) ' must have the same number of rows (observations)!'])
        end
    end
else
    for i=1:length(varlist0)
       eval(['data = [data, ' varlist0{i} '];']) 
    end
    varlist = varlist0;
end

%@test:1
%$ % Create a data m-file
%$ fid = fopen('data_m_file.m','w');
%$ fprintf(fid,'FREQ__ = 4;');
%$ fprintf(fid,'INIT__ = ''1938Q4'';');
%$ fprintf(fid,'NAMES__ = {''azert'';''yuiop''};');
%$ fprintf(fid,'TEX__ = {''azert'';''yuiop''};');
%$ fprintf(fid,'azert = [1; 2; 3; 4; 5];');
%$ fprintf(fid,'yuiop = [2; 3; 4; 5; 6];');
%$ fclose(fid);
%$
%$ % Try to read the data m-file
%$ try
%$     datafile = 'data_m_file';
%$     [freq,init,data,varlist,tex] = load_m_file_data(datafile);
%$     t(1) = 1;
%$ catch exception
%$     t(1) = 0;
%$     T = all(t);
%$     LOG = getReport(exception,'extended');
%$     return
%$ end
%$
%$ % Check the results.
%$ t(2) = dyn_assert(freq,4);
%$ t(3) = dyn_assert(isa(init,'dynDate'),1);
%$ t(4) = dyn_assert(init.freq,4);
%$ t(5) = dyn_assert(init.time,[1938 4]);
%$ t(6) = dyn_assert(varlist,{'azert';'yuiop'});
%$ t(7) = dyn_assert(tex,{'azert';'yuiop'});
%$ t(8) = dyn_assert(data(:,1),[1;2;3;4;5]);
%$ t(9) = dyn_assert(data(:,2),[2;3;4;5;6]);
%$ T = all(t);
%@eof:1