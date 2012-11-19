function [list_of_variables, data, time] = readcsv(file, withtime, withnames, noemptycell)

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

if ~withtime && ~withname && noemptycell
    % Use matlab builtin routine!
    data = csvread(file);
end

if ~( isequal(withtime,0) || isequal(withtime,1) )
    error('readcsv:: Second input argument has to be equal to 1 or 0!')
end

if ~( isequal(withnames,0) || isequal(withnames,1) )
    error('readcsv:: Third input argument has to be equal to 1 or 0!')
end

% Output initialization 
time = []; data = []; list_of_variables = [];

% Check if file exists.
if check_file_extension(file,'csv')
    try
        fid = fopen(file,'r');
    catch
        error(['readcsv: I can''t find file ' file '!'])
    end
else
    error('readcsv: Wrong file extension!')
end

% bfile contains a vector of ascii codes.
bfile = fread(fid);

% Close (csv) file. 
fclose(fid);

% Set newline code (ok for *nix, check for mac and windows)
if isunix
    newline_code = 10;
else
    error('readcsv:: Not implemented for your OS!')
end

% Get the positions of the end-of-line code;
end_of_line_locations = find(bfile==newline_code);
tmp = find(bfile==newline_code);

% Get the number of lines in the file.
ndx = length(tmp);                       

% Create a cell of indices for each line.
b = [1; end_of_line_locations+1];
c = [end_of_line_locations-1; length(bfile)+1];
b = b(1:end-1);
c = c(1:end-1);

linea = 1;

if withnames
    % Get the first line of the csv file (names of the variables).
    linee = char(transpose(bfile(b(linea):c(linea))));
    % Get the content of the first line and determine the number of variables and their names.
    [B,C] = get_cells_id(linee,',');
    if withtime
        B = B(2:end);
        C = C(2:end);
    end
    list_of_variables = cell(length(B),1);
    number_of_variables = length(list_of_variables);
    for i=1:number_of_variables
        list_of_variables(i) = {linee(B(i):C(i))};
    end
    linea = linea+1;
end

% Get following line (number 1 or 2 depending on withnames flag)
linee = char(transpose(bfile(b(linea):c(linea))));
comma_locations = transpose(strfind(linee,','));
B = 1;
C = comma_locations(1)-1;

if withtime
    tmp = linee(B:C);
    % Check the dates formatting
    if ~(isyearly(tmp) || isquaterly(tmp) || ismonthly(tmp) || isweekly(tmp))
        error('readcsv:: Formatting error. I can''t read the dates!')
    end
    if isyearly(tmp)==2
        % Remove the Y (gpm/iris date format) if necessary
        tmp = { tmp(1:end-1) };
    end
    initial_date = dynDate(tmp);
    first = 2;
else
    initial_date = dynDate(1);
    first = 1;
end

if ~withnames
    number_of_variables = length(tmp)-withtime;
end

% Initialization of matrix data.
data = zeros(ndx,number_of_variables);

% Populate data.
for linea = 1+withnames:ndx
    linee = char(transpose(bfile(b(linea):c(linea))));
    [B,C] = get_cells_id(linee,',');
    for i=first:length(B)
        if isequal(B(i),C(i))
            data(linea,i-withtime) = NaN;
        else
            data(linea,i-withtime) = str2double(linee(B(i):C(i)));
        end
    end
end