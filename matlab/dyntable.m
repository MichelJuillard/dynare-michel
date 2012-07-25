function dyntable(title,headers,labels,values,label_width,val_width, ...
                  val_precis)

% Copyright (C) 2002-2011 Dynare Team
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

global options_

if options_.noprint
    return
end

%label_width = max(size(deblank(char(headers(1,:),labels)),2)+2, ...
%                  label_width);
label_width = max(size(deblank(char(headers(1,:),labels)),2))+2;

label_fmt = sprintf('%%-%ds',label_width);

values_length = max(ceil(max(max(log10(abs(values(isfinite(values))))))),1)+val_precis+1;
if any(values) < 0
    values_length = values_length+1;
end
headers_length = max(size(deblank(headers(2:end,:)),2));
%val_width = max(values_length,val_width);
val_width = max(headers_length,values_length)+2;
%header_fmt = sprintf('%%-%ds',val_width);
val_fmt = sprintf('%%%d.%df',val_width,val_precis);

headers_offset = 0;
values_offset = 0;
if headers_length > values_length
    %    values_offset = ceil((val_width-values_length)/2);
end

if length(title) > 0
    disp(sprintf('\n\n%s\n',title));
end
if length(headers) > 0
    hh = sprintf(label_fmt,headers(1,:));
    for i=2:size(headers,1)
        hla = size(deblank(headers(i,:)),2);
        hlb = ceil((val_width - hla)/2);
        hla = val_width - hla - hlb;
        hh  = [hh char(32*ones(1,hlb)) deblank(headers(i,:)) ...
              char(32*ones(1,hla))];
    end
    disp(hh);
end
for i=1:size(values,1)
    disp([sprintf(label_fmt,deblank(labels(i,:))) char(32*ones(1,values_offset)) sprintf(val_fmt,values(i,:))]);
end

% 10/30/02 MJ