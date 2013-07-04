function str = build_report_summary(report, printonscreen, mailreport)

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

if isequal(nargin,1)
    printonscreen = 1;
end

if nargin<3
    mailreport = 0;
else
    if ischar(mailreport)
        mailto = mailreport;
        mailreport = 1;
    else
        error('build_report_summary:: Third argument must be an adress email!')
    end
        
end

system('git show --pretty=format:"Last commit %H by %an, %ar %n-> %s" HEAD > git.info');

str = 'Hi,';
str = char(str,'');
str = char(str,'This is a summary report for the unitary tests in Dynare.');
str = char(str,'');

fid = fopen('git.info');
str = char(str,fgetl(fid));
str = char(str,fgetl(fid));
fclose(fid);

str = char(str,'');
str = char(str,'');

str = char(str,['===========================']);
str = char(str,'DYNARE/MATLAB UNITARY TESTS');
str = char(str,'===========================');
str = char(str,['| TOTAL: ' int2str(size(report,1))]);
str = char(str,['|  PASS: ' int2str(length(find(cell2mat(report(:,3)))))]);
str = char(str,['|  FAIL: ' int2str(length(find(~cell2mat(report(:,3)))))]);
str = char(str,'|');
str = char(str,'| LIST OF FAILED TESTS:');
for i=1:size(report,1)
    if ~report{i,3}
        str = char(str,['|    * ' report{i,1} ' (test #'  int2str(report{i,2}) '{' strrep(int2str(transpose(find(~report{i,4}))),' ',',') '})']);
    end
end

if printonscreen
    fprintf('\n\n')
    disp(str)
    fprintf('\n\n')
end

if mailreport
    if exist('~/.matlab-send-mail-info','file')
        fid = fopen('~/.matlab-send-mail-info','r');
    else
       disp(['build_report_summary:: I Cannot send the report to ' mailto ' because the sender and the smtp server are not defined!'])
       disp(['                       You probably need to add a ''.matlab-send-mail-info'' in your root directory...'])
       return
    end
    setpref('Internet','SMTP_Server',fgetl(fid));
    setpref('Internet','SMTP_Username',fgetl(fid));
    setpref('Internet','E_mail',fgetl(fid));
    setpref('Internet','SMTP_Passeword',fgetl(fid));
    fclose(fid);
    STR = [deblank(str(1,:))];
    for i=2:size(str,1)
        STR = [STR 10 deblank(str(i,:)) ];
    end
    sendmail(mailto,'Dynare/Matlab unitary tests',STR); 
end