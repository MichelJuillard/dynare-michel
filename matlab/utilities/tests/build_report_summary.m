function str = build_report_summary(reportfile, printonscreen, mailreport)

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
    mailreport = 0;
end

if nargin<3
    mailreport = 0;
else
    if ischar(mailreport)
        mailto = mailreport;
        mailreport = 1;
        fid = fopen('~/.matlab-send-report');
        system(['scp ' reportfile ' ' fgetl(fid)]);
        fclose(fid);
    else
        if ~isequal(mailreport,0)
            error('build_report_summary:: Third argument must be an adress email!')
        end
    end
end

reportfilecontent = load(reportfile);
reportcell = reportfilecontent.report;
gitinfo = reportfilecontent.gitinfo;
gitlastcommithash = reportfilecontent.gitlastcommithash;

str = 'Hi,';
str = char(str,'');
str = char(str,'This is a summary report for the unitary tests in Dynare. Full report can be');
str = char(str,['found at http://www.dynare.org/stepan/dynare/tests/' reportfile]);
str = char(str,'');

str = char(str,gitinfo(1,:));
str = char(str,gitinfo(2,:));

str = char(str,'');
str = char(str,'');

str = char(str,['===========================']);
str = char(str,'DYNARE/MATLAB UNITARY TESTS');
str = char(str,'===========================');
str = char(str,['| TOTAL: ' int2str(size(reportcell,1))]);
str = char(str,['|  PASS: ' int2str(length(find(cell2mat(reportcell(:,3)))))]);
str = char(str,['|  FAIL: ' int2str(length(find(~cell2mat(reportcell(:,3)))))]);
str = char(str,'|');
str = char(str,'| LIST OF FAILED TESTS:');
for i=1:size(reportcell,1)
    if ~reportcell{i,3}
        str = char(str,['|    * ' reportcell{i,1} ' (test #'  int2str(reportcell{i,2}) '{' strrep(int2str(transpose(find(~reportcell{i,4}))),' ',',') '})']);
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