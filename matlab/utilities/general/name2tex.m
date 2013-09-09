function tex = name2tex(name, info) % --*-- Unitary tests --*--

% Converts plain text name into tex name.
 
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

if nargin<2
    info = 0;
end

if iscell(name)
    nn = length(name);
else
    nn = 1;
end

if iscell(name)
    tex = cell(nn,1);
end

for i=1:nn
    if info
        if iscell(name)
            tex(i) = {name{i}};
            id = findstr(name{i},'_');
            len = length(tex{i});
        else
            tex = name;
            id = findstr(name,'_');
            len = length(tex);
        end
        if isempty(id)
            continue
        end
        n = length(id);
        if id(1)==1
            if iscell(name)
                tex(i) = {['\_', tex{i}(2:end)]};
            else
                tex = ['\_', tex(2:end)];
            end
            if n>1
                id = id(2:end)+1;
                n = n-1;
                len = len+1;
            else
                continue
            end
        end
        if id(end)==len
            if iscell(name)
                tex(i) = {[tex{i}(1:end-1) '\_']};
            else
                tex = [tex(1:end-1) '\_'];
            end
            if n>1
                id = id(1:end-1);
                n = n-1;
            else
                continue
            end
        end
        if n==1
            if iscell(name)
                tex(i) = {[ tex{i}(1:(id-1)) '_{'  tex{i}((id+1):end) '}' ]};
            else
                tex = [ tex(1:(id-1)) '_{'  tex((id+1):end) '}' ];
            end
            continue
        else
            for j=1:n-1
                if iscell(name)
                    tex(i) = {[tex{i}(1:id(j)-1) '\_' tex{i}((id(j)+1):end)]};
                else
                    tex = [tex(1:id(j)-1) '\_' tex((id(j)+1):end)];
                end
                id = id(2:end)+1;
            end
            if iscell(name)
                tex(i) = {[tex{i}(1:(id-1)) '_{'  tex{i}((id+1):end) '}']};
            else
                tex = [tex(1:(id-1)) '_{'  tex((id+1):end) '}'];
            end
        end
    else
        if iscell(name)
            tex(i) = {strrep(name{i}, '_', '\_')};
        else
            tex = strrep(name, '_', '\_');
        end
    end
end

%@test:1
%$ t = zeros(16,1);
%$ t1 = name2tex('_azert');
%$ t2 = name2tex('azert_');
%$ t3 = name2tex('_azert_');
%$ t4 = name2tex('azert_uiop');
%$ t5 = name2tex('azert_uiop_qsdfg');
%$ t6 = name2tex('azert_uiop_qsdfg_');
%$ t7 = name2tex('_azert_uiop_qsdfg');
%$ t8 = name2tex('_azert_uiop_qsdfg_');
%$ t11 = name2tex('_azert',1);
%$ t12 = name2tex('azert_',1);
%$ t13 = name2tex('_azert_',1);
%$ t14 = name2tex('azert_uiop',1);
%$ t15 = name2tex('azert_uiop_qsdfg',1);
%$ t16 = name2tex('azert_uiop_qsdfg_',1);
%$ t17 = name2tex('_azert_uiop_qsdfg',1); 
%$ t18 = name2tex('_azert_uiop_qsdfg_',1);
%$ 
%$ t(1) = dyn_assert(strcmp(t1,'\\_azert'),1);
%$ t(2) = dyn_assert(strcmp(t2,'azert\\_'),1);
%$ t(3) = dyn_assert(strcmp(t3,'\\_azert\\_'),1);
%$ t(4) = dyn_assert(strcmp(t4,'azert\\_uiop'),1);
%$ t(5) = dyn_assert(strcmp(t5,'azert\\_uiop\\_qsdfg'),1);
%$ t(6) = dyn_assert(strcmp(t6,'azert\\_uiop\\_qsdfg\\_'),1);
%$ t(7) = dyn_assert(strcmp(t7,'\\_azert\\_uiop\\_qsdfg'),1);
%$ t(8) = dyn_assert(strcmp(t8,'\\_azert\\_uiop\\_qsdfg\\_'),1);
%$ t(9) = dyn_assert(strcmp(t11,'\\_azert'),1);
%$ t(10) = dyn_assert(strcmp(t12,'azert\\_'),1);
%$ t(11) = dyn_assert(strcmp(t13,'\\_azert\\_'),1); 
%$ t(12) = dyn_assert(strcmp(t14,'azert_{uiop}'),1);
%$ t(13) = dyn_assert(strcmp(t15,'azert\\_uiop_{qsdfg}'),1);
%$ t(14) = dyn_assert(strcmp(t16,'azert\\_uiop_{qsdfg\\_}'),1);
%$ t(15) = dyn_assert(strcmp(t17,'\\_azert\\_uiop_{qsdfg}'),1);
%$ t(16) = dyn_assert(strcmp(t18,'\\_azert\\_uiop_{qsdfg\\_}'),1);
%$
%$ T = all(t);
%@eof:1

%@test:2
%$ t = zeros(16,1);
%$ t1 = name2tex({'_azert'});
%$ t2 = name2tex({'azert_'});
%$ t3 = name2tex({'_azert_'});
%$ t4 = name2tex({'azert_uiop'});
%$ t5 = name2tex({'azert_uiop_qsdfg'});
%$ t6 = name2tex({'azert_uiop_qsdfg_'});
%$ t7 = name2tex({'_azert_uiop_qsdfg'});
%$ t8 = name2tex({'_azert_uiop_qsdfg_'});
%$ t11 = name2tex({'_azert'},1);
%$ t12 = name2tex({'azert_'},1);
%$ t13 = name2tex({'_azert_'},1);
%$ t14 = name2tex({'azert_uiop'},1);
%$ t15 = name2tex({'azert_uiop_qsdfg'},1);
%$ t16 = name2tex({'azert_uiop_qsdfg_'},1);
%$ t17 = name2tex({'_azert_uiop_qsdfg'},1);
%$ t18 = name2tex({'_azert_uiop_qsdfg_'},1);
%$
%$ t(1) = dyn_assert(t1,{'\\_azert'});
%$ t(2) = dyn_assert(t2,{'azert\\_'});
%$ t(3) = dyn_assert(t3,{'\\_azert\\_'});
%$ t(4) = dyn_assert(t4,{'azert\\_uiop'});
%$ t(5) = dyn_assert(t5,{'azert\\_uiop\\_qsdfg'});
%$ t(6) = dyn_assert(t6,{'azert\\_uiop\\_qsdfg\\_'});
%$ t(7) = dyn_assert(t7,{'\\_azert\\_uiop\\_qsdfg'});
%$ t(8) = dyn_assert(t8,{'\\_azert\\_uiop\\_qsdfg\\_'});
%$ t(9) = dyn_assert(t11,{'\\_azert'});
%$ t(10) = dyn_assert(t12,{'azert\\_'});
%$ t(11) = dyn_assert(t13,{'\\_azert\\_'});
%$ t(12) = dyn_assert(t14,{'azert_{uiop}'});
%$ t(13) = dyn_assert(t15,{'azert\\_uiop_{qsdfg}'});
%$ t(14) = dyn_assert(t16,{'azert\\_uiop_{qsdfg\\_}'});
%$ t(15) = dyn_assert(t17,{'\\_azert\\_uiop_{qsdfg}'});
%$ t(16) = dyn_assert(t18,{'\\_azert\\_uiop_{qsdfg\\_}'});
%$
%$ T = all(t);
%@eof:2

%@test:3
%$ t = zeros(4,1);
%$ try
%$     t1 = name2tex({'_azert';'azert_';'_azert_';'azert_uiop';'azert_uiop_qsdfg';'azert_uiop_qsdfg_'});
%$     t(1) = 1;
%$ catch
%$     % Nothing to do here.
%$ end
%$
%$ if t(1)
%$     try
%$         t2 = name2tex({'_azert';'azert_';'_azert_';'azert_uiop';'azert_uiop_qsdfg';'azert_uiop_qsdfg_'},1);
%$         t(2) = 1;
%$     catch
%$         % Nothing to do here.
%$     end
%$ end
%$
%$ if t(1)
%$     t(3) = dyn_assert(t1,{'\\_azert';'azert\\_';'\\_azert\\_';'azert\\_uiop';'azert\\_uiop\\_qsdfg';'azert\\_uiop\\_qsdfg\\_'});
%$ end
%$
%$ if t(2)
%$     t(4) = dyn_assert(t2,{'\\_azert';'azert\\_';'\\_azert\\_';'azert_{uiop}';'azert\\_uiop_{qsdfg}';'azert\\_uiop_{qsdfg\\_}'});
%$ end
%$
%$ T = all(t);
%@eof:3