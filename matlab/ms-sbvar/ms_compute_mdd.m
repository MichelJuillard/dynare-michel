function oo_=ms_compute_mdd(options_,oo_)
%function ms_compute_mdd()
% calls ms sbvar mdd mex file
%
% INPUTS
%    options_
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2011 Dynare Team
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

compute_mdd = create_mdd_commandline(options_);
[err] = ms_sbvar_command_line(compute_mdd);
mexErrCheck('ms_sbvar_command_line mdd',err);

% grab the muller/bridge log mdd from the output file
mull_exp = 'Muller \w+\(\w+\) \= (\d+.\d+e\+\d+)';
bridge_exp = 'Bridge \w+\(\w+\) \= (\d+.\d+e\+\d+)';
bridge_mdd = -1; muller_mdd = -1;
mdd_filename = ['mdd_t' num2str(options_.ms.mdd_proposal_type(1)) '_' options_.ms.output_file_tag '.out'];
if exist(mdd_filename,'file')
    mdd_fid = fopen(mdd_filename);
    tline = fgetl(mdd_fid);
    while ischar(tline)
        mull_tok = regexp(tline,mull_exp,'tokens');
        bridge_tok = regexp(tline,bridge_exp,'tokens');
        if (~isempty(mull_tok))
            muller_mdd = str2double(mull_tok{1}{1});
        end
        if (~isempty(bridge_tok))
            bridge_mdd = str2double(bridge_tok{1}{1});
        end
        tline = fgetl(mdd_fid);
    end
    oo_.ms.mueller_log_mdd = muller_mdd;
    oo_.ms.bridged_log_mdd = bridge_mdd;
end
end

function opt=create_mdd_commandline(options_)

opt = '-mdd';
opt = [opt set_flag_if_exists(options_.ms, 'output_file_tag')];
opt = [opt set_flag_if_exists(options_.ms, 'load_mh_file ')];
opt = [opt set_flag_if_exists(options_.ms, 'mdd_proposal_draws ')];
opt = [opt set_flag_if_exists(options_.ms, 'mdd_use_mean_center')];

if size(options_.ms.mdd_proposal_type, 2) ~= 3
    error('ERROR: mdd_proposal_type takes exactly 3 arguments');
else
    opt = [opt ' -pt ' num2str(options_.ms.mdd_proposal_type(1))];
    opt = [opt ' -l '  num2str(options_.ms.mdd_proposal_type(2))];
    opt = [opt ' -h '  num2str(options_.ms.mdd_proposal_type(3))];
end

if options_.DynareRandomStreams.seed
    opt = [' -seed' num2str(options_.DynareRandomStreams.seed)];
end
end
