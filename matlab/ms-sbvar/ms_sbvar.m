function ms_sbvar(ms_flag, M, options)

dynareroot = strrep(which('dynare'),'dynare.m','');
ms_root = [dynareroot '/ms-sbvar'];

addpath([ms_root '/cstz']);
addpath([ms_root '/identification']);
addpath([ms_root '/switching_specification']);
addpath([ms_root '/mhm_specification']);

clean_ms_files(M.fname);

options.data = read_variables(options.datafile,options.varobs,[],options.xls_sheet,options.xls_range);

if options.forecast == 0
    options.forecast = 4;
end

options.ms.output_file_tag = M.fname;
%options.ms.markov_file = 'specification_2v2c.dat';
%options.ms.mhm_file = 'MHM_input.dat';
%options.ms.restriction_fname = 'ftd_upperchol3v'; 

if options.ms.upper_cholesky
    if options.ms.lower_cholesky
        error(['Upper Cholesky and lower Cholesky decomposition can''t be ' ...
               'requested at the same time!'])
    else
        options.ms.restriction_fname = 'upper_cholesky';
    end
elseif options.ms.lower_cholesky
    options.ms.restriction_fname = 'lower_cholesky';
elseif ~isempty(options.ms.Qi) && ~isempty(options.ms.Ri)
    options.ms.restriction_fname = 'exclusions';
end

if ms_flag == 1
    % temporary fix
    options.ms.markov_file = 'markov_file';
    sz_prd(M,options);
else
    ms_mardd(options);
end
