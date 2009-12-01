function swz_sbvar(ms_flag, M, options)

dynareroot = strrep(which('dynare'),'dynare.m','');
swz_root = [dynareroot '/swz'];

addpath([swz_root '/cstz']);
addpath([swz_root '/identification']);
addpath([swz_root '/switching_specification']);
addpath([swz_root '/mhm_specification']);

options.data = read_variables(options.datafile,options.varobs,[],options.xls_sheet,options.xls_range);

options.ms.output_file_tag = M.fname;
%options.ms.markov_file = 'specification_2v2c.dat';
%options.ms.mhm_file = 'MHM_input.dat';
%options.ms.restriction_fname = 'ftd_upperchol3v'; 


if ms_flag == 1
    sz_prd(options)
else
    swz_mardd(options);
end
