function dyn_data_01=read_variables(file_name_01,var_names_01,dyn_data_01,xls_sheet,xls_range)
  
% function dyn_data_01=read_variables(file_name_01,var_names_01,dyn_data_01,xls_sheet,xls_range)
% Read data
%
% INPUTS
%    file_name_01:    file name
%    var_names_01:    variables name
%    dyn_data_01:     
%    xls_sheet:       Excel sheet name
%    xls_range:       Excel range specification
%
% OUTPUTS
%    dyn_data_01:
%
% SPECIAL REQUIREMENTS
% all local variables have complicated names in order to avoid name
% conflicts with possible user variable names
%  
% part of DYNARE, copyright Dynare Team (2005-2007)
% Gnu Public License.


  dyn_size_01 = size(dyn_data_01,1);
  var_size_01 = size(var_names_01,1);
  if exist(file_name_01)
    file_name_02 = [file_name_01 '.m'];
    dyn_instr_01 = file_name_01;
    eval(dyn_instr_01);
    for dyn_i_01=1:var_size_01
      dyn_tmp_01 = eval(var_names_01(dyn_i_01,:));
      if length(dyn_tmp_01) > dyn_size_01 & dyn_size_01 > 0
	error('data size is too large')
      end
      dyn_data_01(:,dyn_i_01) = dyn_tmp_01;
    end
  elseif exist([file_name_01 '.mat'])
    file_name_02 = [file_name_01 '.mat']
    s = load(file_name_01);
    for dyn_i_01=1:var_size_01
      dyn_tmp_01 = s.(deblank(var_names_01(dyn_i_01,:)));
      if length(dyn_tmp_01) > dyn_size_01 & dyn_size_01 > 0
	error('data size is too large')
      end
      dyn_data_01(:,dyn_i_01) = dyn_tmp_01;
    end
  elseif exist([file_name_01 '.xls'])
    file_name_02 = [file_name_01 '.xls'];
    [num,txt,raw] = xlsread(file_name_01,xls_sheet,xls_range);
    for dyn_i_01=1:var_size_01
      iv = strmatch(var_names_01(dyn_i_01,:),raw(1,:),'exact');
      dyn_tmp_01 = [raw{2:end,iv}]';
      if length(dyn_tmp_01) > dyn_size_01 & dyn_size_01 > 0
	error('data size is too large')
      end
      dyn_data_01(:,dyn_i_01) = dyn_tmp_01;
    end
    
  end
  disp(sprintf('Loading %d observations from %s\n',...
	       size(dyn_data_01,1),file_name_02))
