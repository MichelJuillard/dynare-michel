function set_param_value(pname,value)
  global M_
  
  i = strmatch(pname,M_.param_names,'exact');
  
  if isempty(i)
    error(['Parameter name ' pname ' doesn''t exist'])
  end
  
  M_.params(i) = value;
  assignin('base',pname,value);