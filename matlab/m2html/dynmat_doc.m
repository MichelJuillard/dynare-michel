evalin('base',['cd ' strrep(which('dynare.m'),'dynare.m','') '..'])
m2html('mfiles', 'matlab', 'htmldir','doc/matlab','graph','on');