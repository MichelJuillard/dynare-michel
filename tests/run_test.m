function run_test()
test_files = {
'.'  'ramst';
'.'  'ramst_a';
'.' 'example1';
'.' 'example2';
'.' 't_sgu_ex1';
'arima' 'mod1';
'arima' 'mod1a';
'arima' 'mod1b';
'arima' 'mod1c';
'arima' 'mod2';
'arima' 'mod2a';
'arima' 'mod2b';
'arima' 'mod2c';
'fs2000' 'fs2000';
'fs2000' 'fs2000a';
}

results = cell(length(test_files),1);

for i=1:length(test_files)
     results{i}= run_test1(test_files{i,1},test_files{i,2});
end

for i=1:length(test_files)
  disp(test_files{i,2})
  disp(results{i})
end
function msg=run_test1(path1,mod_file)
     global options_
     clear options_
     old_path = pwd;
     cd(path1);
     msg = 'OK';
     expr = ['disp(''error in ' mod_file ''');msg=lasterr;disp(msg)'];
     eval(['dynare ' mod_file ' noclearall'],'eval(expr)');
     cd(old_path)
     
