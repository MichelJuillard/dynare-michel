function run_test()
test_files = {
 '.'  'ramst';
 '.'  'ramst_a';
'.' 'example1';
%'.' 'example2';
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
  'estimation_options','fs2000A';
  'estimation_options','fs2000B';
             }

results = cell(length(test_files),1);

for i=1:length(test_files)
        tic;
       rand('state',1);
       randn('state',1);
     results{i}= run_test1(test_files{i,1},test_files{i,2});
    
    toc;
    
end


for i=1:length(test_files)
  disp(test_files{i,2})
  disp(results{i})
  disp(toc)
end

end

function msg=run_test1(path1,mod_file)
     global options_ oo_
     options_=struct('trick',1);
     oo_=struct('trick',1);
     old_path = pwd;
     cd(path1);
     msg = 'OK';
     expr = ['disp(''error in ' mod_file ''');msg=lasterr;disp(msg)'];
     
     try
         old_results=load([mod_file '_results']);
     catch
         eval(['dynare ' mod_file ' noclearall'],'eval(expr)');
         msg = 'no previous results'
         cd(old_path)
         return
     end
     eval(['dynare ' mod_file ' noclearall'],'eval(expr)');
     msg = strvcat(msg, comparison(oo_,old_results.oo_,0, 'oo_'));
     cd(old_path)
     end
     
function msg=comparison(A,B,tol,parent)

    msg = '';
    mca = my_class(A);
    mcb = my_class(B);
    if size(mca,2) ~= size(mcb,2) || any(mca ~= mcb)
        msg = ['Types differ in ' parent];
        return
    end

    switch my_class(A)
     case 'numeric'
      if size(A)==size(B)
          if ~all(all(abs(A-B)<tol))
              msg = ['Difference greater than tolerance in ' parent];
          end
      else
          msg = ['Sizes differ in ' parent];
      end
     case 'char'
      if A~=B
          msg = ['Strings differ in ' parent];
      end
     case 'cell' 
      if size(A) == size(B)
          [M,N] = size(A);
          for i = 1:M
              for j = 1:N
                  msg = strvcat(msg, comparison(A{i,j}, B{i,j}, tol, sprintf('%s{%d,%d}', parent, i, j)))
              end
          end
      else
          msg = [ 'Cell sizes differ in ' parent ];
      end
     case 'struct' 
      
      namesA=fieldnames(A);
      namesB=fieldnames(B);

      if length(namesA)~=length(namesB)
          msg = [ 'Number of fields differ in ' parent ]
          return
      end

      namesA = sort(namesA);
      namesB = sort(namesB);

      for i=1:length(namesA)
          if namesA{i}==namesB{i}
              msg = strvcat(msg, comparison(A.(namesA{i}),B.(namesB{i}),tol, sprintf('%s.%s', parent, namesA{i})));
          else
              msg = strvcat(msg, sprintf('Field %s.%s in 1st arg does not exist in 2nd arg', parent, namesA{i}));
          end    
      end
      
     otherwise
      msg = ['Unsupported type in ' parent ];
      
    end
end

function t = my_class(a)
    c = class(a);
    if isequal(c, 'double') | isequal(c, 'single') | isequal(c, 'int8')| isequal(c, 'uint8') | ...
            isequal(c, 'int16') | isequal(c, 'uint16') | isequal(c, 'int32') | isequal(c, 'uint32') | ...
            isequal(c, 'int64') | isequal(c, 'uint64')
        t = 'numeric';
    else
        t = c;
    end
end     
    
     
