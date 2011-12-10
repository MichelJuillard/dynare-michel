function bytecode_debug()
global M_ oo_ options_ AA;
fid = fopen([M_.fname '_options.txt'],'wt');
nfields = fieldnames(options_);
fprintf(fid, '%d %d %d\n',size(nfields,1), size(options_,1), size(options_,2));
for i=1:size(nfields, 1)
  disp(nfields(i));
  if iscell(nfields(i))
    AA = cell2mat(nfields(i));
  else
    AA = nfields(i);
  end;
  if iscell(AA)
    AA = cell2mat(AA);
  end;
  fprintf(fid, '%s\n', AA);
  Z = getfield(options_, AA);
  print_object(fid, Z);
end;
fclose(fid);

fid = fopen([M_.fname '_M.txt'],'wt');
nfields = fields(M_);
fprintf(fid, '%d %d %d\n',size(nfields,1), size(M_,1), size(M_,2));
for i=1:size(nfields, 1)
  disp(nfields(i));
  if iscell(nfields(i))
    AA = cell2mat(nfields(i));
  else
    AA = nfields(i);
  end;
  fprintf(fid, '%s\n', AA);
  print_object(fid, getfield(M_, AA));
end;
fclose(fid);


fid = fopen([M_.fname '_oo.txt'],'wt');
nfields = fields(oo_);
fprintf(fid, '%d %d %d\n',size(nfields,1), size(oo_,1), size(oo_,2));
for i=1:size(nfields, 1)
  disp(nfields(i));
  if iscell(nfields(i))
    AA = cell2mat(nfields(i));
  else
    AA = nfields(i);
  end;
  if iscell(AA)
    AA = cell2mat(AA);
  end;
  fprintf(fid, '%s\n', AA);
  print_object(fid, getfield(oo_, AA));
end;
fclose(fid);

function print_object(fid, object_arg)
 if iscell(object_arg)
   object = cell2mat(object_arg);
 else
   object = object_arg;
 end;
 if isa(object,'float') == 1
   fprintf(fid, '%d ', 0);
   fprintf(fid, '%d %d\n',size(object,1), size(object,2));
   fprintf(fid, '%f\n', object);
   %for i=1:size(object, 2) 
     %for j=1:size(object, 1)
       %fprintf(fid, '%f\n', object(i,j));
     %end;
   %end;
 elseif isa(object,'char') == 1
   fprintf(fid, '%d ', 3);
   fprintf(fid, '%d %d\n',size(object,1), size(object,2));
   %object
   for i=1:size(object, 1)
     %for i=1:size(object, 2)
       fprintf(fid, '%s ', object(i,:));
     %end;
     %fprintf(fid, '\n');
   end;
   fprintf(fid, '\n');
 elseif isa(object,'struct') == 1
   fprintf(fid, '%d ', 5);
   nfields = fields(object);
   fprintf(fid, '%d %d %d\n',size(nfields,1), size(object,1), size(object,2));
   for j=1:size(object, 1) * size(object, 2)
     nfields = fields(object(j));
     for i=1:size(nfields, 1)
       if iscell(nfields(i))
         AA = cell2mat(nfields(i));
       else
         AA = nfields(i);
       end;
       fprintf(fid, '%s\n', AA);
       print_object(fid, getfield(object, AA));
     end;
   end;
 else
   disp(['type ' object  'note handle']);
 end;




