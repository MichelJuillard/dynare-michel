function mydelete(fname,pname)

if nargin ==0,
    disp('mydelete(fname)')
    return
end

if nargin ==1,
    pname='';
end

file_to_delete = dir([pname,fname]);
for j=1:length(file_to_delete),
    delete([pname,file_to_delete(j).name]);
end
