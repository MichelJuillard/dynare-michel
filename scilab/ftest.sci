function []=ftest(s1,s2)
// Copyright (C) 2001 Michel Juillard
// 
 
global('nvx','nvy','x','y','lag1');
 
if size(s1,1)~=2 then
  error('Spécifiez deux fichiers pour la comparaison.');
end
 
for i = 1:2
  if ~(matrix(find(abs(s1(i,:))==46),1,-1)==[]) then
    error('Entrez les noms de fichiers sans extensions.');
  end
end
 
s1 = s1+[' ';' '];
file1 = part(s1(1),1:min(find(abs(part(s1(1),:))==32))-1)+'.BIN';
file2 = part(s1(2),1:min(find(abs(part(s1(2),:))==32))-1)+'.BIN';
 
[fid,%v] = mopen(file1,'r',0)
if %v<0 then fid = -1;end
n1 = mtlb_fread(fid,1,'int');
n2 = mtlb_fread(fid,1,'int');
n3 = mtlb_fread(fid,1,'int');
lag1 = mtlb_fread(fid,4,'int');
nvx = mtlb_fread(fid,[n1,n3],'int');
x = mtlb_fread(fid,[n1,n2],'float64');
mclose(fid);
nvx = ascii(nvx);
 
[fid,%v] = mopen(file2,'r',0)
if %v<0 then fid = -1;end
n1 = mtlb_fread(fid,1,'int');
n2 = mtlb_fread(fid,1,'int');
n3 = mtlb_fread(fid,1,'int');
lag2 = mtlb_fread(fid,4,'int');
nvy = mtlb_fread(fid,[n1,n3],'int');
y = mtlb_fread(fid,[n1,n2],'float64');
mclose(fid);
nvy = ascii(nvy);
 
if size(x,1)~=size(y,1) then
  error('FTEST: The two files don''t have the same number of variables.');
end
 
for i = 1:size(x,1)
  if ~(part(nvx(i),:)==part(nvy(i),:)) then
    error('FTEST: The two files don''t have the same  variables.');
  end
end
 
if nnz(lag1-lag2)>0 then
  error('FTEST: Leads and lags aren''t the same in both files.');
end
 
j = zeros(size(s2,1),1);
for i = 1:size(s2,1)
   
  //!! Unknown function strmatch ,the original calling sequence is used
  k = strmatch(s2(i,:),nvx,'exact');
  if k==[] then
    t = 'FTEST: Variable '+s2(i)+'doesn''t exist';
    error(t);
  else
    j(i,1) = k(:);
  end
end
 
y = y(j,:);
x = x(j,:);
 
//06/18/01 MJ replaced beastr by strmatch
