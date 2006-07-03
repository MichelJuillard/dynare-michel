function []=disp_dr(dr,order,var_list)
// Copyright (C) 2001 Michel Juillard
// 
nx = size(dr('ghx'),2);
nu = size(dr('ghu'),2);
 
k1 = dr('order_var');
kstate = dr('kstate');
k = find(kstate(:,2) <= ykmin_+1);
klag = kstate(k,[1 2]);
  
  
nvar = size(var_list,1);
if nvar == 0 then
  nvar = endo_nbr;
  ivar = [1:nvar];
else
  ivar=zeros(nvar,1);
  for i=1:nvar
    i_tmp = grep_exact(lgy_(k1,:),var_list(i,:));
    if isempty(i_tmp) then
      dyn_disp(var_list(i,:));
      error (['One of the variable specified does not exist']) ;
    else
      ivar(i) = i_tmp;
    end
  end
end

dyn_disp('POLICY AND TRANSITION FUNCTIONS');
// variable names
str = '                    ';
for i = 1:nvar
  str = str+msprintf('%16s',lgy_(k1(ivar(i))));
end
mprintf('%s\n',str);
mfprintf(fh_log,'%s\n',str);
// 
// constant
// 
str = 'Constant            ';
flag = 0;
for i = 1:nvar
  x = dr.ys(k1(ivar(i)));
  if order == 2 then
    x= x+dr.ghs2(ivar(i))/2;
  end
  if abs(x) > .000001 then
    flag = 1;
    str = str+msprintf('%16.6f',x);
  else
    str = str+'               0';
  end
end
if flag then
  mprintf('%s\n',str);
  mfprintf(fh_log,'%s\n',str);
end
if order == 2 then
  str = '(Correction)        ';
  flag = 0;
  for i = 1:nvar
    x = dr.ghs2(ivar(i))/2;
    if abs(x) > 1e-6 then
      flag = 1;
      str = str+msprintf('%16.6f',x)
    else
      str = str+'               0';
    end
  end
  if flag then
    mprintf('%s\n',str);
    mfprintf(fh_log,'%s\n',str);
  end
end
// 
// ghx
// 
for k = 1:nx
  flag = 0;
  str1 = msprintf('%s(%d)',lgy_(k1(klag(k,1))),klag(k,2)-ykmin_-2);
  str = msprintf('%-20s',str1);
  for i = 1:nvar
    x = dr.ghx(ivar(i),k);
    if abs(x)>.000001 then
      flag = 1;
      str = str+msprintf('%16.6f',x);
    else
      str = str+'               0';
    end
  end
  if flag then
    mprintf('%s\n',str);
    mfprintf(fh_log,'%s\n',str);
  end
end
// 
// ghu
// 
for k = 1:nu
  flag = 0;
  str = msprintf('%-20s',lgx_(k));
  for i = 1:nvar
    x = dr('ghu')(ivar(i),k);
    if abs(x)>.000001 then
      flag = 1;
      str = str+msprintf('%16.6f',x);
    else
      str = str+'               0';
    end
  end
  if flag then
    mprintf('%s\n',str);
    mfprintf(fh_log,'%s\n',str);
  end
end
 
if order>1 then
  // ghxx
  for k = 1:nx
    for j = 1:k
      flag = 0;
      str = msprintf('%-20s',lgy_(k1(dr('nstatic')+k),:)+'(-1)'+lgy_(k1(dr('nstatic')+j),:)+'(-1)');
      for i = 1:endo_nbr
        if k==j then
          x = dr.ghxx(i,(k-1)*nx+j)/2;
        else
          x = dr.ghxx(i,(k-1)*nx+j);
        end
        if abs(x)>.000001 then
          flag = 1;
          str = str+msprintf('%16.6f',x);
        else
          str = str+'               0';
        end
      end
      if flag then
        mprintf('%s\n',str);
        mfprintf(fh_log,'%s\n',str);
      end
    end
  end
  // 
  // ghuu
  // 
  for k = 1:nu
    for j = 1:k
      flag = 0;
      str = msprintf('%-20s',lgx_(k,:)+' '+lgx_(j,:));
      for i = 1:endo_nbr
        if k==j then
          x = dr.ghuu(i,(k-1)*nu+j)/2;
        else
          x = dr.ghuu(i,(k-1)*nu+j);
        end
        if abs(x)>.000001 then
          flag = 1;
          str = str+msprintf('%16.6f',x);
        else
          str = str+'               0';
        end
      end
      if flag then
        mprintf('%s\n',str);
        mfprintf(fh_log,'%s\n',str);
      end
    end
  end
  // 
  // ghxu
  // 
  for k = 1:nx
    for j = 1:nu
      flag = 0;
      str = msprintf('%-20s',lgy_(k1(dr('nstatic')+k),:)+'(-1) '+lgx_(j,:));
      for i = 1:endo_nbr
        x = dr.ghxu(i,(k-1)*nu+j);
        if abs(x)>.000001 then
          flag = 1;
          str = str+msprintf('%16.6f',x);
        else
          str = str+'               0';
        end
      end
      if flag then
        mprintf('%s\n',str);
        mfprintf(fh_log,'%s\n',str);
      end
    end
  end
end
 
// 01/08/2001 MJ added test for order in printing quadratic terms
// 02/21/2001 MJ pass all variable names through deblank()
// 02/21/2001 MJ changed from f to g format to write numbers
// 02/21/2003 MJ improved display format

 
 
 
