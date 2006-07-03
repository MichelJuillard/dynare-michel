function []=simk()
// Copyright (C) 2001 Michel Juillard
// 
 
global y_ it_ jacobia_ d1_ it_ 

broyden_ = 0;

nk = ykmin_+ykmax_+1;
ny = size(iy_,2);
icc1 = iy_(nk,:)>0;
 
for i = 1:ykmax_-1
  icc1 = [iy_(nk-i,:)|icc1(1,:);icc1];
end
 
icc1 = find(icc1')';
iy = iy_>0;
isc = cumsum(sum(bool2s(iy'),1))';
iyr0 = matrix(find(iy_'),1,-1);
ncc1 = size(icc1,1);
ncc = ncc1+1;
ncs = size(iyr0,1);
 
ky = zeros(ny,nk);
// indices of variables at each lead or lag
lky = zeros(nk,1);
for i = 1:nk
  j = matrix(find(iy_(i,:)),1,-1)';
  if j==[] then
    lky(i) = 0
  else
    lky(i) = size(j,1)
    ky(1:lky(i),i) = j;
  end
end
 
jwc = matrix(find(iy(2:ykmax_+1,:)'),1,-1);
// indices of columns for
// triangularization
// as many rows as lags in model
 
if jwc==[] then
  jwc = 0;
  ljwc = 0;
  temp = icc1;
else
  ljwc = size(jwc,1);
  // length of each row in jwc
  temp = union(jwc,icc1)';
  // prepares next iteration
end
 
j1 = ky(1:lky(1),1);
lj1 = lky(1);
 
for i = 2:ykmin_
  [j1,lj1] = ffill(j1,lj1,selif(temp+(i-1)*ny,temp<=ny));
  if ykmax_==1 then
    if lky(i+ykmax_)>0 then
      [jwc,ljwc] = ffill(jwc,ljwc,ky(1:lky(i+ykmax_),i+ykmax_)+(ykmax_-1)*ny);
      if ljwc(i)==0 then
        temp = icc1;
      else
        temp = union(jwc(1:ljwc(i),i),icc1)';
      end
    else
      [jwc,ljwc] = ffill(jwc,ljwc,[]);
      temp = icc1;
    end
  else
    temp = temp(lj1(i)+1:size(temp,1),:)-ny;
    if lky(i+ykmax_)>0 then
      [jwc,ljwc] = ffill(jwc,ljwc,[temp;ky(1:lky(i+ykmax_),i+ykmax_)+(ykmax_-1)*ny]);
    else
      [jwc,ljwc] = ffill(jwc,ljwc,temp);
    end
     
    temp = union(jwc(1:ljwc(i),i),icc1)';
  end
end
 
[j1,lj1] = ffill(j1,lj1,selif(temp+ykmin_*ny,temp<=ny));
ltemp = zeros(ykmin_,1);
jwc1 = zeros(ncc1,ykmin_);
 
for i = 1:ykmin_
  temp = union(jwc(1:ljwc(i),i),icc1)';
  ltemp(i) = size(temp,1)
  if ljwc(i)>0 then
    jwc(1:ljwc(i),i) = indnv(jwc(1:ljwc(i),i),temp);
  end
  jwc1(:,i) = indnv(icc1,temp);
end
 
%v=getdate()
h1 = %v([1:2 6:9]);
 
dyn_disp('-----------------------------------------------------');
dyn_disp('MODEL SIMULATION');
dyn_disp('');
 
for iter = 1:maxit_
  %v1=getdate()
  h2 = %v1([1:2 6:9]);
  y_ = y_(:);
  err_f = 0;
   
  [fid,%v1] = mopen('dynare.swp','w+b',0)
  if %v1<0 then fid = -1;end
   
  it_ = 1+ykmin_;
  ic = 1:ny;
  iyr = iyr0;
  i = ykmin_+1;
  while (i>1)&(it_<=(iter_+ykmin_)) then
    %v2=getdate()
    h3 = %v2([1:2 6:9]);
    if broyden_&(iter>1) then
      d1_ = -evalstr(fname_+'_ff(y_(iyr))');
    else
      jacob(fname_+'_ff',y_(iyr));
      d1_ = -d1_;
    end
     
    err_f = max(err_f,max(abs(d1_)));
    if lky(i)~=0 then
      j1i = ky(1:lky(i),i);
      w0 = jacobia_(:,isc(i-1)+1:isc(i));
    else
      w0 = [];
    end
    ttemp = iy(i+1:i+ykmax_,:)';
    jwci = find(ttemp)';
    if ~(jwci==[]) then
      w = jacobia_(:,isc(i)+1:isc(i+ykmax_));
    end
    j = i;
    while j<=ykmin_ then
      if ~(w0==[]) then
         
        ofs = (it_-ykmin_-ykmin_+j-2)*ny*ncc*8;
        mseek(ofs,fid,'set');
	c = matrix(mget(ncc*ny,'d',fid)',ny,ncc);
         
	if jwci==[] then
          w = -w0*c(j1i,1:ncc1);
          jwci = icc1;
        else
          iz = union(jwci,icc1)';
          ix = indnv(jwci,iz);
          iy__ = indnv(icc1,iz);
          temp = zeros(size(w,1),size(iz,1));
          temp(:,ix) = w;
          temp(:,iy__) = temp(:,iy__)-w0*c(j1i,1:ncc1);
          w = temp;
          jwci = iz;
          clear('temp','iz','ix','iy__');
	end
        d1_ = d1_-w0*c(j1i,ncc);
        clear('c');
      end
      j = j+1;
      if jwci==[] then
        j1i = [];
        if lky(j+ykmax_)~=0 then
          jwci = ky(1:lky(j+ykmax_),j+ykmax_)+(ykmax_-1)*ny;
          w = jacobia_(:,isc(j+ykmax_-1)+1:isc(j+ykmax_));
        else
          jwci = [];
        end
      else
        j1i = selif(jwci,jwci<(ny+1));
        w0 = w(:,1:size(j1i,1));
        if size(jwci,1)==size(j1i,1) then
          if lky(j+ykmax_)~=0 then
            jwci = ky(1:lky(j+ykmax_),j+ykmax_)+(ykmax_-1)*ny;
            w = jacobia_(:,isc(j+ykmax_-1)+1:isc(j+ykmax_));
          else
            jwci = [];
          end
        else
          jwci = jwci(size(j1i,1)+1:size(jwci,1),:)-ny;
          w = w(:,size(j1i,1)+1:size(w,2));
          if lky(j+ykmax_)~=0 then
            jwci = [jwci;ky(1:lky(j+ykmax_),j+ykmax_)+(ykmax_-1)*ny];
            w = [w,jacobia_(:,isc(j+ykmax_-1)+1:isc(j+ykmax_))];
            //   else
            //     jwci = [] ;
          end
        end
      end
    end
    dyn_disp(jwci)
    dyn_disp(icc1)
    jwci = [indnv(jwci,icc1);ncc];
    w = [w,d1_];
    c = zeros(ny,ncc);
    dyn_disp(jwci)
    c(:,jwci) = w0\w;
    clear('w','w0');
    
  junk = mseek(0,fid,'end');
  
  mput(c,'d',fid);
     
  it_ = it_+1;
  ic = ic+ny;
  iyr = iyr+ny;
  i = i-1;
end
icr0 = (it_-ykmin_-ykmin_-1)*ny;
while it_<=(iter_+ykmin_) then
    if broyden_ then
      d1_ = -evalstr(fname_+'_ff(y_(iyr))');
    else
      jacob(fname_+'_ff',y_(iyr));
      d1_ = -d1_;
    end
    err_f = max(err_f,max(abs(d1_)));
    w0 = jacobia_(:,1:isc(1));
    w = jacobia_(:,isc(1)+1:isc(1+ykmax_));
    j = 1;
    while j<=ykmin_ then
      icr = j1(1:lj1(j),j)-(j-1)*ny;
       
      ofs = (icr0+(j-1)*ny+1-1)*ncc*8;
      junk = mseek(ofs,fid,'set');
      c = matrix(mget(ncc*ny,'d',fid)',ny,ncc);
       
      temp = zeros(ny,ltemp(j));
      if ljwc(j)>0 then
        temp(:,jwc(1:ljwc(j),j)) = w;
      end
      temp(:,jwc1(:,j)) = temp(:,jwc1(:,j))-w0*c(icr,1:ncc1);
      w = temp;
      clear('temp');
      d1_ = d1_-w0*c(icr,ncc);
      clear('c');
      j = j+1;
      w0 = w(:,1:lj1(j));
      if ykmax_==1 then
        w = jacobia_(:,isc(j+ykmax_-1)+1:isc(j+ykmax_));
      else
        w = w(:,lj1(j)+1:size(w,2));
         
        if lky(j+ykmax_)>0 then
          w = [w,jacobia_(:,isc(j+ykmax_-1)+1:isc(j+ykmax_))];
        end
      end
    end
    c = w0\[w,d1_];
    d1_ = [];
    clear('w','w0');
  junk = mseek(0,fid,'end');
  mput(c,'d',fid);
    clear('c');
    it_ = it_+1;
    ic = ic+ny;
    iyr = iyr+ny;
    icr0 = icr0+ny;
  end
  if ct_==1 then
     
    ofs = ((it_-ykmin_-2)*ny+1-1)*ncc*8;
    junk = mseek(ofs,fid,'set');
    c = matrix(mget(ncc*ny,'d',fid),ncc,ny);

    for i = 1:ykmax_
      w = tril(triu(ones(ny,ny+ncc1)));
      w(:,jwc1(:,ykmin_)) = w(:,jwc1(:,ykmin_))+c(:,1:ncc1);
      c = [w(:,ny+1:size(w,2))',c(:,ncc)]/w(:,1:ny);
       
      junk = mseek(0,fid,'end');
      mput(c,'d',fid);
       
      it_ = it_+1;
      ic = ic+ny;
       
    end
  end
  y_ = matrix(y_,ny,iter_+ykmin_+ykmax_);
  if ct_==1 then
    %v2=getdate()
    hbacsup = %v2([1:2 6:9]);
    c = bksupk(ny,fid,ncc,icc1);
    %v2=getdate()
     
    //!! Unknown function etime ,the original calling sequence is used
//    hbacsup = etime(%v2([1:2 6:9]),hbacsup);
    c = matrix(c,ny,iter_+ykmax_)';
     
    y_(:,1+ykmin_:iter_+ykmax_+ykmin_) = y_(:,1+ykmin_:iter_+ykmax_+ykmin_)+slowc_*c';
  else
    %v2=getdate()
    hbacsup = %v2([1:2 6:9]);
    c = bksupk(ny,fid,ncc,icc1);
    %v2=getdate()
     
    //!! Unknown function etime ,the original calling sequence is used
//    hbacsup = etime(%v2([1:2 6:9]),hbacsup);
    c = matrix(c,ny,iter_)';
    y_(:,1+ykmin_:iter_+ykmin_) = y_(:,1+ykmin_:iter_+ykmin_)+slowc_*c';
  end
   
  mclose(fid);
   
  %v1=getdate()
   
  //!! Unknown function etime ,the original calling sequence is used
//  h2 = etime(%v1([1:2 6:9]),h2);
   
  [junk,i1] = max(abs(c),'r');
  [junk,i2] = max(junk);

  err = max(max(abs(c ./ scalv_'),'r'));
  dyn_disp(string(iter)+'- err = '+string(err));
  dyn_disp('err_f = '+string(err_f));
  dyn_disp(' Time of this iteration : '+string(h2));
  if timing_ then
    dyn_disp(' Back substitution  : '+string(hbacsup));
  end
  if err<dynatol_ then
    %v2=getdate()
     
    //!! Unknown function etime ,the original calling sequence is used
//    h1 = etime(%v2([1:2 6:9]),h1);
    mtlb_fprintf('\n');
    dyn_disp(' Total time of simulation : '+string(h1));
    mtlb_fprintf('\n');
    dyn_disp(' Convergence achieved.');
    dyn_disp('-----------------------------------------------------');
    mtlb_fprintf('\n');
    return
     
  end
end
dyn_disp('WARNING : the maximum number of iterations is reached.');
mtlb_fprintf('\n');
dyn_disp('-----------------------------------------------------');
return
 
 
// 2/11/99 MJ took out reshapel
 
 
 
 
 
 
