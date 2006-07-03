function [x]=bseastr(s1,s2)
x=[];
// Copyright (C) 2001 Michel Juillard
// 
 
m = size(s1,1);
x = zeros(m,1);
 
//!! Unknown function deblank ,the original calling sequence is used
s1 = convstr(deblank(s1),'u');
 
//!! Unknown function deblank ,the original calling sequence is used
s2 = convstr(deblank(s2),'u');
 
for im = 1:m
  key = part(s1(im),:);
  h = size(s2,1);
  l = 1;
  while l<=h then
    mid = round((h+l)/2);
    temp = part(s2(mid),:);
    if ~(key==temp) then
      for i = 1:min(length(key),length(temp))
        if part(temp,i)>part(key,i) then
          h = mid-1;
          break
           ;
           
        else
          l = mid+1;
          break
           ;
           
        end
      end
    else
      x(im) = mid
      break
       ;
       
    end
  end
end
 
