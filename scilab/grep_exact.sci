function [row]=grep_exact(str1,str2)
row = [];
r = grep(str1,str2);
for i=1:length(r)
  if str1(r(i)) == str2 then
    row = r(i);
    break
  end
end
