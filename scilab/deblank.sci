function strout = deblank(str)
  // Function removes trailing blanks from input string
  index = find(ascii(str) ~= 32);
  strout = part(str,1:index($));
endfunction
