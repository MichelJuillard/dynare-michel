function y=dy_date(year,period)
  global M_
  
  y = M_.freq*(year-M_.start_date(1))+period-M_.start_date(2)+1;
  
  