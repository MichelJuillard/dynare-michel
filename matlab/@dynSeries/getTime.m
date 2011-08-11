function time = getTime(ts)
    
if ts.freq==1
    time = ts.time(:,1);
    return
end

time = [];

switch ts.freq
  case 4
    for i=1:ts.nobs
        time = char(time,[num2str(ts.time(i,1)) 'Q' num2str(ts.time(i,2))]);
    end
  case 12
    for i=1:ts.nobs
        time = char(time,[num2str(ts.time(i,1)) 'M' num2str(ts.time(i,2))]);
    end
  case 52
    for i=1:ts.nobs
        time = char(time,[num2str(ts.time(i,1)) 'W' num2str(ts.time(i,2))]);
    end
  otherwise
    error('dynSeries::getTime: Unknown type of frequency!')
end

time = time(2:end,:);