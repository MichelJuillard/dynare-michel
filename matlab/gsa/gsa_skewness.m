function s=gsa_skewness(y),

% y=stand_(y);
% s=mean(y.^3);
    m2=mean((y-mean(y)).^2);
    m3=mean((y-mean(y)).^3);
    s=m3/m2^1.5;