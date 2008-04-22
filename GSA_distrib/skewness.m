function s=skewness(y)
    s=sum((y-mean(y).^3));