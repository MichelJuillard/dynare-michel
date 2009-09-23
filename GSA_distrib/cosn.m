function co = cosn(H);

% function co = cosn(H);
% computes the cosine of the angle between the H(:,1) and its
% projection onto the span of H(:,2:end)

% Not the same as multiple correlation coefficient since the means are not
% zero

y = H(:,1);
X = H(:,2:end);

% y = H(:,1);
% X = H(:,2:end);

yhat =  X*(X\y);
co = y'*yhat/sqrt((y'*y)*(yhat'*yhat));

yhat =  X*(X\y);
co = y'*yhat/sqrt((y'*y)*(yhat'*yhat));



