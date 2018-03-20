function [beta, Rsq, NRMSE, varargout] = lsFitting(x,y)
%The comptuer the LS solution for LS (line or curve) fitting of y = ax + b
%Inputs:
%      x: dependent variable measurements (one variable considered here)
%      y: indpendent variable measurements
%Outputs:
%      beta: (a,b) %for line fitting
%      Rsq and NRSME: determination coefficient and normalized rmse
%YL, 03/16/2018

x = x(:);
y = y(:); 

A = [x, ones(length(x),1)];
beta = (A'*A)\A'*y;

yhat = A*beta;
res = y - yhat;
sse = res'*res;
sst = (y - mean(y))'*(y-mean(y));
Rsq = 1 - sse/sst;
NRMSE = sqrt(sse/length(x))/range(y); %using range(y) for mean-centered data
varargout{1} = yhat;
end

