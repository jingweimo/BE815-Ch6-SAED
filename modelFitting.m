%Ordinary least square (OLS) linear model fitting
%textbook (2nd edition) section 6.6.2, example 6.17

%input and output data
L = [0, 0.5, 1, 1.5, 2, 2.5]';
V = [0.05, 0.52, 1.03, 1.50, 2.00, 2.56]';
figure;
plot(L,V,'b+');
xlabel('L');
ylabel('V');

%% LS fitting
%1) linear model: V = a*L + b = [L, 1]*[a,b]' = A*x'
%OLS solution = inv(A'*A)*A'*V
A = [L,ones(length(L),1)];
%x = inv(A'*A)*A'*V;
x = (A'*A)\A'*V;
hold on;
plot(L,A*x,'r-');
text(0.1,0.9,['V=',num2str(x(1)),'*L + ',num2str(x(2))],'Units','normalized','FontSize',12')

%2) qr factorization:A =q*r
[q,r] = qr(A);
r = r(1:2,1:2);
temp = q'*V;
x = r\temp(1:2);

%3)cftool: y = ax+b

%--------------------------------------------------------------------------
%4)determination coefficient
res = A*x - V; %residual
sse = res'*res; %sum of squares of error
sst = (V-mean(V))'*(V-mean(V)); %total sum of squres 
Rsq = 1 - sse/sst;

%5)confidence intervals of parameters:[beta-t*se, beta+t*se]
ssx = (L-mean(L))'*(L-mean(L));
df = length(L) - 2;
mse = sse/df; 
deltaParam1 = abs(tinv(0.025,df)*sqrt(mse/ssx));
deltaParam2 = abs(tinv(0.025,df)*sqrt(mse*sum(L.^2)/length(L)/ssx));
param1CI = [x(1)-deltaParam1,x(1)+deltaParam1];
param2CI = [x(2)-deltaParam2,x(2)+deltaParam2];
%compare the CIs with those obtained by cftool

%6)confidence and prediction intervals 
ypred = A*x;
ypred_CI_low = ypred - abs(tinv(0.025,df)*sqrt(mse)*sqrt(1/length(L) + ((L-mean(L)).^2)/ssx));
ypred_CI_high = ypred + abs(tinv(0.025,df)*sqrt(mse)*sqrt(1/length(L) + ((L-mean(L)).^2)/ssx));
ypred_PI_low = ypred - abs(tinv(0.025,df)*sqrt(mse)*sqrt(1 + 1/length(L) + ((L-mean(L)).^2)/ssx));
ypred_PI_high = ypred + abs(tinv(0.025,df)*sqrt(mse)*sqrt(1 + 1/length(L) + ((L-mean(L)).^2)/ssx));
figure(2);
plot(L,V,'b+','MarkerSize',7); hold on;
plot(L,A*x,'r-');
plot(L,ypred_CI_low,'g--',L,ypred_CI_high,'g--');
plot(L,ypred_PI_low,'k--',L,ypred_PI_high,'k--');

%% Model performance comparisons
%1) high-covariance data
mu = [0 0];
sigma = [1 0.95; 0.95 1];
rng default  % For reproducibility
R = mvnrnd(mu,sigma,100);
figure
plot(R(:,1),R(:,2),'o');
xlabel('x','FontSize',12);
ylabel('y','FontSize',12);
set(gca,'xTick',[],'yTick',[]);
box off;
%lsFitting
[beta,R2,nrmse,yhat] = lsFitting(R(:,1),R(:,2));
hold on;
plot(R(:,1),yhat,'r-');
text(0.1,0.9,['R^{2}=',num2str(R2,'%1.2f')],'Units','normalized')
text(0.1,0.8,['NRMSE=',num2str(100*nrmse,'%2.1f'),'%'],'Units','normalized')

%2) add artificial outliers 
mu2 = [0.2 1.5];
sigma2 = [1 0.35; 0.35 1];
rng default  % For reproducibility
R2 = mvnrnd(mu2,sigma2,15);
%data overlay
% figure
% plot(R(:,1),R(:,2),'o');
% xlabel('x','FontSize',12);
% ylabel('y','FontSize',12);
% set(gca,'xTick',[],'yTick',[]);
% box off; hold on;
% plot(R2(:,1),R2(:,2),'ro');
Rnew = cat(1,R,R2);
[~,R2New,nrmseNew,yhatNew] = lsFitting(Rnew(:,1),Rnew(:,2));
figure
plot(Rnew(:,1),Rnew(:,2),'o');
xlabel('x','FontSize',12);
ylabel('y','FontSize',12);
set(gca,'xTick',[],'yTick',[]);
box off; hold on;
plot(Rnew(:,1),yhatNew,'r-');
text(0.1,0.9,['R^{2}=',num2str(R2New,'%1.2f')],'Units','normalized')
text(0.1,0.8,['NRMSE=',num2str(100*nrmseNew,'%2.1f'),'%'],'Units','normalized')

%3) MD outlier detection
n = size(Rnew,1);
chisquare_val = sqrt(chi2inv(0.95,2)); %95% confidence
temp = (Rnew - repmat(mean(Rnew,1),[n,1]))/cov(Rnew)*(Rnew - repmat(mean(Rnew,1),[n,1]))';
distMaha = sqrt(diag(temp));
pntIdx = 1:n;
figure;
plot(pntIdx, distMaha,'o'); hold on;
plot(pntIdx, chisquare_val*ones(1,n),'r--');
xlabel('Data point');
ylabel('Mahalanobis distance');
set(gca,'xTick',[]);
%highlight outliers 
plot(pntIdx(distMaha>chisquare_val),distMaha(distMaha>chisquare_val),'r.');

%4) re LS fitting
newIdx = distMaha<=chisquare_val;
[beta,R2,nrmse,yhat] = lsFitting(Rnew(newIdx,1),Rnew(newIdx,2));
figure
plot(Rnew(newIdx,1),Rnew(newIdx,2),'o');
xlabel('x','FontSize',12);
ylabel('y','FontSize',12);
set(gca,'xTick',[],'yTick',[]);
box off;
hold on;
plot(Rnew(newIdx,1),yhat,'r-');
text(0.1,0.9,['R^{2}=',num2str(R2,'%1.2f')],'Units','normalized')
text(0.1,0.8,['NRMSE=',num2str(100*nrmse,'%2.1f'),'%'],'Units','normalized')

