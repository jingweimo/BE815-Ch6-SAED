%Multivariate simulations: mvnrnd.m

%2-D data
mu = [2 3];
sigma = [1 1.5; 1.5 3];
rng default  % For reproducibility
R = mvnrnd(mu,sigma,100);
figure
plot(R(:,1),R(:,2),'o');
xlabel('x');
ylabel('y');
set(gca,'xTick',[],'yTick',[]);
box off;
% %lsFitting
% [beta,R2,nrmse,yhat] = lsFitting(R(:,1),R(:,2));
% hold on;
% plot(R(:,1),yhat,'r-');
% text(0.1,0.9,['R^{2}=',num2str(R2,'%1.2f')],'Units','normalized')
% text(0.1,0.8,['NRMSE=',num2str(100*nrmse,'%2.1f'),'%'],'Units','normalized')

%3-D data 
mu = [2 3 1];
sigma = [1 0.5 2;0.5 3 0;2 0 5]; %[vec,val] =eig(sigma)
rng default  % For reproducibility
R = mvnrnd(mu,sigma,100);
figure
plot3(R(:,1),R(:,2),R(:,3),'o');
xlabel('x1');
ylabel('x2');
zlabel('x3');
set(gca,'xTick',[],'yTick',[],'zTick',[]);
box off;

%% ED vs MD 
mu = [0 0];
sigma = [1 0.8; 0.8 1];
rng default  % For reproducibility
R = mvnrnd(mu,sigma,200);
% figure
% plot(R(:,1),R(:,2),'o');
% xlabel('x1');
% ylabel('x2');
% set(gca,'xTick',[],'yTick',[]);
% box off;

%MD ellipse
md_cutoff = 3;
drawMDellipse(R,md_cutoff,'dispOutliers',0);
%ED circle
data_center = mean(R);
hold on;
rectangle('Position',[data_center-md_cutoff,2*md_cutoff,2*md_cutoff],'Curvature',[1,1],'EdgeColor','g');
xlabel('x1'); ylabel('x2');
set(gca,'xTick',[],'yTick',[]);
box off;

%% PCA applications
%1) dimension reduction for 2-D visualization
mu1 = [0.5 3 5];
sigma1 = [1 0.6 0.8;0.6 3 0.6;0.8 0.6 1]; %[vec,val] =eig(sigma1)
mu2 = [2 6 10];
sigma2 = [1 0.9 0.8;0.9 3 0.7;0.8 0.7 3]; %[vec,val] =eig(sigma2)
rng default  % For reproducibility
R1 = mvnrnd(mu1,sigma1,50);
R2 = mvnrnd(mu2,sigma2,50);
figure
plot3(R1(:,1),R1(:,2),R1(:,3),'bo'); hold on;
plot3(R2(:,1),R2(:,2),R2(:,3),'ro');
xlabel('x1');
ylabel('x2');
zlabel('x3');
set(gca,'xTick',[],'yTick',[],'zTick',[]);

%mixed together
R = cat(1,R1,R2); 
[~,T,~] = pca(R); %cumsum(L)/sum(L);
figure;
plot(T(1:50,1),T(1:50,2),'bo'); hold on;
plot(T(51:100,1),T(51:100,2),'ro'); 
xlabel('PC1');
ylabel('PC2');

%3) data diagnostics 
mu3 = [2 0 1];
sigma3 = [1 -0.5 0;-0.5 2 -0.6;0 -0.6 1]; %[vec,val] =eig(sigma1)
rng default  % For reproducibility
R3 = mvnrnd(mu3,sigma3,100);
figure
plot3(R3(:,1),R3(:,2),R3(:,3),'o'); hold on;
xlabel('x1');
ylabel('x2');
zlabel('x3');
set(gca,'xTick',[],'yTick',[],'zTick',[]);
%pca
[~,T,~] = pca(R3); %cumsum(L)/sum(L);
% figure;
% plot(T(:,1),T(:,2),'bo'); hold on;
% xlabel('PC1');
% ylabel('PC2');
%draw an error ellipse
drawMDellipse(T(:,1:2),sqrt(chi2inv(0.975,2)),'dispOutliers',1);
xlabel('PC1');
ylabel('PC2');

%calculation of score distance
sd = zeros(size(R3,1),1);
for ii = 1:length(sd)
    d = 0;
    for k = 1:2 %two PCs 
        d = d + T(ii,k)^2/L(k);
    end
    sd(ii) = sqrt(d);
end
dataIdx = 1:length(sd);
figure;
plot(dataIdx,sd,'o');
xlabel('Data points');
ylabel('Score distance');
set(gca,'xtick',[]);
hold on; 
plot(sqrt(chi2inv(0.975,2))*ones(size(sd)),'r--');
outIdx = sd>sqrt(chi2inv(0.975,2));
plot(dataIdx(outIdx),sd(outIdx),'r.');

