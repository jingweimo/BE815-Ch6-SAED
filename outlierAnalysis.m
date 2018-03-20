%error_ellipse: http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
close all;

%create data
mu = [0, 0]; %pre-mean-centered
sigma = [1,0.5;0.5,1];
n = 200;
data = mvnrnd(mu,sigma,n);
y1 = data(:,1);
y2 = data(:,2);
figure; 
plot(y1,y2,'o');
mindata = min(min(data));
maxdata = max(max(data));
xlim([mindata-3, maxdata+3]);
xlim([mindata-3, maxdata+3]);
xlabel('Variable 1');
ylabel('Variable 2');
grid;

%% eig docomposition
[eigenvec, eigenval] = eig(cov(data));
dataVar = diag(eigenval);
% data variances
[var_large, idx_large] = max(dataVar);
[var_small, idx_small] = min(dataVar);
% eigen vectors
eignvec_large = eigenvec(:,idx_large);
eignvec_small = eigenvec(:,idx_small);

% chi-squre cutoff value for a 97.5% error ellipse
% note the Mahalanobis distance follows a chi-square distribution
chisquare_val = sqrt(chi2inv(0.975,2));
theta_grid = linspace(0,2*pi);

% define x and y coordinates (based on the Mahalanobis concept)
ellipse_x_r  = chisquare_val*sqrt(var_large)*cos(theta_grid);
ellipse_y_r  = chisquare_val*sqrt(var_small)*sin(theta_grid);

% ellipse orientation (used for ellipse drawing)
phi = atan2(eignvec_large(2), eignvec_large(1));
if(phi < 0)
    phi = phi + 2*pi;
end
R = [cos(phi), sin(phi); -sin(phi), cos(phi)];
%rotatation wrt the x-axis alligned ellipse
r_ellipse = [ellipse_x_r; ellipse_y_r]' * R;

% draw the error ellipse
data_center = mean(data);
x0=data_center(1);
y0=data_center(2);
hold on;
plot(r_ellipse(:,1) + x0,r_ellipse(:,2) + y0,'r-')

% plot the ellipe axes
a = chisquare_val*sqrt(var_large);
b = chisquare_val*sqrt(var_small);
diff_large = [eignvec_large(1)*a, eignvec_large(2)*a] - data_center;
diff_small = [eignvec_small(1)*b, eignvec_small(2)*b] - data_center;
quiver(x0, y0, diff_large(1),diff_large(2),0,'-k','MaxHeadSize',0.5);
quiver(x0, y0, diff_small(1),diff_small(2),0,'-k','MaxHeadSize',0.5);

%% identify outliers
for ii = 1:n
    tempX = y1(ii);
    tempY = y2(ii);
    if ((tempX-x0)*cos(phi) + (tempY-y0)*sin(phi))^2/a^2 + ((tempX-x0)*sin(phi) - (tempY-y0)*cos(phi))^2/b^2 > 1
        plot(tempX,tempY,'r.');
    end
end

%% Mahalanobis distance for outlier identification 
% MD = (x-mean(x))*inv(C)*(x-mean(x))' where each row in x represent a obvervation 
distMaha = ones(1,n);
for ii = 1:n
    tempX = y1(ii);
    tempY = y2(ii);
    tempDist = ([tempX, tempY] - data_center)/(cov(data))*([tempX, tempY] - data_center)';
    distMaha(ii) = sqrt(tempDist);
end
%matrix notation 
%temp = (data - repmat(data_center,[n,1]))/cov(data)*(data - repmat(data_center,[n,1]))';
%distMaha = sqrt(diag(temp));
pntIdx = 1:n;
figure;
plot(pntIdx, distMaha,'o'); hold on;
plot(pntIdx, chisquare_val*ones(1,n),'r--');
xlabel('Point index');
ylabel('Mahalanobis distance');

%highlight outliers
plot(pntIdx(distMaha>chisquare_val),distMaha(distMaha>chisquare_val),'r.');

