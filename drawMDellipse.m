function varargout = drawMDellipse(data,chisquare_val,varargin)
%Draw an error ellipse based on MD distances from the center of 2-D data
% Inputs: m*2 matrix in which each column corresponds to one variable
% chisquare_value: cut-off value for define the ellipse size, and generally,
%                  chisquare_val = sqrt(chi2inv(0.975,2))
%          
%03/16/2018

%eig decomposition
[eigenvec, eigenval] = eig(cov(data));
dataVar = diag(eigenval);

%Variance and vectors along the major and minor axies
[var_large, idx_large] = max(dataVar);
[var_small, idx_small] = min(dataVar);
eignvec_large = eigenvec(:,idx_large);
%eignvec_small = eigenvec(:,idx_small);

%define coordinates 
theta_grid = linspace(0,2*pi);
ellipse_x_r  = chisquare_val*sqrt(var_large)*cos(theta_grid);
ellipse_y_r  = chisquare_val*sqrt(var_small)*sin(theta_grid);

%ellipse orientation
phi = atan2(eignvec_large(2), eignvec_large(1));
if(phi < 0)
    phi = phi + 2*pi;
end
R = [cos(phi), sin(phi); -sin(phi), cos(phi)];
r_ellipse = [ellipse_x_r; ellipse_y_r]' * R;

%drawing
figure;
plot(data(:,1),data(:,2),'o');
grid on;
hold on;

data_center = mean(data);
x0=data_center(1);
y0=data_center(2);
hold on;
plot(r_ellipse(:,1) + x0,r_ellipse(:,2) + y0,'r-')

%identify outliers
ip = inputParser;
ip.addOptional('dispOutliers',0, @(n) validateattributes(n,{'double'},{'integer','real','scalar'}));
ip.parse(varargin{:});
res = ip.Results;
if res.dispOutliers == 1
    idx = [];
    y1 = data(:,1);
    y2 = data(:,2);
    a = chisquare_val*sqrt(var_large);
    b = chisquare_val*sqrt(var_small);
    for ii = 1:size(data,1)
        tempX = y1(ii);
        tempY = y2(ii);
        if ((tempX-x0)*cos(phi) + (tempY-y0)*sin(phi))^2/a^2 + ((tempX-x0)*sin(phi) - (tempY-y0)*cos(phi))^2/b^2 > 1
            idx = [idx; ii];
            plot(tempX,tempY,'r.');
        end
    end
varargout{1} = idx;
end
end

