function out = ThompsonOutliers(in)
%Find the single-variable outliers by the Thompson's rule
%YL 03/25/2018

%Thompson's value
tau = [1.15,1.393,1.572,1.656,1.711,1.749,1.777,1.798,1.815,1.829,1.840,1.849,...
    1.858,1.865,1.871,1.876,1.881,1.885,1.889,1.893,1.896,1.899,1.902,1.904,1.906,1.908,...
    1.910,1.911,1.913,1.914,1.916,1.917,1.919,1.920,1.921,1.922,1.923,1.924];
%or using the analytical formula for the tau values
%tau = tinv(0.975,N-2)*(N-1)/sqrt(N)/sqrt(N-2+tinv(0.975,N-2)^2)

out = [];
N = length(in);
if N>=3 && N<=40
    done = false;
    while ~done
        mu = mean(in);
        sd = std(in);
        ind = find(abs(in-mu)==max(abs(in-mu)),1); %one at a time
        if abs(in(ind)-mu)>tau(N-2)*sd
            out = [out, in(ind)];
            %updating
            in(ind) = [];
            N = length(in);
            if N<3
                break;
            end           
        else
            done = true;
        end
    end
else
end
end

