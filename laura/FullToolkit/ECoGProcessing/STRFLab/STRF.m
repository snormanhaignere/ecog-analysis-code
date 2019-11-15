function w = STRF(x,y,lags,lambda,method,lag_flag)
% w = STRF(x,y,lags,lamda)
% Function to calculate a spectrotemporal receptive field (STRF)
% or Reconstruction Filter via linear regression.

if ~exist('lag_flag','var')
    lag_flag = 1;
end

if size(x,2) > size(x,1)
    x = x';
end
if size(y,2) > size(y,1)
    y = y';
end

if lag_flag
    x_lag = lagmatrix(x,lags); % create lag matrix
else
    x_lag = x;
    clear x
end

x_lag = [ones(size(x_lag,1),1) x_lag]; % Add constant parameter

XY = x_lag'*y; % cross covariance
XXt = x_lag'*x_lag; % auto covariance

w = zeros(length(lambda),size(x_lag,2),size(y,2));

switch method
    case 'svd'
        
        [u,d,v] = svd(XXt);
        D = diag(d);
        
        tmp = D/sum(D);
        
        for i = 1:length(lambda)
            for cnt1 = 1:length(tmp)
                if sum(tmp(1:cnt1))>lambda(i)
                    break;
                end
            end
            D2 = 1./D;
            D2(cnt1+1:end) = 0;
            RR_Inv = (v*diag(D2)*u');
            w(i,:,:) = RR_Inv*XY;
        end
        
    case 'ridge'        
        for i = 1:length(lambda)
            w(i,:,:) = (XXt + lambda(i)*eye(size(XXt,1)))\XY;
        end
end

w = squeeze(w);

end