function H = strf_1D_timedomain(D, F, max_lag)

% function H = strf_1D_timedomain(D, F, max_lag)
% 
% Calculates a 1D STRF by regressing shifted copies of the feature vectors in matrix F against data
% vectors in D.
% 
% -- Input --
% 
% F: Number of time-points x number of features matrix
% 
% D: Number of time-points x number of data points
% 
% max_lag: Maximum lag allowed in the kernel
% 
% -- Output -- 
% 
% H: lag x number of features x number of data points
% 
% -- Example --
% 
% % feature vector
% N = 100;
% f = rand(N,1); % feature vector
% 
% % kernel
% M = 10;
% h = rand(M,1);
% 
% % convolve feature vector with kernel 
% hpad = [h; zeros(N-M,1)]; % kernel
% d = ifft(fft(hpad) .* fft(f)); % data vector
% 
% % infer kernel using strf estimation 
% h_est = strf_1D_timedomain(d, f, M);
% 
% % compare inferred and true kernel
% h'
% h_est'
% 
% 2016-02-16: Documented by Sam NH

% dimensionality of input matrix
n_d = size(D,2);
n_f = size(F,2);

H = nan(max_lag, n_f, n_d);
for j = 1:n_f
    
    Fshift = nan(length(F(:,j)), max_lag);
    for k = 1:max_lag
        Fshift(:,k) = circshift(F(:,j),k-1);
    end
    
    Finv = pinv(Fshift);
    
    for i = 1:n_d
                
        H(:,j,i) = Finv * D(:,i);
        
    end
end

H = squeeze(H);