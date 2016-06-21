function H = strf_2D_timedomain(D, F, max_lag)

% function H = strf_2D_timedomain(D, F, max_lag)
% 
% Calculates a 2D STRF by regressing shifted copies of the feature vectors in matrix F against data
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
% P = 3;
% F = rand(N,P); % feature vector
% 
% % kernel
% M = 10;
% H = rand(M,P);
% 
% % convolve feature vector with kernel 
% Hpad = [H; zeros(N-M,P)]; % kernel
% D = zeros(N,1);
% for i = 1:P
%     D = ifft(fft(Hpad(:,i)) .* fft(F(:,i))) + D; % data vector
% end
% 
% % infer kernel using strf estimation 
% H_est = strf_2D_timedomain(D, F, M);
% 
% % compare inferred and true kernel
% subplot(2,1,1);
% imagesc(H')
% subplot(2,1,2);
% imagesc(H_est')
% 
% 2016-02-16: Documented by Sam NH

% dimensionality of input matrix
n_d = size(D,2);
n_f = size(F,2);
n_t = size(F,1);

% created shifted copies of each column of F
Fshift = nan(n_t, max_lag, n_f);
for j = 1:n_f
    for k = 1:max_lag
        Fshift(:,k,j) = circshift(F(:,j),k-1);
    end
end

% unwrap all lags
Fshift = reshape(Fshift, [n_t, max_lag * n_f]);
Finv = pinv(Fshift);

% estimate the strf
H = nan(max_lag * n_f, n_d);
for i = 1:n_d
    H(:,i) = Finv * D(:,i);
end

% reshape to separate lags and features
H = reshape(H, [max_lag, n_f, n_d]);