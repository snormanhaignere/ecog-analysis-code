function H = strf_2D_freqdomain(D, F, max_lag, frac_eig)

% function H = strf_2D_freqdomain(D, F, max_lag, frac_eig)
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
% frac_eig: Number of eigen-vectors of the feature matrix to use when performing the analysis
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
% frac_eig = 1;
% H_est = strf_2D_freqdomain(D, F, M, frac_eig);
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
n_t = size(D,1);

% FFTs of input matrices
D_fft = fft(D);
F_fft = fft(F);

% Correlation matrix
Cff = nan(n_f * max_lag);
for j = 1:n_f
    fprintf('%d\n',j); drawnow;
    for k = 1:j

        % cross-correlation between a pair of features
        Cff_pair = ifft(conj(F_fft(:,j)) .* F_fft(:,k));
        
        % create the toeplitz matrix from the beginning and end of the full cross-corr function
        toepmat = toeplitz(Cff_pair(1:max_lag), Cff_pair([1, n_t-(0:max_lag-2)]));
        
        % toeplitz sub matrix
        xi = (1:max_lag) + (j-1)*max_lag;
        yi = (1:max_lag) + (k-1)*max_lag;
        Cff(xi, yi) = toepmat;
        Cff(yi, xi) = toepmat';
        
    end
end

% eigen value decomposition
% of correlation matrix
[Q,S] = eig(Cff);

% sort by eigen values
eigvals = diag(S);
[~,xi] = sort(eigvals,'descend');
Q = Q(:,xi);
S = diag(eigvals(xi));
eigvals = eigvals(xi); %#ok<NASGU>

% select top N% of eigen vectors
ncomp = round(max_lag * n_f * frac_eig);
Q = Q(:,1:ncomp);
S = S(1:ncomp,1:ncomp);
invCff = Q * (S \ Q');

% feature-data correlation matrix
H = nan(max_lag * n_f, n_d);
for j = 1:n_d
    Cfd = nan(n_f * max_lag,1);
    for i = 1:n_f
        Cfd_single_feature = ifft(conj(F_fft(:,i)) .* D_fft(:,j));
        Cfd((1:max_lag) + (i-1)*max_lag) = Cfd_single_feature(1:max_lag);
    end
    H(:,j) = invCff * Cfd;
end

H = reshape(H, [max_lag, n_f, n_d]);






