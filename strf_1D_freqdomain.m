function H = strf_1D_freqdomain(D, F, max_lag, frac_eig)

% function H = strf_1D_freqdomain(D, F, max_lag, frac_eig)
% 
% Calculates a 1D STRF using frequency-domain operations
% (compare with strf_1D_timedomain.m) 
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
% frac_eig = 1;
% h_est = strf_1D_freqdomain(d, f, M, frac_eig);
% 
% % compare inferred and true kernel
% h'
% h_est'
% 
% 2016-02-16: Documented by Sam NH

% dimensionality of input matrix
n_d = size(D,2);
n_f = size(F,2);

% FFTs of input matrices
D_fft = fft(D);
F_fft = fft(F);

H = nan(max_lag, n_f, n_d);
for j = 1:n_f
    
    % feature autocorrelation
    Cff = ifft(abs(F_fft(:,j)).^2);
    Cff = Cff(1:max_lag);
    
    % number of components
    % eigen value decomposition
    [Q,S] = eig(toeplitz(Cff));
    
    % sort by eigen values
    eigvals = diag(S);
    [~,xi] = sort(eigvals,'descend');
    Q = Q(:,xi);
    S = diag(eigvals(xi));
    eigvals = eigvals(xi); %#ok<NASGU>
    
    % select top N% of eigen vectors
    ncomp = round(max_lag * frac_eig);
    Q = Q(:,1:ncomp);
    S = S(1:ncomp,1:ncomp);
    invCfd = Q * (S \ Q');
    
    for i = 1:n_d
        
        % cross-correlation of feature and data
        Cfd = ifft(conj(F_fft(:,j)) .* D_fft(:,i)); % Dft of the cross-correlation function
        Cfd = Cfd(1:max_lag);
        
        % invert
        H(:,j,i) = invCfd * Cfd;
        
    end
end

H = squeeze(H);