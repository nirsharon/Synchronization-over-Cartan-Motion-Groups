% script name: "test_sync_SEk_by_SVD"
%
% We test the "sync_SEk_by_SVD" procedure
%
% N.S, March 2016

n = 80;
k = 5;
s = (k+1)*n;  % the size of matrices


%---- synthetic data in SE(k) ------
SEk_array = zeros(k+1,k+1,n);
R = 1;  % scale for transitive part data. can ignore.

for l=1:n
    A = 3*rand(k);
    [q, ~] = qr(A);
    % positive determinant
    q(:,k) = det(q)*q(:,k);
    % transational part
    b = R*rand(k,1);
    % embed in the matrix
    SEk_array(1:k,1:k,l) = q;
    SEk_array(1:k,k+1,l) = b;
    SEk_array(k+1,k+1,l) = 1;
end


p_values = .21;%.25;% 0.1:0.3:.9;  % the probabilities of non-outliers
number_p_values = numel(p_values);
error_rates = zeros(number_p_values,1);

for q=1:number_p_values
    
    % ================= Choosing randomly places for outliers ======
    p = p_values(q);
    %p = 1; %0.6;  % the probability of non-outliers
    m = n*(n-1)/2;    % full graph
    non_outliers = floor(p*m); %number of non-outliers
    y = randsample(m,non_outliers);
    
    % the outliers places
    prob_arr = sparse(n,n);
    idx = find(~tril(ones(n)));  % indices of upper side blocks
    prob_arr(idx(y))=1;          % mark only the relevant, non-outliers
    confidence_weights = eye(n)+prob_arr+prob_arr';

    % added noise. set to zero sig1 and sig2 for no noise
   parms.d = k; 
   parms.sig1 = 0;%.3; 
   parms.sig2 = 0;%.1;
   noise_func = @naive_random_SE_d;
   % another option for no noise...
   % parms = [];
   % noise_func = @(x) eye(k+1);
   Affin_mat = MakeAffinityMatrix(SEk_array, prob_arr, noise_func, parms);

    %---- calling the functions -----
    estimations = sync_SEk_by_SVD( triu(Affin_mat), confidence_weights, k );
    error_rates(q) = error_calc_SE_k( estimations, SEk_array );
end
error_rates'
if numel(error_rates)>3
    plot(p_values,error_rates,'r','LineWidth',3.5);
end
% drawing the spectrum if needed
% spec = eig(Affin_mat); h1 = histogram(abs(spec),n);

