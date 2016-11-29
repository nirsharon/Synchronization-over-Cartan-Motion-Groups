function [errors, avg_snr] = MakeFigure_MissingDataWvaryNoise(n, d)
% Ground for comparison: Full graph, ourlier, fixed noise level
% Compared methods     : Spectral (SVD), contraction, ASAP
%
%
% NS, June 16

% initialization
s          = (d+1)*n;               % the size
AvailData  = 0.25;
if n<150                            % the different noise levels
    noiseValues    = .6:-0.08:.28;%.3;        
else
    noiseValues    = .8:-0.08:.4;%.3;        
end

iterations = numel(noiseValues);
errors     = zeros(iterations,2);   % initiliaze errors array
repeatedIters = 3;                  % repeating trials

% generate ground truth (synthetic) data
GT_data = generate_data_SE_k(n, d, 2);
avg_snr = zeros(iterations,1);
% main loop
for q=1:iterations
    % setting measurements entries (non-outliers/outliers)
    p = AvailData; 
    m = n*(n-1)/2;               % full graph (upper/lower side)
    non_outliers = floor(p*m);   % number of nonoutliers
    y = randsample(m,non_outliers);
    [i,j] = find(triu(ones(n),1));  % indices of upper side blocks
    I = i(y); J=j(y);
    prob_arr = sparse(I,J,ones(numel(I),1),n,n);
    
    confidence_weights = eye(n)+prob_arr+prob_arr';
    
    parms.d    = d; 
    parms.sig1 = noiseValues(q); 
    parms.sig2 = noiseValues(q); 
    noise_func = @naive_random_SE_d;                 
    
    % calling the functions    
    SPEC_error = zeros(repeatedIters,1);
    CONT_error = zeros(repeatedIters,1);
    ASAP_error = zeros(repeatedIters,1);
    snr_level = zeros(repeatedIters,1);

    for r=1:repeatedIters
        [Affin_mat, snr_level(r)] = MakeAffinityMatrix(GT_data, prob_arr, noise_func, parms, 0);
        % Spectral method based on SVD
        estimations   = sync_SEk_by_SVD( triu(Affin_mat), confidence_weights, d );
        SPEC_error(r) = error_calc_SE_k( estimations, GT_data );
        % Contraction method
        lambda = 200;                                         % TO DO: SMART Lambda CHOICE
        estimations2  =  SyncSEbyContraction(Affin_mat, confidence_weights, d, lambda);
        CONT_error(r) = error_calc_SE_k( estimations2, GT_data );       
        %     Matrix contraction (SVD)                                        
        %     estimations3 = sync_SEk_by_PD_contraction( triu(Affin_mat), confidence_weights, d, lambda );
        %     errors(q,3) = error_calc_SE_k( estimations3, SEk_GT_data );
        
        estimations4 = sync_SEk_by_ASAP( triu(Affin_mat), confidence_weights, d );
        ASAP_error(r)  = error_calc_SE_k( estimations4, GT_data );
    end
    errors(q,1) = mean(SPEC_error);
    errors(q,2) = mean(CONT_error);
    errors(q,3) = mean(ASAP_error);
    avg_snr(q) = mean(snr_level);
end

%avg_snr
figure;
hold on;
plot(avg_snr, errors(:,3),'--k','LineWidth',4);
plot(avg_snr, errors(:,1),'-.','Color','red','LineWidth',4.7);
plot(avg_snr, errors(:,2),'Color','blue','LineWidth',5);
legend('Separation','Spectral','Contraction');%,'Location','northwest');
fn = sprintf('n = %d, d = %d', n, d);
title(fn);
xlabel('SNR');
ylabel('MSE');
hold off;
end

