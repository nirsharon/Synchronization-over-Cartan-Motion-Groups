function [ estimations ] = Rounding_by_LS(U)
% We estimate in least-squares fashion the right rounding


[ln,l] = size(U);
k = l-1;
n = ln/l;

B = zeros(n,l*n);
for i=1:n
    B(i,l*i) = 1;
end

A = B*U;

%least squares estimations
[nu, s, v] = svd(A);  
alpha =v(:,2:end);                           % exact alpha = null(A);
beta = (1/s(1))*v(:,1)*nu(:,1)'*ones(n,1);   % exact pinv(A)*ones(n,1);


initial_estimations = [U*alpha,U*beta];
if det(initial_estimations(1:k,1:k))<0
    alpha(:,1) = alpha(:,1)*(-1);
    initial_estimations = [U*alpha,U*beta];
end


estimations = projecting_SE(initial_estimations, k); %, Affinity_mat, confidence_weights)


end

