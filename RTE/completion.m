function [errApprox, errGlobal, ratio, n, m] = completion(data, p)
% data is the matrix to approximate
% p is the probability that an entry is used 


% size of the 
[n,m] = size(data);

% computing the mask
E = (rand(n,m) < p);

% indices of the mask
indE = find(E);
ratio = length(indE)/(n*m);

fprintf(strcat("number of data points is ", num2str(length(indE))))

% we start the optimiziation
cvx_precision best
cvx_begin
variable Xe(n,m)
minimize norm_nuc(Xe)
subject to
Xe(indE) == data(indE)
cvx_end

% disp(num2str(norm_nuc(Xe)))

errApprox = norm(Xe(indE)-data(indE))/norm(data(indE));
errGlobal = norm(Xe(:)-data(:))/norm(data(:));
