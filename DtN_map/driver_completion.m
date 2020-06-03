%driver
% load cvx 
run ../../../Software/cvx/cvx_setup.m
cvx_solver Mosek 


% different probabilities
p = (1:10)/10;
% number of trials per
Nit = 1;

Nx = 64;
load(['test_DE/DtN_',num2str(Nx),'.mat'])

data = flux_d_total;

% we take the off-diagonal block
data = data(1:end/4, end/2+1:end/4*3);

errApproxArray = zeros(length(p), Nit);
errGlobalArray = zeros(length(p), Nit);
ratioArray = zeros(length(p), Nit);
nArray = zeros(length(p), Nit);

% we save the data in a file 
name = ['completion_',num2str(Nx),'.mat'];
 

for ii= 1:length(p)
    for jj = 1:Nit
        [errApprox, errGlobal, ratio, n, m] = completion(data, p(ii));
    
        errApproxArray(ii, jj) = errApprox;
        errGlobalArray(ii, jj) = errGlobal;
        ratioArray(ii, jj) = ratio;
        nArray(ii, jj) = n;
        
    end
    % we save at each iteration
    save(name, 'errApproxArray', 'errGlobalArray', 'ratioArray', 'nArray')
end


