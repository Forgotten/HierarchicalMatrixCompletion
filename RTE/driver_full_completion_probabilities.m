%driver
% load cvx 
run /home/lzepeda/Software/cvx/cvx_startup.m
%cvx_solver Mosek 


% different probabilities
p = (1:20)/20;
% number of trials per
Nit = 10;

%number of levels
Nlevls = 2;

for kk = 1:4

    Nr = 50*2^(kk-1);
    Nv = 60*2^(kk-1);
    load(['data/albedo_ep_5_freq_2_Nr_', num2str(Nr),'_Nv_',num2str(Nv),'.mat'])

    data = f;
    data = data(1:end/2, end/2+1:end);

    errApproxArray = zeros(length(p), Nit);
    errGlobalArray = zeros(length(p), Nit);
    ratioArray = zeros(length(p), Nit);
    nArray = zeros(length(p), Nit);

    % we save the data in a file 
    name = ['completion_',num2str(Nv),'.mat'];
 
    % looping over the different probabilities
    for ii= 1:length(p)
        % looping over each realization
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


end