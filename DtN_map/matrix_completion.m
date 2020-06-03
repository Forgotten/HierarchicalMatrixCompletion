Nx = 512;
load(['test_DE/DtN_',num2str(Nx),'.mat'])

data = flux_d_total;

% flagg to only use matrix completion in an off-diagonal block
offDiagFlag = 1;

if offDiagFlag == 1
    data = data(1:end/4, end/2+1:end/4*3);
end

% flag to activate a low rank restriction
lowRankFlag = 0;
rankMax = 1;

if lowRankFlag == 1
    [U, S, V] = svd(data);
    data = U(:, 1:rankMax)*(S(1:rankMax,1:rankMax)*V(:,1:rankMax).');
end

% probability of succes for the distribution
p = 0.5;

% size of the 
[n,m] = size(data);

% computing the mask
E = (rand(n,m) < p);

% indices of the mask
indE = find(E);
fprintf(strcat("number of data points is ", num2str(length(indE))))

% we start the optimiziation
cvx_precision best
cvx_begin
variable Xe(n,m)
minimize norm_nuc(Xe)
subject to
Xe(indE) == data(indE)
cvx_end

fprintf(['Error on the constraint :', num2str(norm(Xe(indE)-data(indE))), '\n'])
fprintf(['Error on the reconstruction :', num2str(norm(Xe(:)-data(:))), '\n'])
