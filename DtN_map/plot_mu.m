

Mu_U = zeros(2,5);
Mu_V = zeros(2,5);
MaxUV = zeros(2,5);

for ii = 5:9
    Nx = 2^ii;
    load(['test_DE/DtN_',num2str(Nx),'.mat'])

    data = flux_d_total;

    data = data(1:end/4, end/2+1:end/4*3);
    
    % size of the matrix to be recovered
    [n,m] = size(data);

    % we compute the coherence of the subspaces
    [U, S, V] = svd(data);

    threshold = 1e-6;
    s = S(S>threshold) ;
    r = length(s);

    mu_U = coherence(U(:,1:r));
    mu_V = coherence(V(:,1:r));

    fprintf(strcat("mu(U) is ", num2str(mu_U), "\n"))
    fprintf(strcat("mu(V) is ", num2str(mu_V), "\n"))
 
    Ur = U(:, 1:r);
    Vr = V(:, 1:r);
    MaxUV(1, ii-4) = max(max(abs(Ur*Vr')));
    Mu_U(1,ii-4) = mu_U;
    Mu_V(1,ii-4) = mu_V; 
end 


for ii = 5:9
    Nx = 2^ii;
    load(['test_DE/DtN_',num2str(Nx),'.mat'])

    data = flux_d_total;
    
    data = data(end/8+1:end/4, end/2+1:end/8*5);
    % size of the matrix to be recovered
    [n,m] = size(data);

    % we compute the coherence of the subspaces
    [U, S, V] = svd(data);

    threshold = 1e-6;
    s = S(S>threshold) ;
    r = length(s);

    mu_U = coherence(U(:,1:r));
    mu_V = coherence(V(:,1:r));

    fprintf(strcat("mu(U) is ", num2str(mu_U), "\n"))
    fprintf(strcat("mu(V) is ", num2str(mu_V), "\n"))
 
    Ur = U(:, 1:r);
    Vr = V(:, 1:r);
    MaxUV(2, ii-4) = max(max(abs(Ur*Vr')));
    
    Mu_U(2,ii-4) = mu_U;
    Mu_V(2,ii-4) = mu_V; 
end 



figure(1); clf();
semilogx(2.^[5:9], Mu_U, 'LineWidth', 3); 
hold on; 
semilogx(2.^[5:9], Mu_V, '--', 'LineWidth', 3); 
legend('a) $\mu_U$','b) $\mu_U$','a) $\mu_V$','b) $\mu_V$','Interpreter','latex')
set(gca,'FontSize',24)

xlabel('n')
ylabel('coherence')


figure(2); clf();
loglog(2.^[5:9], MaxUV, 'LineWidth', 3); 
hold on; 
loglog(2.^[5:9], 8./(2.^[5:9]), '--', 'LineWidth', 3); 
legend('a) $max(UV^*)$','b) $max(UV^*)$','$\qquad \mathcal{O}(r/n)$','Interpreter','latex')
set(gca,'FontSize',24)


xlabel('n')
ylabel('Max value UV^*')
