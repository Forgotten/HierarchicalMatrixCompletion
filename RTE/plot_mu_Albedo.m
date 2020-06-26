n_trials = 4; 
Mu_U = zeros(2,n_trials);
Mu_V = zeros(2,n_trials);
MaxUV = zeros(2,n_trials);

for ii = 1:4
    Nr = 50*2^(ii-1);
    Nv = 60*2^(ii-1);
    load(['data/albedo_ep_5_freq_2_Nr_', num2str(Nr),'_Nv_',num2str(Nv),'.mat'])

    data = f;

    data = data(1:end/4, end/2+1:end*3/4);
    
    % size of the matrix to be recovered
    [n,m] = size(data);

    % we compute the coherence of the subspaces
    [U, S, V] = svd(data);

    threshold = 1e-3;
    s = S(S>threshold) ;
    r = length(s);

    mu_U = coherence(U(:,1:r));
    mu_V = coherence(V(:,1:r));

    fprintf(strcat("mu(U) is ", num2str(mu_U), "\n"))
    fprintf(strcat("mu(V) is ", num2str(mu_V), "\n"))
 
    Ur = U(:, 1:r);
    Vr = V(:, 1:r);
    MaxUV(1, ii) = max(max(abs(Ur*Vr')));
    Mu_U(1,ii) = mu_U;
    Mu_V(1,ii) = mu_V; 
end 


for ii = 1:4
    Nr = 50*2^(ii-1);
    Nv = 60*2^(ii-1);
    load(['data/albedo_ep_5_freq_2_Nr_', num2str(Nr),'_Nv_',num2str(Nv),'.mat'])

    data = f;

    data = data(1:end/2, 1:end/2);
    % size of the matrix to be recovered
    [n,m] = size(data);

    % we compute the coherence of the subspaces
    [U, S, V] = svd(data);

    threshold = 1e-3;
    s = S(S>threshold) ;
    r = length(s);

    mu_U = coherence(U(:,1:r));
    mu_V = coherence(V(:,1:r));

    fprintf(strcat("mu(U) is ", num2str(mu_U), "\n"))
    fprintf(strcat("mu(V) is ", num2str(mu_V), "\n"))
 
    Ur = U(:, 1:r);
    Vr = V(:, 1:r);
    MaxUV(2, ii) = max(max(abs(Ur*Vr')));
    
    Mu_U(2,ii) = mu_U;
    Mu_V(2,ii) = mu_V; 
end 



figure(1); clf();
semilogx(60*2.^[0:3], Mu_U, 'LineWidth', 3); 
hold on; 
semilogx(60*2.^[0:3], Mu_V, '--', 'LineWidth', 3); 
legend('$a\,\, \mu_U$','$b\,\, \mu_U$','$a\,\, \mu_V$','$b\,\, \mu_V$','Interpreter','latex', 'FontSize',32)
set(gca,'FontSize',22)

xlabel('n')
ylabel('coherence')


figure(2); clf();
loglog(60*2.^[0:3], MaxUV, 'LineWidth', 3); 
hold on; 
loglog(60*2.^[0:3], 8./(60*2.^[1:4]), '--', 'LineWidth', 3); 
legend('$a\,\, {\max}(UV^*)$','$b\,\, \max(UV^*)$','$\qquad \mathcal{O}(r/n)$','Interpreter','latex', 'FontSize',32)
set(gca,'FontSize',24)


xlabel('n')
ylabel('Max value UV^*')