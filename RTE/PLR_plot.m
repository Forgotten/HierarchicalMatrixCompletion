Nr = 100;
Nv = 120;

load(['data/albedo_ep_0_freq_2_Nr_', num2str(Nr),'_Nv_',num2str(Nv),'.mat'])

data = f(:, end:-1:1);

plr_matrix = PLR_matrix(data, 5, 0.00001, 'l2');

[hist_size szs nums mydisp] = plr_matrix.histogram('draw');
set(gca,'FontSize',24)