Nx = 64;
load(['test_DE/DtN_',num2str(Nx),'.mat'])

data = flux_d_total;

plr_matrix = PLR_matrix(data, 5, 0.000001, 'l2');

[hist_size szs nums mydisp] = plr_matrix.histogram('draw');
set(gca,'FontSize',48)

% removing ticks and labels
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])