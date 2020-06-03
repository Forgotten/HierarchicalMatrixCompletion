epsilon_min = 1e-4

load data_completion/completion_32.mat

err_32 = sum(errGlobalArray < epsilon_min, 2)/size(errGlobalArray,2);

load data_completion/completion_64.mat

err_64 = sum(errGlobalArray < epsilon_min, 2)/size(errGlobalArray,2);

load data_completion/completion_128.mat

err_128 = sum(errGlobalArray < epsilon_min, 2)/size(errGlobalArray,2);

load data_completion/completion_512.mat

err_512 = ones(10,1);
err_512(1:3) = sum(errGlobalArray < epsilon_min, 2)/size(errGlobalArray,2);


Err = [err_32, err_64, err_128, err_256, err_512];

Err = fliplr(Err)
cmap = gray(256);
colormap(cmap)
imagesc(Err.')


xticklabels = 0.1:0.1:1.;
xticks = linspace(1, size(Err.', 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)


yticklabels = 2.^(9:-1:5);
yticks = linspace(1, size(Err.', 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
set(gca, 'FontSize', 18)