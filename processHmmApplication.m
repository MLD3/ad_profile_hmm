rocScore        = cellfun(@(seq)seq{1}, auroc);

rocScore        = reshape(rocScore, nRestarts, length(kGrid), ...
    length(trajGrid), nDurations, nClfSeriesOpt);

% if ~exist('pClassifiers', 'var')
%     pClassifiers    = zeros(nHmmExp, nHmmExp, nDurations);
%     
%     for d=1:nDurations
%         mask        = ~censored(:, d);
%         decMat      = cell2mat(cat(2, clfPredictions{:, d}));
%         [~, ~, ~, ~, ~, ~, pClassifiers(:, :, d)]   = rocComp(decMat, ...
%             labels(mask, d));
%     end
% end

rocStats            = quantile(rocScore, [0.025 0.5 0.975], 1);

modalityNames       = {'BASELINE+DEMOG+GENETIC', 'BASELINE'};

dd              = 3;
clfSi           = 2;

figure; hold on;
% plot(kGrid, rocScore(:, :, dd), 'o--', 'linewidth', 2);
errorbar(kGrid, squeeze(rocStats(2, :, 1, dd, clfSi)), ...
    squeeze(rocStats(2, :, 1, dd, clfSi) - rocStats(1, :, 1, dd, clfSi)), ...
    squeeze(rocStats(3, :, 1, dd, clfSi) - rocStats(2, :, 1, dd, clfSi)), ...
    'o--', 'linewidth', 2);
errorbar(kGrid, squeeze(rocStats(2, :, 2, dd, clfSi)), ...
    squeeze(rocStats(2, :, 2, dd, clfSi) - rocStats(1, :, 2, dd, clfSi)), ...
    squeeze(rocStats(3, :, 2, dd, clfSi) - rocStats(2, :, 2, dd, clfSi)), ...
    'o--', 'linewidth', 2);
grid on; box on;
ax              = gca;
ax.XLim         = [min(kGrid)-0.8 max(kGrid)+0.8];
legend('Baseline+Genetic+Demographic', 'Baseline');
ylabel('AUROC');
xlabel('Number of states');
ax.FontSize     = 30;
ax.XTick        = min(kGrid):2:max(kGrid);

if false
modalityIdx     = reshape(colvec(1:nHmmExp), size(kOpt));
row             = 1;
col             = 3;
figure;
imagesc(pClassifiers(modalityIdx(:, row), modalityIdx(:, col), dd));
axis square; colorbar;
ax                  = gca;
ax.XTick            = 1:length(kGrid);
ax.YTick            = 1:length(kGrid);
ax.XTickLabel       = kGrid;
ax.YTickLabel       = kGrid;
ax.YLabel.String    = modalityNames{row};
ax.XLabel.String    = modalityNames{col};

sumLL           = reshape(colvec(squeeze(sum(loglik, 1))), size(kOpt));
sumLLZ          = reshape(colvec(squeeze(sum(loglikZ, 1))), size(kOpt));

figure; hold on;
plot(kGrid, sumLL, 'o--', 'linewidth', 2);
grid on; box on;
ax              = gca;
ax.XLim         = [min(kGrid)-0.8 max(kGrid)+0.8];
legend('Baseline+Genetic+Demographic', 'Baseline+Demographic', ...
    'Baseline+Genetic', 'Baseline');
ylabel('AUROC');
xlabel('Number of states');
ax.FontSize     = 30;
ax.XTick        = min(kGrid):2:max(kGrid);
end

if false

regErrors       = cell2mat(cat(1, testError{:}));
regErrors       = arrayfun(@(seq)seq.all(1), regErrors);
regErrors       = reshape(regErrors, size(kOpt));

figure; hold on;
plot(kGrid, regErrors, 'o--', 'linewidth', 1.5);
grid on; box on;
ax              = gca;
ax.XLim         = [min(kGrid)-0.8 max(kGrid)+0.8];
legend('Baseline+Genetic+Demographic', 'Baseline+Demographic', ...
    'Baseline+Genetic', 'Baseline');
ylabel('RMSE');
xlabel('Number of states');
ax.FontSize     = 30;
ax.XTick        = min(kGrid):2:max(kGrid);

end