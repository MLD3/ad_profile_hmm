if ~exist('scoreBaseline', 'var')
    prop        = load('profile-hmm-application-random.mat');
    scBase      = load('profile-hmm-score-application-no-clinical');
    dxBase      = load('profile-hmm-convtime-application-no-clinical');
end

kOpt            = prop.kOpt;
tIdx            = 1;
cIdx            = 2;

prop.rocScore   = cellfun(@(seq)seq{1}, prop.auroc);
prop.rocScore   = reshape(prop.rocScore, prop.nRestarts, ...
    length(prop.kGrid), length(prop.trajGrid), ...
    prop.nDurations, prop.nClfSeriesOpt);
prop.rocScore   = squeeze(prop.rocScore(:, :, tIdx, :, cIdx));
prop.rocStat    = quantile(prop.rocScore, [0.025 0.5 0.975], 1);

prop.prob       = reshape(prop.clfPredictions, prop.nRestarts, ...
    length(prop.kGrid), length(prop.trajGrid), ...
    prop.nDurations, prop.nClfSeriesOpt);
prop.prob       = squeeze(prop.prob(:, tIdx, :));

dxBase.rocScore         = cellfun(@(seq)seq{1}, dxBase.auroc);
dxBase.rocScore         = reshape(dxBase.rocScore, ...
    size(dxBase.kOpt, 1), size(dxBase.kOpt, 2), []);
dxBase.rocScore         = squeeze(dxBase.rocScore(:, tIdx, :));
dxBase.prob             = reshape(dxBase.clfPredictions, ...
    size(dxBase.kOpt, 1), size(dxBase.kOpt, 2), []);
dxBase.prob             = squeeze(dxBase.prob(:, tIdx, :));

scBase.rocScore      = cellfun(@(seq)seq{1}, scBase.auroc);
scBase.rocScore      = reshape(scBase.rocScore, ...
    size(scBase.kOpt, 1), size(scBase.kOpt, 2), []);
scBase.rocScore      = squeeze(scBase.rocScore(:, tIdx, :));
scBase.prob          = reshape(scBase.clfPredictions, ...
    size(scBase.kOpt, 1), size(scBase.kOpt, 2), []);
scBase.prob          = squeeze(scBase.prob(:, tIdx, :));

if false
prop.regError           = reshape(cell2mat(cat(1, ...
    prop.testError{:})), size(kOpt));
prop.rmse               = arrayfun(@(seq)seq.all(1), prop.regError);
prop.ci                 = arrayfun(@(seq)seq.all(2), prop.regError);

scBase.regError      = reshape(cell2mat(cat(1, ...
    scBase.testError{:})), size(kOpt));
scBase.rmse          = arrayfun(@(seq)seq.all(1), ...
    scBase.regError);
scBase.ci            = arrayfun(@(seq)seq.all(2), ...
    scBase.regError);

dxBase.regError       = reshape(cell2mat(cat(1, ...
    dxBase.testError{:})), size(kOpt));
dxBase.rmse           = arrayfun(@(seq)seq.all(1), ...
    dxBase.regError);
dxBase.ci             = arrayfun(@(seq)seq.all(2), ...
    dxBase.regError);

rmseMat         = permute(cat(3, scoreBaseline.rmse, convBaseline.rmse, ...
    proposed.rmse), [1 3 2]); 
ciMat           = permute(cat(3, scoreBaseline.ci, convBaseline.ci, ...
    proposed.ci), [1 3 2]);

figure; hold on;
plot(kGrid, rmseMat(:, :, tIdx), 'o--', 'linewidth', 2);
grid on; box on;
ax              = gca;
ax.XLim         = [min(kGrid)-0.8 max(kGrid)+0.8];
legend('ADAS-Cog Baseline', 'Conv. Time Baseline', ...
    'HMM Model');
ylabel('AUROC');
xlabel('Number of states');
ax.FontSize     = 30;
ax.XTick        = min(kGrid):2:max(kGrid);
end

prop.sumLL          = sum(prop.loglik);
prop.sumLL          = reshape(prop.sumLL, prop.nRestarts, ...
    length(prop.kGrid), length(prop.trajGrid));
prop.llStat         = quantile(prop.sumLL, [0.025 0.5 0.975], 1);

scBase.sumLL        = reshape(colvec(squeeze(sum(...
    scBase.loglik, 1))), size(scBase.kOpt));
dxBase.sumLL          = reshape(colvec(squeeze(sum(...
    dxBase.loglik, 1))), size(dxBase.kOpt));

kGrid           = prop.kGrid;
dd              = 6;

figure; hold on;
plot(scBase.kGrid, scBase.sumLL(:, tIdx), 'o--', 'linewidth', 2);
plot(dxBase.kGrid, dxBase.sumLL(:, tIdx), 'o--', 'linewidth', 2);
errorbar(prop.kGrid, prop.llStat(2, :, tIdx), ...
    prop.llStat(2, :, tIdx)-prop.llStat(1, :, tIdx), ...
    prop.llStat(3, :, tIdx)-prop.llStat(2, :, tIdx), ...
    'd--', 'linewidth', 2);
grid on; box on;
ax              = gca;
ax.XLim         = [min(prop.kGrid)-0.8 max(prop.kGrid)+0.8];
legend('ADAS-Cog Baseline', 'Conv. Time Baseline', ...
    'HMM Model');
ylabel('Log-likelihood');
xlabel('Number of states');
ax.FontSize     = 30;
ax.XTick        = min(prop.kGrid):2:max(prop.kGrid);

figure; hold on;
plot(scBase.kGrid, scBase.rocScore(:, dd), 'o--', 'linewidth', 2);
plot(dxBase.kGrid, dxBase.rocScore(:, dd), 'o--', 'linewidth', 2);
errorbar(prop.kGrid, prop.rocStat(2, :, dd), ...
    prop.rocStat(2, :, dd)-prop.rocStat(1, :, dd), ...
    prop.rocStat(3, :, dd)-prop.rocStat(2, :, dd), ...
    'd--', 'linewidth', 2);
grid on; box on;
ax              = gca;
ax.XLim         = [min(prop.kGrid)-0.8 max(prop.kGrid)+0.8];
legend('ADAS-Cog Baseline', 'Conv. Time Baseline', ...
    'HMM Model');
ylabel('AUROC');
xlabel('Number of states');
ax.FontSize     = 30;
ax.XTick        = min(prop.kGrid):2:max(prop.kGrid);



decMatScore     = cell2mat(cat(2, scBase.prob{:, dd}));
decMatConv      = cell2mat(cat(2, dxBase.prob{:, dd}));
decMatProposed  = cell2mat(cat(2, prop.prob{:, dd}));
lab             = prop.labels(~prop.censored(:, dd), dd);

[area, ~, ~, ~, ~, ~, pval]     = rocComp([decMatConv decMatProposed], ...
    lab);