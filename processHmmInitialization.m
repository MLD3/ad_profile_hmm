sumLL           = reshape(sum(loglik), size(kOpt));

convtimeLL      = squeeze(sumLL(:, strcmpi(initGrid, 'convtime'), :, :));
scoreLL         = squeeze(sumLL(:, strcmpi(initGrid, 'score'), :, :));

mixtureLL       = sumLL(:, strcmpi(initGrid, 'score-perturbed'), :, :);
randomLL        = sumLL(:, strcmpi(initGrid, 'random'), :, :);

mixtureLL       = permute(mixtureLL, [2 1 3 4]);
randomLL        = permute(randomLL, [2 1 3 4]);

if false
    sz              = size(mixtureLL);
    
    mixtureBoot     = bootstrp(1000, @max, mixtureLL);
    statMix         = quantile(mixtureBoot, [0.25 0.5 0.75], 1);
    statMix         = reshape(statMix, [3 sz(2:end)]);
    
    randomBoot      = bootstrp(1000, @max, randomLL);
    statRandom      = quantile(randomBoot, [0.25 0.5 0.75], 1);
    statRandom      = reshape(statRandom, [3 sz(2:end)]);
else
%     sampleIdx       = randsample(nRandomRestarts, nRandomRestarts);
%     mixtureLL       = mixtureLL(sampleIdx, :, :, :);
%     randomLL        = randomLL(sampleIdx, :, :, :);
    
    statMix         = quantile(mixtureLL, [0.0 0.5 0.99], 1);
    statRandom      = quantile(randomLL, [0.0 0.5 0.99], 1);
end

muMix           = mean(mixtureLL, 1);
seMix           = std(mixtureLL, 0, 1)/sqrt(nRandomRestarts);

muRandom        = mean(randomLL, 1);
seRandom        = std(randomLL, 0, 1)/sqrt(nRandomRestarts);

tIdx        = 4;

figure('Name', 'Log-likelihood'); hold on;
errorbar(kGrid, statMix(2, :, tIdx, 1), ...
    statMix(2, :, tIdx, 1) - statMix(1, :, tIdx, 1), ...
    statMix(3, :, tIdx, 1) - statMix(2, :, tIdx, 1), 'o--', ...
    'linewidth', 1.5);
errorbar(kGrid, statMix(2, :, tIdx, 2), ...
    statMix(2, :, tIdx, 2) - statMix(1, :, tIdx, 2), ...
    statMix(3, :, tIdx, 2) - statMix(2, :, tIdx, 2), 's--', ...
    'linewidth', 1.5);
errorbar(kGrid, statRandom(2, :, tIdx, 1), ...
    statRandom(2, :, tIdx, 1) - statRandom(1, :, tIdx, 1), ...
    statRandom(3, :, tIdx, 1) - statRandom(2, :, tIdx, 1), 'o--', ...
    'linewidth', 1.5);
errorbar(kGrid, statRandom(2, :, tIdx, 2), ...
    statRandom(2, :, tIdx, 2) - statRandom(1, :, tIdx, 2), ...
    statRandom(3, :, tIdx, 2) - statRandom(2, :, tIdx, 2), 's--', ...
    'linewidth', 1.5);
plot(kGrid, scoreLL(:, tIdx, 1), 'o-.', 'linewidth', 1.5);
plot(kGrid, convtimeLL(:, tIdx, 1), 's-.', 'linewidth', 1.5);
legend('ADAS-Random: General Model', 'ADAS-Random: LTR Model', ...
    'Random: General Model', 'Random: LTR Model', ...
    'ADAS-Cog', 'Time-to-conversion');
grid on; box on;
ax          = gca;
ax.XTick    = kGrid;
ax.XLim     = [min(kGrid)-0.8 max(kGrid)+0.8];

figure('Name', 'Log-likelihood'); hold on;
plot(kGrid, statMix(3, :, tIdx, 1), 'o--', 'linewidth', 1.5);
plot(kGrid, statMix(3, :, tIdx, 2), 'o--', 'linewidth', 1.5);
plot(kGrid, statRandom(3, :, tIdx, 1), 'd--', 'linewidth', 1.5);
plot(kGrid, statRandom(3, :, tIdx, 2), 'd--', 'linewidth', 1.5);
plot(kGrid, scoreLL(:, tIdx, 1), '^-.', 'linewidth', 1.5);
plot(kGrid, convtimeLL(:, tIdx, 1), 's-.', 'linewidth', 1.5);
legend('ADAS-Random: General Model', 'ADAS-Random: LTR Model', ...
    'Random: General Model', 'Random: LTR Model', ...
    'ADAS-Cog', 'Time-to-conversion');
grid on; box on;
ax          = gca;
ax.XTick    = kGrid;
ax.XLim     = [min(kGrid)-0.8 max(kGrid)+0.8];


% grpRandom       = 1:length(kGrid);
% grpMixture      = max(grpRandom)+(1:length(kGrid));
% posRandom       = kGrid-0.2;
% posMixture      = kGrid+0.2;
% colorRandom     = ones(size(kGrid))*3;
% colorMixture    = ones(size(kGrid))*4;
% colorIdx        = interweave(colorRandom, colorMixture);
% 
% figure; hold on;
% ax              = gca;
% boxplot([mixtureLL randomLL], [grpMixture grpRandom], 'positions', ...
%     [posMixture posRandom]);
% 
% h       = findobj(ax, 'Tag', 'Box');
% for j=1:length(h)
%     patch(get(h(j),'XData'), get(h(j),'YData'), ...
%         ax.ColorOrder(colorIdx(j), :), 'FaceAlpha', 0.5);
% end
% plot(kGrid, convtimeLL, 'o--', 'linewidth', 2);
% plot(kGrid, scoreLL, 's--', 'linewidth', 2);
% 
% boxes   = findobj(ax, 'Type', 'Patch');
% lines   = findobj(ax, 'Type', 'Line', 'Tag', '');
% 
% ax.XTick        = kGrid;
% ax.XTickLabel   = kGrid;
% legend([boxes(1) boxes(2) lines'], 'Random', 'Mixture', ...
%     'ADAS-Cog', 'Time to conversion');
% grid on; box on;