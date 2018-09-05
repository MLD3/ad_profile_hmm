hasMissing      = cellfun(@(seq)any(all(isnan(seq), 2)), trajectories);
isMissing       = isnan(rawPhi);

[aic, bic, aicRatio]    = getHmmAic(colvec(hmmModels), loglik);

aic         = reshape(aic, size(kOpt));

idx         = 7;

model       = models{idx};

paths       = cellfun(@(seq)colvec(hmmMap(model, seq')), trajectories, ...
    'UniformOutput', false);
stackedStates   = colvec(cat(1, paths{:}));
[~, ~, alphas]  = cellfun(@(seq)hmmInferNodes(model, seq'), trajectories, ...
    'UniformOutput', false);
stackedAlphas   = cat(2, alphas{:})';

index           = sub2ind(size(stackedAlphas), ...
    colvec(1:length(stackedStates)), stackedStates);

%% find the unique paths amongst uncensored data
letters         = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', ...
    'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', ...
    'x', 'y', 'z'};
stackedLetters  = colvec(letters(stackedStates));
letterPaths     = mat2cell(stackedLetters, counts);
letterPaths     = cellfun(@(seq)cat(2, seq{:}), letterPaths, ...
    'UniformOutput', false);
letterTrans     = cellfun(@unique, letterPaths, 'uniformoutput', false);

uncensored      = cellfun(@(label,yupper)any(label == 2) && ...
    ~any(yupper == 999), dx, yu);

yuMci           = arrayfun(@(ii)colvec(allData.data.CMCI_U(...
    allData.trajIdx{ii})), 1:length(allData.trajIdx), ...
    'uniformoutput', false);
% uncnormal       = cellfun(@(label, yupper)any(label == 1) && ...
%     ~any(yupper == 999), dx, yuMci);
uncnormal       = cellfun(@(label)any(label == 1) && any(label == 2), dx);

mciLetterTrans  = letterTrans(uncensored);
nlLetterTrans   = letterTrans(uncnormal);

[allLetterMciTrans, ~, allMciTransGroups]   = unique(mciLetterTrans);
mciTrans            = cell(size(allLetterMciTrans));
mciTransCounts      = histc(allMciTransGroups, 1:length(allLetterMciTrans));

[allLetterNlTrans, ~, allNlTransGroups]     = unique(nlLetterTrans);
nlTrans             = cell(size(allLetterNlTrans));
nlTransCounts       = histc(allNlTransGroups, 1:length(allLetterNlTrans));

for ii=1:length(allLetterMciTrans)
    letterPath      = allLetterMciTrans{ii};
    zPath           = zeros(size(letterPath));
    for l=1:length(letterPath)
        zPath(l)    = find(strcmpi(letters, letterPath(l)));
    end
    mciTrans{ii}    = zPath;
end

for ii=1:length(allLetterNlTrans)
    letterPath      = allLetterNlTrans{ii};
    zPath           = zeros(size(letterPath));
    for l=1:length(letterPath)
        zPath(l)    = find(strcmpi(letters, letterPath(l)));
    end
    nlTrans{ii}     = zPath;
end

[~, order]          = sort(mciTransCounts, 'descend');
mciTransCounts      = mciTransCounts(order);
allLetterMciTrans   = allLetterMciTrans(order);
mciTrans            = mciTrans(order);
allMciTransGroups   = order(allMciTransGroups);

[~, order]          = sort(nlTransCounts, 'descend');
nlTransCounts       = nlTransCounts(order);
allLetterNlTrans    = allLetterNlTrans(order);
nlTrans             = nlTrans(order);
allNlTransGroups    = order(allNlTransGroups);

kk              = 10;
topKCounts      = mciTransCounts(1:kk);
topKLabels      = cell(kk, 1);

for ii=1:kk
    nSymbols        = length(mciTrans{ii});
    label           = sprintf('%d', mciTrans{ii}(1));
    for s=2:nSymbols
        label       = cat(2, label, sprintf(',%d', mciTrans{ii}(s)));
    end
    topKLabels{ii}  = label;
end

if true
figure('units', 'normalize', 'outerposition', [0 0 1 1]);
bb      = bar(topKCounts);
ax              = gca;
ax.XTick        = 1:kk;
ax.XTickLabel   = topKLabels;
ax.FontSize     = 32;
ax.XTickLabelRotation   = 90;
ax.XLim                 = [0 kk+1];
bb.FaceColor            = [0.4660 0.6740 0.1880];
grid on; box on;
tightfig;
end

stateSeq        = [];
distance        = [];

stateIdx    = 9;
pos         = nan(size(paths));
hasTrans    = false(size(paths));

for ii=1:length(paths)
    if ismember(stateIdx, paths{ii})
        pos(ii)     = find(paths{ii} == stateIdx, 1, 'first');
    end
    if ismember(2, dx{ii}) && ismember(3, dx{ii}) 
        adPos       = find(dx{ii} == 3, 1, 'first');
        
        if all(dx{ii}(adPos+1:end) == 3)
            transPos    = find(diff(paths{ii}));
            
            for t=1:length(transPos)
                pos         = transPos(t);
                stateSeq    = cat(1, stateSeq, ...
                    [paths{ii}(pos) paths{ii}(pos+1)]);
                distance    = cat(1, distance, adPos - pos - 1);
            end
            
            if length(transPos) == 1
                stateSeq    = cat(1, stateSeq, ...
                    [paths{ii}(adPos-1) paths{ii}(adPos)]);
                distance    = cat(1, distance, 0);
            end
        end
    end
    if length(unique(paths{ii})) > 1
        hasTrans(ii)    = true;
    end
end
[unqOnsetCodes, ~, onsetProfile]    = unique(stateSeq, 'rows');
onsetCounts         = histc(onsetProfile, 1:size(unqOnsetCodes, 1));
[uDist, ~, distProfile]     = unique(distance);

heatmap     = accumarray([onsetProfile distProfile], 1);

yLabels     = cell(size(unqOnsetCodes));
for ii=1:length(yLabels)
    yLabels{ii}     = sprintf('%d->%d', unqOnsetCodes(ii, 1), ...
        unqOnsetCodes(ii, 2));
end

if false
figure;
imagesc(heatmap(diff(unqOnsetCodes, 1, 2) ~= 0, :)); colorbar;
ax          = gca;
ax.XTick    = 1:length(uDist);
ax.XTickLabel   = uDist;
size(unqOnsetCodes, 1);
ax.YTick        = 1:sum(diff(unqOnsetCodes, 1, 2) ~= 0);
ax.YTickLabel   = yLabels(diff(unqOnsetCodes, 1, 2) ~= 0);
end

uncA        = countTransitions(paths(uncensored), model.nstates);
empiricalA  = countTransitions(paths, model.nstates);

uncModel    = hmmProfileFitFullyObs(paths(uncensored), ...
    trajectories(uncensored), 'discrete', 'nStates', model.nstates);

edgeMatrices    = {[7 9; 9 11; 11 12], [8 10; 10 12]};
edgeColors      = [0.4660 0.6740 0.1880; 0 0.4470 0.7410; ...
    0.9290 0.6940 0.1250];
edgeLines       = {'-.', '--', '-.'};

emStates    = 1:model.nstates;

% visualizeHmm(allData, trajectories, rawTrajectories, dx, age, ...
%     model, 'showAnart', true, 'showApoe', false, 'showTrans', false, ...
%     'showDx', false, 'showDurations', false, 'showEmission', false, ...
%     'showEmpiricalTrans', false, 'emissionParams', ...
%     {'emissionPath', emStates, ...
%     'edgeMatrices', edgeMatrices, 'edgeColors', edgeColors, ...
%     'edgeLines', edgeLines});

visualizeHmm(allData, trajectories, rawTrajectories, dx, age, ...
    model, 'showAnart', true, 'showApoe', false, 'showTrans', false, ...
    'showReading', true, ...
    'showDx', true, 'showDurations', false, 'showEmission', true, ...
    'showEmpiricalTrans', true, 'showAge', true, ...
    'filterMask', true(size(uncensored)), 'emissionParams', ...
    {'emissionPath', emStates, ...
    'edgeMatrices', [], 'edgeColors', edgeColors, ...
    'edgeLines', edgeLines, 'sampleObs', false, 'nBootSamples', 1});

%% Visualize the AIC vs K
if false
figure; hold on;
plot(kGrid, aic, 'o--', 'linewidth', 2);
plot(kGrid, bic, 's--', 'linewidth', 2);
legend('Akaike Information Criterion (AIC)', 'BIC');
xlabel('Number of states (K)')
grid on; box on;
ax              = gca;
ax.FontSize     = 32;
ax.XTick        = 0:2:max(kGrid);
ax.XLim         = [min(kGrid)-0.8 max(kGrid)+0.8];
end

%% Train a mixture-model on cross-sectional data to compare 

if false
K               = model.nstates;
initParams      = struct('mixWeight', normalize(ones(1, K)), ...
    'T', model.emission.T);
mixModel        = mixDiscreteFit(discretePhi, K, 'initParams', ...
    initParams, 'verbose', true);

hmmPath         = cellfun(@(seq)hmmMap(model, seq'), trajectoriesNan, ...
    'UniformOutput', false);
hmmStates       = cat(2, hmmPath{:});
[~, mixStates]  = max(mixDiscreteInferLatent(mixModel, discretePhi), [], 2);    
end