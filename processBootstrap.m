hasMissing      = cellfun(@(seq)any(all(isnan(seq), 2)), trajectories);
isMissing       = isnan(rawPhi);

rawTraj     = cellfun(@(seq)seq(:, :), rawTrajectories, ...
    'UniformOutput', false);

% [model, modelOrder]     = sortHmmModel(model, 'sortMethod', 'adas');

paths       = cellfun(@(seq)colvec(hmmMap(model, seq')), traj, ...
    'UniformOutput', false);
stackedStates   = colvec(cat(1, paths{:}));
[~, ~, alphas]  = cellfun(@(seq)hmmInferNodes(model, seq'), traj, ...
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

kk              = 20;
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
ax.XLabel.String        = 'States in Path';
ax.YLabel.String        = 'Count';
ax.XTickLabelRotation   = 90;
ax.XLim                 = [0 kk+1];
bb.FaceColor            = [0.4660 0.6740 0.1880];
ax.YGrid                = 'on';
box on;
% tightfig;
end

uncA        = countTransitions(paths(uncensored), model.nstates);
empiricalA  = countTransitions(paths, model.nstates);

if true
subNodes    = unique(cat(2, mciTrans{1:kk}));
% activeMask  = any(triu(uncA, 1) > 5, 2) | any(triu(uncA, 1) > 5, 1);
% subNodes    = subNodes(ismember(subNodes, find(activeMask)));
subGraph    = empiricalA(subNodes, subNodes);
subGraph(subGraph < 5)     = 0;
subGraph(1:length(subGraph)+1:end)  = 0;

% subNodes        = [7  12  14  15  16  17  18  19  20  21  22];
% subNodes        = [3 4 6 7 8 9 10 11 12];
% yPos            = [2 -1.5 -2  2.5 1.5 -1  1.5 1   -1  1   0 ];
% colorIndices    = [1  1   1   1   1   1   2   2   2   2   2 ];

G           = digraph(subGraph);

fig         = figure;
ax          = gca;

colormap(fig, ax.ColorOrder([end-1 2], :));
lWidth      = 30*G.Edges.Weight/max(G.Edges.Weight);
graphPlot   = plot(G, 'edgelabel', G.Edges.Weight, ...
    'linewidth', lWidth, 'nodelabel', arrayfun(@num2str, subNodes, ...
    'uniformoutput', false), 'layout', 'layered');
%     'NodeCData', colorIndices);

ax.FontSize     = 24;
ax.XTick        = [];
ax.YTick        = [];

graphPlot.MarkerSize    = 80;
% graphPlot.NodeColor     = [0.9290 0.6940 0.1250];
graphPlot.ArrowSize     = 25;
% graphPlot.NodeLabel     = {};
% graphPlot.EdgeLabel     = {};
graphPlot.EdgeColor     = ax.ColorOrder(3, :);
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

uncModel    = hmmProfileFitFullyObs(paths(uncensored), ...
    trajectories(uncensored), 'discrete', 'nStates', model.nstates);

emStates        = 1:model.nstates;
% emStates        = [16 17 18 19 20 21 22];
edgeMatrices    = {[16 18; 18 21; 21 22], [17 20; 20 22], ...
    [19 21; 21 22]};
edgeColors      = [0.4660 0.6740 0.1880; ...
    0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250; 0 0.4470 0.7410];
edgeLines       = {'--', '-.', '-.'};

plotIndices     = {{1 6 7 10}, {8 9 12 13}, {11, 14, 15, 16}};

visualizeHmm(allData, traj, rawTraj, dx, age, ...
    model, 'showAnart', true, 'showApoe', true, 'showTrans', false, ...
    'showReading', true, ...
    'showDx', true, 'showDurations', true, 'showEmission', true, ...
    'showEmpiricalTrans', true, 'showAge', true, ...
    'filterMask', true(size(uncensored)), 'emissionParams', ...
    {'emissionPath', emStates, ...
    'edgeMatrices', [], 'edgeColors', edgeColors, ...
    'edgeLines', edgeLines, 'sampleObs', false, 'nBootSamples', 1, ...
    'plotIndices', plotIndices, 'quantileLims', ...
    permute(biomarkerCi, [3 1 2])});
