hasMissing      = cellfun(@(seq)any(all(isnan(seq), 2)), trajectoriesNan);

pvals           = zeros(3, 3, length(kGrid));

for k=1:length(kGrid)
    pvals(:, :, k)      = pairwiseHypothesisTest([missingLoglikZ(:, k) ...
        paddedLoglikZ(:, k) nanLoglikZ(:, k)], 'wilcoxon');
end

[aicNan, bicNan, aicRatio]  = getHmmAic(nanModels, nanLoglik);
[aicPad, bicPad]            = getHmmAic(paddedModels, paddedLoglik);
[aicMiss, bicMiss]          = getHmmAic(missingModels, missingLoglik);

idx                 = 3;
paddedModel         = paddedModels{idx};
nanModel            = nanModels{idx};
missingModel        = missingModels{idx};

target              = 'nan';

if strcmpi(target, 'nan')
    model           = nanModel;
    aic             = aicNan;
    bic             = bicNan;
elseif strcmpi(target, 'padded')
    model           = paddedModel;
    aic             = aicPad;
    bic             = bicPad;
elseif strcmpi(target, 'missing')
    model           = missingModel;
    aic             = aicMiss;
    bic             = bicMiss;
end

% paths       = cellfun(@(seq)hmmMap(model, seq'), trajectoriesNan, ...
%     'UniformOutput', false);
% model       = sortHmmModel(model, paths, dx);

paths       = cellfun(@(seq)hmmMap(model, seq'), trajectoriesNan, ...
    'UniformOutput', false);
stackedPaths    = colvec(cat(2, paths{:}));

%% find the unique paths amongst uncensored data
letters         = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', ...
    'k', 'l', 'm', 'n'};
stackedLetters  = colvec(letters(stackedPaths));
letterPaths     = mat2cell(stackedLetters, counts);
letterPaths     = cellfun(@(seq)cat(2, seq{:}), letterPaths, ...
    'UniformOutput', false);

uncensored      = cellfun(@(seq)ismember(2, seq) && ismember(3, seq), dx);

uncLetterPaths  = letterPaths(uncensored);

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
        transPos    = find(diff(paths{ii}));
        
        for t=1:length(transPos)
            pos         = transPos(t);
            stateSeq    = cat(1, stateSeq, ...
                [paths{ii}(pos) paths{ii}(pos+1)]);
            distance    = adPos - pos - 1;
        end
    end
    if length(unique(paths{ii})) > 1
        hasTrans(ii)    = true;
    end
end
uncA        = countTransitions(paths(uncensored), model.nstates);

starter     = find(pos == 1);
seenStates  = unique(cat(2, paths{starter}));

% edgeMatrix                  = false(size(empiricalA));
% edgeMatrix(7:end, 7:end)    = empiricalA(7:end, 7:end) > 5;

src         = 1:model.nstates-1;
dest        = 2:model.nstates;
edgeMatrix  = full(adjacency(digraph(src, dest)));
% emStates        = [7 8 11 13 14];
emStates        = 1:model.nstates;
visualizeHmm(allData, trajectoriesNan, rawTrajectoriesPad, dxPad, ...
    model, 'showAnart', true, 'showApoe', true, 'showTrans', false, ...
    'showDx', true, 'showDurations', false, 'showEmission', true, ...
    'emissionPath', emStates, 'showEmpiricalTrans', true, ...
    'edgeMatrix', edgeMatrix(emStates, emStates));

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