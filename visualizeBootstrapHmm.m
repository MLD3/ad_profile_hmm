function visualizeBootstrapHmm(allData, baseTrajectories, ...
    baseRawTrajectories, baseDx, baseAge, baseYl, baseYu, ...
    bootModels, bootSamples, bootStates, varargin)

[   showAnart                                           , ...
    showReading                                         , ...
    showAge                                             , ...
    showApoe                                            , ...
    showEmpiricalTrans                                  , ...
    showDx                                              , ...
    showDurations                                       , ...
    showPathDist                                        , ...
    topKCount                                           , ...
    showEmission                                        , ...
    quantileValues                                      , ...
    emissionParams          ] = process_options(varargin, ...
    'showAnart'             , true                      , ...
    'showReading'           , true                      , ...
    'showAge'               , true                      , ...
    'showApoe'              , true                      , ...
    'showEmpiricalTrans'    , true                      , ...
    'showDx'                , true                      , ...
    'showDurations'         , false                     , ...
    'showPathDist'          , true                      , ...
    'topKCount'             , 15                        , ...
    'showEmission'          , true                      , ...
    'quantileValues'        , [0.25 0.5 0.75]           , ...
    'emissionParams'        , {}                        );

nBoot           = size(bootSamples, 2);

instIndex       = mat2cell(colvec(1:size(allData.data, 1)), ...
    allData.counts);

allAnart        = cell(nBoot, 1);
allAge          = cell(nBoot, 1);

allDxFreq       = cell(nBoot, 1);
allGeneFreq     = cell(nBoot, 1);
allReadFreq     = cell(nBoot, 1);

allFirstState   = cell(nBoot, 1);
allPaths        = cell(nBoot, 1);

allLetterStateSets  = cell(nBoot, 1);
allEmpiricalTrans   = cell(nBoot, 1);

nStates             = bootModels{1}.nstates;
stateIdx            = colvec(1:nStates);

letters         = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', ...
    'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', ...
    'x', 'y', 'z'};

for bb=1:nBoot
    sampleIdx           = bootSamples(:, bb);
    
    bootCounts          = allData.counts(sampleIdx);
    bootInstances       = cell2mat(instIndex(sampleIdx));
    
    stackedStates       = bootStates{bb};
    model               = bootModels{bb};
    
    dx                  = baseDx(sampleIdx);
    age                 = baseAge(sampleIdx);
    yu                  = baseYu(sampleIdx);
    
    uncensored          = cellfun(@(label, yupper)any(label == 2) && ...
        ~any(yupper == 999), dx, yu);
    
    stackedLetters      = colvec(letters(stackedStates));
    letterPaths         = mat2cell(stackedLetters, bootCounts);
    letterPaths         = cellfun(@(seq)cat(2, seq{:}), letterPaths, ...
        'uniformoutput', false);
    letterStateSet      = cellfun(@unique, letterPaths, 'uniformoutput', ...
        false);
    
    allLetterStateSets{bb}      = letterStateSet(uncensored);
        
    firstIdx            = cumsum([1; bootCounts(1:end-1)]);
    allFirstState{bb}   = stackedStates(firstIdx);
    
    allAnart{bb}        = allData.data.ANARTERR(bootInstances);
    allAnart{bb}        = allAnart{bb}(firstIdx);
    allAge{bb}          = cell2mat(age);
    
    stackedDx       = cell2mat(dx);
    knownDx         = ~isnan(stackedDx);
    allDxFreq{bb}   = accumarray([stackedStates(knownDx) ...
        stackedDx(knownDx)], 1);
    
    stackedGenes    = getGeneProfile(allData.data(bootInstances, :));
    knownGene       = ~isnan(stackedGenes);
    allGeneFreq{bb} = accumarray([stackedStates(knownGene) ...
        stackedGenes(knownGene)], 1);
    
    reading             = discretize(allData.data.PTEDUCAT(bootInstances), ...
        [0 15 18 inf]);
    knownReading        = ~isnan(reading);
    allReadFreq{bb}     = accumarray([stackedStates(knownReading) ...
        reading(knownReading)], 1);
    
    allPaths{bb}    = mat2cell(bootStates{bb}, bootCounts);
    
    uncPaths                = allPaths{bb}(uncensored);
    allEmpiricalTrans{bb}   = countTransitions(uncPaths, nStates);
end

if showAnart
    figure; hold on;
    boxplot(cell2mat(allAnart), cell2mat(allFirstState));
    ax                  = gca;
    ax.XLabel.String    = 'HMM State';
    ax.YLabel.String    = '# ANART Errors';
    ax.FontSize         = 30;
    ax.XLim             = [0.2 nStates+0.8];
    ax.YGrid            = 'on';
end

if showAge
    figure; hold on;
    boxplot(cell2mat(allAge), cell2mat(colvec(bootStates)));
    ax                  = gca;
    ax.XLabel.String    = 'HMM State';
    ax.YLabel.String    = 'Age';
    ax.FontSize         = 30;
    ax.YGrid            = 'on';
end

if showReading
    normFreq        = cellfun(@(seq)normalize(seq, 2), allReadFreq, ...
        'uniformoutput', false);
    freqMat         = cat(3, normFreq{:});
    freqStats       = quantile(freqMat, quantileValues, 3);
    
    figure; hold on;
    [hBars, hErrors]    = barwitherr(cat(3, freqStats(:, :, 2) - ...
        freqStats(:, :, 1), freqStats(:, :, 3) - freqStats(:, :, 2)), ...
        freqStats(:, :, 2));
    ax                  = gca;
    ax.XLabel.String    = 'HMM State';
    ax.YLabel.String    = 'Probability';
    legend('Low', 'Intermediate', 'High');
    ax.FontSize         = 30;
    ax.XLim             = [0.2 nStates+0.8];
    ax.YGrid            = 'on';
    hBars(1).FaceAlpha  = 0.6;
end

if showApoe
    normFreq        = cellfun(@(seq)normalize(seq, 2), allGeneFreq, ...
        'uniformoutput', false);
    freqMat         = cat(3, normFreq{:});
    freqStats       = quantile(freqMat, quantileValues, 3);
    
    figure; hold on;
    [hBars, hErrors]    = barwitherr(cat(3, freqStats(:, :, 2) - ...
        freqStats(:, :, 1), freqStats(:, :, 3) - freqStats(:, :, 2)), ...
        freqStats(:, :, 2));
    ax                  = gca;
    ax.XLabel.String    = 'HMM State';
    ax.YLabel.String    = 'Probability';
    legend('ApoE 2,2', 'ApoE 3,3', 'ApoE 4,4', 'ApoE 2,3', ...
        'ApoE 3,4', 'ApoE 2,4');
    ax.FontSize         = 30;
    ax.XLim             = [0.2 nStates+0.8];
    ax.YGrid            = 'on';
    hBars(1).FaceAlpha  = 0.6;
end

if showDx
    normFreq        = cellfun(@(seq)normalize(seq, 2), allDxFreq, ...
        'uniformoutput', false);
    freqMat         = cat(3, normFreq{:});
    freqStats       = quantile(freqMat, quantileValues, 3);
    
    figure; hold on;
    [hBars, hErrors]    = barwitherr(cat(3, freqStats(:, :, 2) - ...
        freqStats(:, :, 1), freqStats(:, :, 3) - freqStats(:, :, 2)), ...
        freqStats(:, :, 2));
    ax                  = gca;
    ax.XLabel.String    = 'HMM State';
    ax.YLabel.String    = 'Probability';
    legend('NL', 'MCI', 'AD');
    ax.FontSize         = 30;
    ax.XLim             = [0.2 nStates+0.8];
    ax.YGrid            = 'on';
    hBars(1).FaceColor  = ax.ColorOrder(end-2, :);
    hBars(2).FaceColor  = ax.ColorOrder(1, :);
    hBars(3).FaceColor  = ax.ColorOrder(2, :);
end

if showEmpiricalTrans
    medianTransCount    = quantile(cat(3, allEmpiricalTrans{:}), 0.5, 3);
    
    figure;
    title('Empirical Transition Matrix');    
    imagesc(medianTransCount);
    ax          = gca;
    axis square;
    ax.XLabel.String    = 'Destination Stage';
    ax.YLabel.String    = 'Source Stage';
    ax.FontSize         = 32;
    ax.XTick            = min(stateIdx):2:max(stateIdx);
    ax.YTick            = min(stateIdx):2:max(stateIdx);
    colorbar;
end

if showDurations
    maxDuration         = 20;
    
    durationProb        = zeros(nStates, maxDuration, nBoot);
    stateNames          = cell(nStates, 1);
    
    for bb=1:nBoot
        model       = bootModels{bb};
        for kk=1:nStates
            subMatrix       = stateIdx(kk);
        
            probs           = model.A(subMatrix, subMatrix)*...
                ones(maxDuration, 1);
        
            durationProb(kk, :, bb)     = (1-probs).*...
                cumprod([1; probs(1:end-1)]);
        
            stateNames{kk}  = sprintf('K=%d', kk);
        end
    end
    
    durationStats       = quantile(durationProb, quantileValues, 3);
    
    figure; hold on;
    [hBars, hErrors]    = barwitherr(cat(3, durationStats(:, :, 2) - ...
        durationStats(:, :, 1), durationStats(:, :, 3) - ...
        durationStats(:, :, 2)), durationStats(:, :, 2));
    ax                  = gca;
    ax.XLabel.String    = 'HMM State';
    ax.YLabel.String    = 'Probability of Staying';
    legend(arrayfun(@num2str, (1:maxDuration)*0.5, 'uniformoutput', false));
    ax.FontSize         = 30;
    ax.XLim             = [0.2 nStates+0.8];
    ax.YGrid            = 'on';
    hBars(1).FaceAlpha  = 0.6;
end

if showPathDist
    [uniqueLetterSets, ~, letterSetGroup]   = unique(...
        cat(1, allLetterStateSets{:}));
    
    pathCounts      = cellfun(@length, allLetterStateSets);
    bootSetIdx      = arrayfun(@(ii)ones(pathCounts(ii), 1)*ii, 1:nBoot, ...
        'uniformoutput', false);
    
    setFreqMat      = accumarray([letterSetGroup cat(1, bootSetIdx{:})], 1);
    
    uniqueSets      = cell(size(uniqueLetterSets));
    setLabels       = cell(size(uniqueLetterSets));
    
    for ii=1:length(uniqueLetterSets)
        letterPath      = uniqueLetterSets{ii};
        uniqueSets{ii}  = arrayfun(@(ll)find(strcmpi(letters, ...
            letterPath(ll))), 1:length(letterPath));
        
        label           = sprintf('%d', uniqueSets{ii}(1));
        for ll=2:length(letterPath)
            label       = cat(2, label, sprintf(',%d', uniqueSets{ii}(ll)));
        end
        setLabels{ii}   = label;
    end
    
    setFreqStats    = quantile(setFreqMat, quantileValues, 2);
    
    [~, order]      = sort(setFreqStats(:, 2), 'descend');
    topKIdx         = order(1:topKCount);
    
    figure('units', 'normalize', 'outerposition', [0 0 1 1]); 
    [hBars, hErrors]    = barwitherr(cat(3, ...
        setFreqStats(topKIdx, 2) - setFreqStats(topKIdx, 1), ...
        setFreqStats(topKIdx, 3) - setFreqStats(topKIdx, 2)), ...
        setFreqStats(topKIdx, 2));
    ax                      = gca;
    ax.XTick                = 1:topKCount;
    ax.XTickLabel           = setLabels(topKIdx);
    ax.FontSize             = 32;
    ax.XLabel.String        = 'States in Path';
    ax.YLabel.String        = 'Count';
    ax.XTickLabelRotation   = 90;
    ax.XLim                 = [0 topKCount+1];
    hBars.FaceColor         = [0.4660 0.6740 0.1880];
    ax.YGrid                = 'on';
    box on;

    subNodes        = unique(cat(2, uniqueSets{:}));
    
    medianA         = quantile(cat(3, allEmpiricalTrans{:}), 0.5, 3);
    
    subGraph                            = medianA(subNodes, subNodes);
    subGraph(subGraph < 10)             = 0;
    subGraph(1:length(subGraph)+1:end)  = 0;
    
    G           = digraph(subGraph);
    
    fig         = figure;
    ax          = gca;
    
    colormap(fig, ax.ColorOrder([end-1 2], :));
    lWidth      = 30*G.Edges.Weight/max(G.Edges.Weight);
    graphPlot   = plot(G, 'edgelabel', G.Edges.Weight, ...
        'linewidth', lWidth, 'nodelabel', arrayfun(@num2str, subNodes, ...
        'uniformoutput', false), 'layout', 'layered');
    %     'NodeCData', colorIndices);
    
    ax.FontSize     = 32;
    ax.XTick        = [];
    ax.YTick        = [];
    
    graphPlot.MarkerSize    = 80;
    % graphPlot.NodeColor     = [0.9290 0.6940 0.1250];
    graphPlot.ArrowSize     = 25;
    % graphPlot.NodeLabel     = {};
    % graphPlot.EdgeLabel     = {};
    graphPlot.EdgeColor     = ax.ColorOrder(3, :);
end

if showEmission
    [   emissionPath                                            , ...
        quantileValues                                          , ...
        nSamples                                                , ...
        edgeMatrices                                            , ...
        edgeColors                                              , ...
        plotIndices                                             , ...
        edgeLines           ] = process_options(emissionParams  , ...
        'emissionPath'      , 1:nStates                         , ...
        'quantileValues'    , [0.25 0.5 0.75]                   , ...
        'nSamples'          , 1000                              , ...
        'edgeMatrices'      , []                                , ...
        'edgeColors'        , []                                , ...
        'plotIndices'       , []                                , ...
        'edgeLines'         , {}                                );
    
    d               = size(baseRawTrajectories{1}, 2);
    
    allSamples      = zeros(nStates, d, nBoot*nSamples);
    
    for bb=1:nBoot
        rawPhi      = cell2mat(baseRawTrajectories(bootSamples(:, bb)));
        discPhi     = cell2mat(baseTrajectories(bootSamples(:, bb)));
        
        model       = bootModels{bb};
        startIdx    = (bb-1)*nSamples+1;
        endIdx      = bb*nSamples;
        
        for ii=1:d
            rawVals     = rawPhi(:, ii);
            discVals    = discPhi(:, ii);
            
            isObs       = ~isnan(discVals);
            
            for kk=1:nStates
                weights         = model.emission.T(kk, discVals(isObs), ii);
                allSamples(kk, ii, startIdx:endIdx)     = ...
                    randsample(rawVals(isObs), nSamples, true, weights);
            end
        end
    end
    
    obsStats        = quantile(allSamples, quantileValues, 3);
    obsStats        = permute(obsStats, [3, 1, 2]);
    clear allSamples;
    
    stateVector     = stateIdx(emissionPath);
    
    if isempty(plotIndices)
        nPlots          = ceil(d/4);
        plotIndices     = cell(nPlots, 1);
        for p=1:nPlots
            axesIdx         = (p-1)*4 + 1 : p*4;
            axesIdx         = axesIdx(axesIdx <= d & axesIdx ~= 17);
            plotIndices{p}  = mat2cell(colvec(axesIdx), ...
                ones(length(axesIdx), 1));
        end
    end
    
    nPlots          = length(plotIndices);
    nEdges          = length(edgeMatrices);
    
    for p=1:nPlots
        axes        = plotIndices{p};
        nAxes       = length(axes);
        
        fig         = figure('units', 'normalize', 'outerposition', ...
            [0 0 1 1]);
        subPlots    = cell(nAxes, 1);
        
        for aa=1:nAxes
            hSub    = subplot(nAxes, 1, aa);
            hold on;
            subPlots{aa}    = hSub;
            
            dimIdx          = axes{aa};
            nLines          = length(dimIdx);
            
            y       = reshape(obsStats(2, :, dimIdx), nStates, []);
            yNeg    = y - reshape(obsStats(1, :, dimIdx), nStates, []);
            yPos    = reshape(obsStats(3, :, dimIdx), nStates, []) - y;
            
            eb  = errorbar(repmat(stateVector, 1, nLines), ...
                y(emissionPath, :), yNeg(emissionPath, :), ...
                yPos(emissionPath, :), 'o', 'LineWidth', 1.5, ...
                'color', [0.4940 0.1840 0.5560]);
            
            for ee=1:nEdges
                edges       = edgeMatrices{ee};
                color       = edgeColors(ee, :);
                lineStyle   = edgeLines{ee};
                
                for ii=1:size(edges, 1)
                    stSeq   = edges(ii, :);
                    plot(stSeq, y(stSeq, :), lineStyle, 'Color', color, ...
                        'linewidth', 1.5);
                end
            end

            grid on; box on;
            
            hSub.FontSize   = 24;
        
            if aa < nAxes
                hSub.XTickLabel     = [];
            end
            
            hSub.Position(4)    = hSub.Position(4) + 0.05;
            hSub.Position(3)    = hSub.Position(3) + 0.04;
            
            if mod(aa, 2) == 0
                hSub.YAxisLocation  = 'right';
            end
            
            res         = max([yNeg; yPos])*0.1;
            lims        = [min(y-yNeg)-res max(y+yPos)+res];

            hSub.YLim       = lims;
            hSub.XTick      = min(stateVector):max(stateVector);
            
            legend(allData.phiNames(dimIdx), 'Interpreter', 'none');
            hSub.Legend.Location    = 'best';
        end
        
        linkaxes(cat(2, subPlots{:}), 'x');
        hSub.XLim       = [min(stateVector)-0.5 max(stateVector)+0.5];
    end
end

end