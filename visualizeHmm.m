function visualizeHmm(allData, trajectories, rawTrajectories, dx, ...
    age, model, varargin)

[   showAnart                                           , ...
    showReading                                         , ...
    showAge                                             , ...
    showApoe                                            , ...
    showTrans                                           , ...
    showEmpiricalTrans                                  , ...
    showDx                                              , ...
    showDurations                                       , ...
    showEmission                                        , ...
    filterMask                                          , ...
    emissionParams          ] = process_options(varargin, ...
    'showAnart'             , true                      , ...
    'showReading'           , true                      , ...
    'showAge'               , true                      , ...
    'showApoe'              , true                      , ...
    'showTrans'             , true                      , ...
    'showEmpiricalTrans'    , true                      , ...
    'showDx'                , true                      , ...
    'showDurations'         , false                     , ...
    'showEmission'          , true                      , ...
    'filterMask'            , true(size(trajectories))  , ...
    'emissionParams'        , {}                        );

trajectories        = trajectories(filterMask);
rawTrajectories     = rawTrajectories(filterMask);
dx                  = dx(filterMask);
age                 = age(filterMask);

stackedPhi      = cat(1, trajectories{:});
stackedDx       = cat(1, dx{:});
stackedAge      = cat(1, age{:});

paths           = cellfun(@(seq)hmmMap(model, seq'), trajectories, ...
    'UniformOutput', false);

nStates         = model.nstates;

if ~isfield(model.emission, 'parentState')
    parentState     = colvec(1:nStates);
else
    parentState     = model.emission.parentState;
end

nExpansions     = nStates/max(parentState) - 1;
[~, parentIdx]  = unique(parentState);
nParents        = length(parentIdx);

counts          = allData.counts;

stackedStates   = colvec(cat(2, paths{:}));
stackedStates   = parentState(stackedStates);
firstState      = colvec(cellfun(@(seq)seq(1), paths));
firstIdx        = cumsum([1; counts(1:end-1)]);

if showAnart
    anart           = allData.data.ANARTERR;
    anartTraj       = mat2cell(anart, counts);
    anartTraj       = anartTraj(filterMask);
    anart           = cat(1, anartTraj{:});
    
    figure; hold on;
    boxplot(anart(firstIdx), firstState);
    ax                  = gca;
    ax.XLabel.String    = 'HMM State';
    ax.YLabel.String    = '# ANART Errors';
    ax.FontSize         = 30;
    ax.XLim             = [0.2 nParents+0.8];
    ax.YGrid            = 'on';
end

if showReading
    reading         = allData.data.PTEDUCAT;
    reading         = mat2cell(reading, counts);
    reading         = reading(filterMask);
    reading         = cat(1, reading{:});
    
    reading         = discretize(reading, [0 15 18 inf]);
    
    freqMat         = accumarray([stackedStates(~isnan(reading)) ...
        reading(~isnan(reading))], 1);
    
    figure; hold on;
%     boxplot(reading(firstIdx), stackedStates(firstIdx));
    bar(freqMat);
    ax                  = gca;
    ax.XLabel.String    = 'HMM State';
    ax.YLabel.String    = 'Count';
    legend('Low', 'Intermediate', 'High');
    ax.FontSize         = 30;
    ax.XLim             = [0.2 nParents+0.8];
    ax.YGrid            = 'on';
end

if showAge
    figure; hold on;
    boxplot(stackedAge, stackedStates);
    ax                  = gca;
    ax.XLabel.String    = 'HMM State';
    ax.YLabel.String    = 'Age';
    ax.FontSize         = 30;
    ax.YGrid            = 'on';
end

if showApoe
    apgen           = allData.data{:, {'CAT_APGEN1_2_0', ...
        'CAT_APGEN1_3_0', 'CAT_APGEN1_4_0', 'CAT_APGEN2_2_0', ...
        'CAT_APGEN2_3_0', 'CAT_APGEN2_4_0'}};
    apgen           = mat2cell(apgen, counts);
    apgen           = apgen(filterMask);
    apgen           = cat(1, apgen{:});
    
    [apgen1, ~]     = find(apgen(:, 1:3)');
    [apgen2, ~]     = find(apgen(:, 4:6)');
    apgen           = [apgen1+1 apgen2+1];

    geneProfile     = nan(size(apgen, 1), 1);

    for ii=1:length(geneProfile)
        if isequal(apgen(ii, :), [2 2])
            geneProfile(ii)     = 1;
        elseif isequal(apgen(ii, :), [3 3])
            geneProfile(ii)     = 2;
        elseif isequal(apgen(ii, :), [4 4])
            geneProfile(ii)     = 3;
        elseif isequal(apgen(ii, :), [2 3])
            geneProfile(ii)     = 4;
        elseif isequal(apgen(ii, :), [3 4])
            geneProfile(ii)     = 5;
        elseif isequal(apgen(ii, :), [2 4])
            geneProfile(ii)     = 6;
        end
    end
    
    heatmap         = accumarray([stackedStates geneProfile], 1);
    
    figure; hold on;
    bar(heatmap);
    legend('ApoE 2,2', 'ApoE 3,3', 'ApoE 4,4', 'ApoE 2,3', ...
        'ApoE 3,4', 'ApoE 2,4')
    grid on
    box on
    ax  = gca;
    ax.XLim = [0.2 nParents+0.8];
    ax.FontSize = 30;
    ax.XLabel.String = 'HMM State';
    ax.YLabel.String = 'Count';
end

if showDx
    knownDx         = ~isnan(stackedDx);
    freqMat         = accumarray([stackedStates(knownDx) ...
        stackedDx(knownDx)], 1);
    
    figure; hold on;
    dxBar               = bar(freqMat);
    legend('NL', 'MCI', 'AD');
    grid on; box on;
    ax              = gca;
    dxBar(1).FaceColor  = ax.ColorOrder(end-2, :);
    dxBar(2).FaceColor  = ax.ColorOrder(1, :);
    dxBar(3).FaceColor  = ax.ColorOrder(2, :);
    ax.XLim         = [0.2 nParents+0.8];
    ax.FontSize     = 30;
    ax.XLabel.String    = 'HMM State';
    ax.YLabel.String    = 'Count';
end

if showTrans
    figure;
    title('Learned Transition Matrix');
    imagesc(model.A);
    ax          = gca;
    axis square;
    ax.XLabel.String    = 'Destination Stage';
    ax.YLabel.String    = 'Source Stage';
    colorbar;
end

if showEmpiricalTrans
    empiricalA          = countTransitions(paths, model.nstates);
    
    figure;
    title('Empirical Transition Matrix');    
    imagesc(empiricalA);
    ax          = gca;
    axis square;
    ax.XLabel.String    = 'Destination Stage';
    ax.YLabel.String    = 'Source Stage';
    ax.FontSize         = 32;
    ax.XTick            = min(parentIdx):2:max(parentIdx);
    ax.YTick            = min(parentIdx):2:max(parentIdx);
    colorbar;
end

if showDurations
    maxDuration         = 20;
    
    durationProb        = zeros(nParents, maxDuration);
    stateNames          = cell(nParents, 1);
    
    for ii=1:nParents
        subMatrix       = parentIdx(ii) + (0:nExpansions);
        
        if nExpansions > 0
            probs       = [diag(model.A(subMatrix, subMatrix), 1); ...
                model.A(subMatrix(end), subMatrix(end))*...
                ones(maxDuration-nExpansions, 1)];
        else
            probs       = model.A(subMatrix, subMatrix)*ones(maxDuration, 1);
        end
        
        durationProb(ii, :)     = (1-probs).*cumprod([1; probs(1:end-1)]);
        
        stateNames{ii}  = sprintf('K=%d', ii);
    end
    
    figure;
    imagesc(durationProb);
    ax          = gca;
    box on;
    xlabel('Duration');
    ylabel('Probability');
    ax.XTick        = 1:2:maxDuration;
    ax.XTickLabel   = ax.XTick*6;
    ax.YTick        = 1:nParents;
    ax.YTickLabel   = parentIdx;
    ax.FontSize     = 30;
    colorbar;
    
%     figure; hold on;
%     plot((1:maxDuration)*6, durationProb(1:end-1, :)', 'o--', ...
%         'LineWidth', 1.5);
%     legend(stateNames(1:end-1));
%     grid on; box on;
%     xlabel('Duration');
%     ylabel('Probability');
%     ax.FontSize     = 30;
end

if showEmission
    [   emissionPath                                            , ...
        quantileLims                                            , ...
        nBootSamples                                            , ...
        sampleObs                                               , ...
        edgeMatrices                                            , ...
        edgeColors                                              , ...
        plotIndices                                             , ...
        edgeLines           ] = process_options(emissionParams  , ...
        'emissionPath'      , 1:nStates                         , ...
        'quantileLims'      , []                                , ...
        'nBootSamples'      , 10000                             , ...
        'sampleObs'         , false                             , ...
        'edgeMatrices'      , []                                , ...
        'edgeColors'        , []                                , ...
        'plotIndices'       , []                                , ...
        'edgeLines'         , {}                                );
    
    rawPhi          = cat(1, rawTrajectories{:});
    
    d               = size(rawPhi, 2);
    stateVector     = parentIdx(emissionPath);
    
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
    
    if isempty(quantileLims)
        quantileLims    = getSampleQuantiles(model.emission, stackedPhi, ...
            rawPhi, parentIdx, stackedStates, 'nBootSamples', nBootSamples, ...
            'quantileVals', [0.25 0.5 0.75], ...
            'sampleObs', sampleObs);
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
            
            y       = reshape(quantileLims(2, :, dimIdx), nParents, []);
            yNeg    = y - reshape(quantileLims(1, :, dimIdx), nParents, []);
            yPos    = reshape(quantileLims(3, :, dimIdx), nParents, []) - y;
            
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