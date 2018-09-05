showDx              = false;
showDurations       = false;
showApoe            = false;
showAnart           = false;
showAge             = false;
showReading         = false;
showTrans           = false;
showEmission        = false;
showPathDist        = false;
showPathGraph       = false;
showTrajectories    = true;
showPathTimes       = false;
showPathConfusion   = false;
showSurvival        = false;
showSlopes          = false;
showInteraction     = false;
showDistributions   = false;
showBackwardElim    = false;

computeStats        = false;
sampleEmissions     = true;

quantileValues  = [0.25 0.5 0.75];
topKCount       = 20;

emStates        = 8:12;
edgeMatrices    = {[16 18; 18 21; 21 22], [17 20; 20 22], ...
    [19 21; 21 22]};
edgeColors      = [0.4660 0.6740 0.1880; ...
    0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250; 0 0.4470 0.7410];
edgeLines       = {'--', '-.', '-.'};

% plotIndices     = {[1 6 7 10], [8 9 11 13], [14, 15, 16, 18]};
plotIndices     = {[1 10], [8 18]};

isContinuous    = {[true, true, true, true], [true, true, false, false], ...
    [true, true, true, true]};
resolution      = {[-1, -1, -1, -1], [-1 -1 3 1], [4 2 2 1]};

emissionParams  = {'emissionPath', emStates, 'quantileValues', ...
    [0.005 0.5 0.995], 'nSamples', 1000, 'edgeMatrices', [], ...
    'edgeColors', edgeColors, 'edgeLines', edgeLines, ...
    'plotIndices', plotIndices};

if computeStats

baseYl                  = yl;
baseYu                  = yu;
baseTrajectories        = trajectories;
baseRawTrajectories     = rawTrajectories;
baseDx                  = dx;
baseAge                 = age;

nBoot           = size(bootSamples, 2);

instIndex       = mat2cell(colvec(1:size(allData.data, 1)), ...
    allData.counts);

allAnart        = cell(nBoot, 1);
allAge          = cell(nBoot, 1);
allReading      = cell(nBoot, 1);
allGenes        = cell(nBoot, 1);

allDxFreq       = cell(nBoot, 1);

allFirstState   = cell(nBoot, 1);
allPaths        = cell(nBoot, 1);

allLetterStateSets  = cell(nBoot, 1);
allUncTrans         = cell(nBoot, 1);
allEmpiricalTrans   = cell(nBoot, 1);

nStates             = hmmModels{1}.nstates;
stateIdx            = colvec(1:nStates);

letters         = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', ...
    'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', ...
    'x', 'y', 'z'};

transDuration       = cell(nBoot, 1);
adDuration          = cell(nBoot, 1);
uncensStatus        = cell(nBoot, 1);

maxCount            = max(allData.counts);

allFullStates       = cell(nBoot, 1);
    
parfor bb=1:nBoot
    sampleIdx           = bootSamples(:, bb);
    
    bootCounts          = allData.counts(sampleIdx);
    bootInstances       = cell2mat(instIndex(sampleIdx));
    
    stackedStates       = hmmStates{bb};
    model               = hmmModels{bb};
    
    dx                  = baseDx(sampleIdx);
    age                 = baseAge(sampleIdx);
    yu                  = baseYu(sampleIdx);
    
    uncensored          = cellfun(@(label, yupper)any(label == 2) && ...
        ~any(yupper == 999), dx, yu);
    
    uncensStatus{bb}    = arrayfun(@(ii)uncensored(ii)*...
        ones(bootCounts(ii), 1), colvec(1:length(bootCounts)), ...
        'uniformoutput', false);
    
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
    allAge{bb}          = age;
    
    allReading{bb}      = allData.data.PTEDUCAT(bootInstances);
    allReading{bb}      = allReading{bb}(firstIdx);
    
    stackedDx       = cell2mat(dx);
    knownDx         = ~isnan(stackedDx);
    allDxFreq{bb}   = accumarray([stackedStates(knownDx) ...
        stackedDx(knownDx)], 1);
    
    allGenes{bb}    = getGeneProfile(allData.data(bootInstances, :));
    allGenes{bb}    = allGenes{bb}(firstIdx);
    
    allPaths{bb}    = mat2cell(hmmStates{bb}, bootCounts);
    
    allEmpiricalTrans{bb}   = countTransitions(allPaths{bb}, nStates);
    allUncTrans{bb}         = countTransitions(allPaths{bb}(uncensored), ...
        nStates);
    
    paths               = allPaths{bb};
    timeToTrans         = cell(length(paths), 1);
    timeToAd            = cell(length(paths), 1);
        
    for p=1:length(paths)
        seq             = paths{p};
        transPoints     = find(diff(seq) > 0);
        
        if ~isempty(transPoints)
            transPoints     = [0; transPoints];
        end
        
        timeTillTrans       = nan(size(seq));
        
        for tp=1:length(transPoints)-1
            beginIdx        = transPoints(tp)+1;
            endIdx          = transPoints(tp+1);
            
            timeTillTrans(beginIdx:endIdx)  = (endIdx-beginIdx+1):-1:1;
        end
        
        timeToTrans{p}  = colvec(timeTillTrans);
        
        firstAd         = find(dx{p} == 3, 1, 'first');
        timeTillAd      = nan(size(seq));
        
        if ~isempty(firstAd) && all(dx{p}(firstAd+1:end) == 3)
            timeTillAd(1:firstAd-1)     = firstAd-1:-1:1;
            timeTillAd(firstAd:end)     = 0:-1:firstAd-length(seq);
        end
        
        timeToAd{p}                 = timeTillAd;
    end
        
    transDuration{bb}           = cat(1, timeToTrans{:});
    adDuration{bb}              = cat(1, timeToAd{:});
    
    allFullStates{bb}           = cell2mat(cellfun(@(seq)[seq; ...
        ones(maxCount-length(seq), 1)*(nStates+1)], paths, ...
        'uniformoutput', false));
end

allEmpiricalTrans           = cat(3, allEmpiricalTrans{:});
allUncTrans                 = cat(3, allUncTrans{:});

bootRawPhi      = arrayfun(@(bb)cell2mat(...
    baseRawTrajectories(bootSamples(:, bb))), colvec(1:nBoot), ...
    'uniformoutput', false);
bootDx          = arrayfun(@(bb)cell2mat(baseDx(bootSamples(:, bb))), ...
    colvec(1:nBoot), 'uniformoutput', false);

stackedBootPhi          = cat(1, bootRawPhi{:}); clear bootRawPhi;
stackedBootStates       = cat(1, hmmStates{:});
stackedBootDx           = cat(1, bootDx{:}); clear bootDx;
stackedBootTransTime    = cat(1, transDuration{:}); clear transDuration;
stackedBootAdTime       = cat(1, adDuration{:}); clear adDuration;
stackedBootUncStatus    = cell2mat(cat(1, uncensStatus{:})); ...
stackedFullStates       = cell2mat(allFullStates);
clear uncensStatus;

stackedBootTimeIndex    = repmat(colvec(1:maxCount), nTraj*nBoot, 1);
stateCountMatrix        = zeros(nStates+1, maxCount);
stateTransMatrix        = zeros(nStates+1, nStates+1, maxCount);

for cc=1:maxCount
    sourceMask          = stackedBootTimeIndex == cc;
    stateCountMatrix(:, cc)     = accumarray(...
        stackedFullStates(sourceMask), 1, [nStates+1 1]);
    
    if cc < maxCount
        destMask        = stackedBootTimeIndex == (cc+1);
        stateTransMatrix(:, :, cc)  = accumarray(...
            [stackedFullStates(sourceMask) stackedFullStates(destMask)], ...
            1, [nStates+1, nStates+1]);
    end
end

pathMembership      = cell(nBoot, 1);
stateMembership     = cell(nBoot, 1);

parfor bb=1:nBoot
    paths           = allPaths{bb};
    isPathA         = @(seq)any(ismember(seq, [8 10]));
    isPathB         = @(seq)any(ismember(seq, [9 11]));
%     isPathC         = @(seq)any(ismember(seq, [11 13]));
    membership      = ones(size(paths))*-1;
    membership(cellfun(isPathA, paths))     = 1;
    membership(cellfun(isPathB, paths))     = 2;
%     membership(cellfun(isPathC, paths))     = 3;
    
    pathMembership{bb}      = membership;
    
    stateMembership{bb}     = arrayfun(@(ii)ones(size(paths{ii}))*...
        membership(ii), colvec(1:length(membership)), 'uniformoutput', false);
end

end

if showPathTimes
    startStates     = [8 9];
    nPaths          = length(startStates);
    pathTimes       = cell(nBoot, nPaths);
    
    nSamples        = 200;
    T               = 50;
    
    parfor bb=1:nBoot
        model       = hmmModels{bb};
        for pp=1:nPaths
            times   = getHmmStepCounts(model, startStates(pp), ...
                model.nstates, 'nSamples', nSamples, ...
                'sampleTillReached', false, 'T', T);
            
            pathTimes{bb, pp}   = times;
        end
    end
    
    sampledTimes        = cell(nPaths, 1);
    
    for pp=1:nPaths
        sampledTimes{pp}    = cell2mat(pathTimes(:, pp));
    end
    
    times           = cat(2, sampledTimes{:});
    groups          = repmat(1:nPaths, nSamples*nBoot, 1);
    bIdx            = repmat(colvec(repmat(1:nBoot, nSamples, 1)), 1, nPaths);
    
    frequencies     = accumarray([colvec(times) colvec(groups) ...
        colvec(bIdx)], 1);
    
    probStats       = quantile(cumsum(frequencies/nSamples, 1), ...
        [0.25 0.5 0.75], 3);
    
    tLim            = 25;
    
    jitter          = repmat(linspace(-0.06, 0.06, nPaths), tLim, 1);
    
    figure; hold on;
    errorbar(repmat(colvec(1:tLim)*0.5, 1, nPaths) + jitter, ...
        probStats(1:tLim, :, 2), ...
        probStats(1:tLim, :, 2)-probStats(1:tLim, :, 1), ...
        probStats(1:tLim, :, 3)-probStats(1:tLim, :, 2), '--', ...
        'linewidth', 1.5);
    grid on; box on;
    ax              = gca;
    ylim([0 1]);
    xlim([0 (tLim)*0.5]);
    ax.XTick        = 0:2:26;
    ax.XLim         = [0 (tLim+1)/2];
    ax.FontSize     = 32;
    xlabel('Time in years');
    ylabel('Cumulative Probability');
    legend('Stage 8 -> Terminal Stage (Path A)', ...
        'Stage 9 -> Terminal Stage (Path B)');
end

if showAnart
    figure;
    sidedViolin(cell2mat(allAnart), cell2mat(colvec(allFirstState)), [], ...
        'showData', false, 'maxWidth', 0.35, ...
        'isContinuous', true, 'normType', 'area', ...
        'normWithinSplit', false);
    ax                  = gca;
    ax.XLabel.String    = 'Disease Stage';
    ax.YLabel.String    = 'Number of ANART Errors';
    ax.FontSize         = 32;
    ax.XLim             = [0.2 nStates+0.8];
    ax.YGrid            = 'on';
    box on;
end

if showPathConfusion
    allInstanceIndices          = cell2mat(instIndex(colvec(bootIdx)));
    
    freqMat         = accumarray([allInstanceIndices stackedBootStates], 1);
    
    [frequency, mostFrequent]       = max(freqMat, [], 2);
    
    confusionMat        = zeros(nStates);
    
    for k=1:nStates
        confusionMat(k, :)      = mean(normalize(...
            freqMat(mostFrequent == k, :), 2), 1);
    end
    
    figure;
    imagesc(confusionMat);
    colorbar;
    axis square;
    ax                  = gca;
    ax.FontSize         = 32;
end

if showAge
    age         = cell2mat(cat(1, allAge{:}));
    group       = nan(size(stackedBootStates));
    group(ismember(stackedBootStates, [9 10 11]))   = 1;
    group(ismember(stackedBootStates, [12 13]))     = 2;
    group(ismember(stackedBootStates, [14]))        = 3;
    
    side        = nan(size(stackedBootStates));
    
    pathMem             = cell2mat(cat(1, stateMembership{:}));
    
    side(ismember(pathMem, [1 2]))  = 1;
    side(pathMem == 3)              = 2;
    
    mask        = ~isnan(group) & ~isnan(age) & ~isnan(side);
    
    figure;
    ageViolins  = sidedViolin(age(mask), group(mask), side(mask), ...
        'showData', false, 'maxWidth', 0.35, ...
        'isContinuous', true, 'normType', 'area', ...
        'normWithinSplit', false);
    ax                  = gca;
    ax.YLabel.String    = 'Age';
    ax.FontSize         = 32;
    ax.XLim             = [0.4 3.6];
    ax.XTick            = 1:3;
    ax.YGrid            = 'on';
    ax.XTickLabel       = {'9/10  |   11  ', '12  |  13', '14A  |  14B'};
    box on;
    
    legend([ageViolins{1, 1} ageViolins{1, 2}], 'Path A', 'Path B');
    
    for vv=1:size(ageViolins, 1)
        ageViolins{vv, 1}.FaceColor     = ax.ColorOrder(end-2, :);
        ageViolins{vv, 2}.FaceColor     = ax.ColorOrder(2, :);
    end
end

if showReading
    reading         = cell2mat(allReading);
    
    pMem            = cell2mat(pathMembership);
    bootSetIdx      = colvec(repmat(1:nBoot, 1, nSeries));
    
    group           = nan(size(pMem));
    group(ismember(pMem, [1]))        = 1;
    group(ismember(pMem, [2]))        = 2;
    
    mask            = ~isnan(reading) & ~isnan(group) & (group ~= -1);
    
    discreteReading     = discretize(reading, [0 15 18 inf]);
    readingFreq         = accumarray([group(mask) discreteReading(mask)...
        bootSetIdx(mask)], 1);
    freqStats       = quantile(normalize(readingFreq, 2), ...
        [0.25 0.5 0.75], 3);
    
    figure; hold on;
    [hBars, hErrors]    = barwitherr(freqStats(:, :, 2), ...
        cat(3, freqStats(:, :, 2) - freqStats(:, :, 1), ...
        freqStats(:, :, 3) - freqStats(:, :, 2)));
    ax                  = gca;
    ax.YLabel.String    = 'Probability';
    legend('Low', 'Intermediate', 'High');
    ax.FontSize         = 32;
    ax.YGrid            = 'on';
    hBars(1).FaceColor  = ax.ColorOrder(2, :);
    hBars(2).FaceColor  = ax.ColorOrder(1, :);
    hBars(3).FaceColor  = ax.ColorOrder(5, :);
    hBars(1).FaceAlpha  = 0.7;
    hBars(2).FaceAlpha  = 0.7;
    hBars(3).FaceAlpha  = 0.7;
    box on;
    ax.XTick            = [1 2];
    ax.XTickLabel       = {'Path A', 'Path B'};
    
    allPval             = nan(nBoot, 1);
    bootMedianRead      = nan(nBoot, 2);
    
    parfor bb=1:nBoot
        mask            = pathMembership{bb} ~= -1 & ~isnan(allReading{bb});
        discReading     = discretize(allReading{bb}, [0 15 18 inf]);
        [~, ~, allPval(bb)]     = crosstab(pathMembership{bb}(mask), ...
            discReading(mask));
        
        bootMedianRead(bb, :)   = accumarray(pathMembership{bb}(mask), ...
            allReading{bb}(mask), [], @median);
    end
    
    [table, chi2, pvalue]   = crosstab(group(mask), discreteReading(mask));
    [~, medianP, medianChi2]    = chi2cont(median(readingFreq, 3));
end

if showApoe
    geneProfile     = cell2mat(allGenes);
    e4Count         = nan(size(geneProfile));
    % store count+1 for accumarray because count >= 0
    e4Count(ismember(geneProfile, [1 2 4]))     = 1;
    e4Count(ismember(geneProfile, [5 6]))       = 2;
    e4Count(ismember(geneProfile, 3))           = 3;
    
    pMem            = cell2mat(pathMembership);
    bootSetIdx      = colvec(repmat(1:nBoot, 1, nSeries));
    
    group           = nan(size(pMem));
    group(ismember(pMem, [1]))        = 1;
    group(ismember(pMem, 2))            = 2;
    
    mask            = ~isnan(e4Count) & ~isnan(group) & (group >= 0);
    
    geneFreq        = accumarray([group(mask) e4Count(mask)...
        bootSetIdx(mask)], 1);
    freqStats       = quantile(normalize(geneFreq, 2), ...
        [0.25 0.5 0.75], 3);
    
    figure; hold on;
    [hBars, hErrors]    = barwitherr(freqStats(:, :, 2), ...
        cat(3, freqStats(:, :, 2) - freqStats(:, :, 1), ...
        freqStats(:, :, 3) - freqStats(:, :, 2)));
    ax                  = gca;
    ax.YLabel.String    = 'Probability';
    legend('\epsilon4 Non-Carrier', '\epsilon4 Single Carrier', ...
        '\epsilon4 Double Carrier');
    ax.FontSize         = 32;
    ax.YGrid            = 'on';
    hBars(3).FaceColor  = ax.ColorOrder(2, :);
    hBars(2).FaceColor  = ax.ColorOrder(1, :);
    hBars(1).FaceColor  = ax.ColorOrder(5, :);
    hBars(1).FaceAlpha  = 0.7;
    hBars(2).FaceAlpha  = 0.7;
    hBars(3).FaceAlpha  = 0.7;
    box on;
    ax.XTick            = [1 2];
    ax.XTickLabel       = {'Path A', 'Path B'};
    
    allPval             = nan(nBoot, 1);
    
    parfor bb=1:nBoot
        eCount          = nan(size(allGenes{bb}));
        % store count+1 for accumarray because count >= 0
        eCount(ismember(allGenes{bb}, [1 2 4]))     = 1;
        eCount(ismember(allGenes{bb}, [5 6]))       = 2;
        eCount(ismember(allGenes{bb}, 3))           = 3;
        
        mask            = pathMembership{bb} ~= -1 & ~isnan(eCount);
        [~, ~, allPval(bb)]     = crosstab(pathMembership{bb}(mask), ...
            eCount(mask));
    end
    
    [table, chi2, pvalue]   = crosstab(group(mask), e4Count(mask));
    [~, medianP, medianChi2]    = chi2cont(median(geneFreq, 3));
end

if showInteraction
    allPval         = nan(nBoot, 1);
    
    parfor bb=1:nBoot
        eCount          = nan(size(allGenes{bb}));
        % store count+1 for accumarray because count >= 0
        eCount(ismember(allGenes{bb}, [1 2 4]))     = 1;
        eCount(ismember(allGenes{bb}, [5 6]))       = 2;
        eCount(ismember(allGenes{bb}, 3))           = 3;
        
        discReading     = discretize(allReading{bb}, [0 15 18 inf]);
        
        mask            = pathMembership{bb} == 2 & ~isnan(eCount) & ...
            ~isnan(discReading);
        
        [tbl, ~, allPval(bb)]     = crosstab(eCount(mask), ...
            discReading(mask));
    end
end

if showDx
    normFreq        = cellfun(@(seq)normalize(seq, 1), allDxFreq, ...
        'uniformoutput', false);
    freqMat         = cat(3, normFreq{:});
    freqStats       = quantile(freqMat, quantileValues, 3);
    
    figure; hold on;
    [hBars, hErrors]    = barwitherr(freqStats(:, :, 2), ...
        cat(3, freqStats(:, :, 2) - freqStats(:, :, 1), ...
        freqStats(:, :, 3) - freqStats(:, :, 2)));
    ax                  = gca;
    ax.XLabel.String    = 'Disease Stage';
    ax.YLabel.String    = 'Probability of assignment to stage';
    legend('NL', 'MCI', 'AD');
    ax.FontSize         = 32;
    ax.XLim             = [0.2 nStates+0.8];
    ax.YGrid            = 'on';
    hBars(1).FaceColor  = ax.ColorOrder(end-2, :);
    hBars(2).FaceColor  = ax.ColorOrder(1, :);
    hBars(3).FaceColor  = ax.ColorOrder(2, :);
    box on;
    
    nlIdx       = cell2mat(arrayfun(@(ii)ones(round(...
        freqStats(ii, 1, 2)*1000), 1)*ii, colvec(1:size(freqStats, 1)), ...
        'UniformOutput', false));
    mciIdx      = cell2mat(arrayfun(@(ii)ones(round(...
        freqStats(ii, 2, 2)*1000), 1)*ii, colvec(1:size(freqStats, 1)), ...
        'UniformOutput', false));
    adIdx       = cell2mat(arrayfun(@(ii)ones(round(...
        freqStats(ii, 3, 2)*1000), 1)*ii, colvec(1:size(freqStats, 1)), ...
        'UniformOutput', false));
    
    figure; hold on;
    bNl         = bar(freqStats(:, 1, 2));
    bMci        = bar(freqStats(:, 2, 2));
    bAd         = bar(freqStats(:, 3, 2));
    ax          = gca;
    bNl.FaceColor       = ax.ColorOrder(end-2, :);
    bMci.FaceColor      = ax.ColorOrder(1, :);
    bAd.FaceColor       = ax.ColorOrder(2, :);
    bNl.EdgeColor       = ax.ColorOrder(end-2, :);
    bMci.EdgeColor      = ax.ColorOrder(1, :);
    bAd.EdgeColor       = ax.ColorOrder(2, :);
    bNl.FaceAlpha       = 0.6;
    bMci.FaceAlpha      = 0.6;
    bAd.FaceAlpha       = 0.6;
    bNl.BarWidth        = 1;
    bMci.BarWidth       = 1;
    bAd.BarWidth        = 1;
    errorbar(repmat(colvec(1:nStates), 1, 3), freqStats(:, :, 2), ...
        freqStats(:, :, 2)-freqStats(:, :, 1), ...
        freqStats(:, :, 3)-freqStats(:, :, 2), '.k', 'linewidth', 1.2);
    ax.XLabel.String    = 'Disease Stage';
    ax.YLabel.String    = 'Probability';
    legend('NL', 'MCI', 'AD');
    ax.FontSize         = 32;
    ax.XLim             = [0.2 nStates+0.8];
    ax.YGrid            = 'on';
    box on;
end

if showPathGraph
    subNodes        = 8:12;
    
    confA           = quantile(allEmpiricalTrans, [0.25 0.5 0.75], 3);
    medianA         = confA(:, :, 2);
    
    confA           = confA(subNodes, subNodes, :);
    
    subGraph                            = medianA(subNodes, subNodes);
    subGraph(subGraph < 25)             = 0;
    subGraph(1:length(subGraph)+1:end)  = 0;
    
    yScore          = stackedBootPhi(:, strcmpi(allData.phiNames, 'ADAS'));
    yMask           = ~isnan(yScore);
    yLow            = accumarray(stackedBootStates(yMask), ...
        yScore(yMask), [], @(s)quantile(s, 0.25));
    yMedian         = accumarray(stackedBootStates(yMask), ...
        yScore(yMask), [], @(s)quantile(s, 0.5));
    yHigh           = accumarray(stackedBootStates(yMask), ...
        yScore(yMask), [], @(s)quantile(s, 0.75));
    
    xScore          = stackedBootPhi(:, strcmpi(allData.phiNames, 'AMYLOID'));
    xMask           = ~isnan(xScore);
    xLow            = accumarray(stackedBootStates(xMask), ...
        xScore(xMask), [], @(s)quantile(s, 0.25));
    xMedian         = accumarray(stackedBootStates(xMask), ...
        xScore(xMask), [], @(s)quantile(s, 0.5));
    xHigh           = accumarray(stackedBootStates(xMask), ...
        xScore(xMask), [], @(s)quantile(s, 0.75));
    
    G               = digraph(subGraph);
    
    MP              = get(0, 'MonitorPositions');
    
    fig                 = figure('units', 'pixels', 'position', MP(1, :));
    hold on;
    ax                  = gca;
    ax.FontSize         = 32;
    ax.XTick            = 1.25:0.05:1.45;
    ax.XLabel.String        = 'Amyloid SUVR';
    ax.YLabel.String        = 'ADAS-Cog Score';
    drawnow;
    
    dxCounts            = sum(cat(3, allDxFreq{:}), 3);
    colorScaler         = normalize(dxCounts(:, 2:3), 2);
    colorScaler         = colorScaler(subNodes, end-1);
    
    cmap                = viridis;
    hsv                 = rgb2hsv(cmap);
    colours             = interp1(linspace(0, 1, size(cmap,1)), hsv, ...
        colorScaler);
    colours             = hsv2rgb(colours);
    
    lWidth      = 30*G.Edges.Weight/max(G.Edges.Weight);
    
    graphPlot   = plot(G, 'edgelabel', G.Edges.Weight, ...
        'linewidth', lWidth, 'nodelabel', arrayfun(@num2str, subNodes, ...
        'uniformoutput', false), 'XData', xMedian(subNodes), ...
        'YData', yMedian(subNodes), 'NodeColor', [1 1 1], ...
        'Marker', 'o');
    
    graphPlot.MarkerSize    = 80;
    graphPlot.ArrowSize     = 25;
    graphPlot.NodeLabel     = {};
    graphPlot.EdgeLabel     = {};
    graphPlot.EdgeColor     = [0.8903 0.3480 0.1983];
    graphPlot.EdgeAlpha     = 0.35;
    
    for nn=1:size(G.Nodes, 1)
        s       = scatter(graphPlot.XData(nn), graphPlot.YData(nn), 6500, ...
            colours(nn, :), 'filled');
        s.MarkerFaceAlpha       = 0.7;
        text(graphPlot.XData(nn), graphPlot.YData(nn), ...
            num2str(subNodes(nn)), 'FontSize', 32, 'horizontalalignment', ...
            'center');
    end
    
    nEdges          = size(G.Edges, 1);
    edgeLabels      = cell(nEdges, 1);
    
    ax.Units        = 'inches';
    
    xScaler         = ax.Position(3)/(ax.XLim(2)-ax.XLim(1));
    yScaler         = ax.Position(4)/(ax.YLim(2)-ax.YLim(1));
    
    for ee=1:nEdges
        src         = G.Edges.EndNodes(ee, 1);
        dest        = G.Edges.EndNodes(ee, 2);
        
        dY          = graphPlot.YData(dest)-graphPlot.YData(src);
        dX          = graphPlot.XData(dest)-graphPlot.XData(src);
        slope       = (dY*yScaler)/(dX*xScaler);
        
        xPos        = mean([graphPlot.XData(src) graphPlot.XData(dest)]);
        yPos        = mean([graphPlot.YData(src) graphPlot.YData(dest)]);
        angle       = atand(slope);
        
        if src == 4 && dest == 5
            xPos    = xPos + 0.003
            yPos    = yPos + 0.001;
        end
        
        t           = text(xPos, yPos+0.26, ...
            sprintf('%0.0f (%0.0f%s%0.0f)', G.Edges.Weight(ee), ...
            confA(src, dest, 1), char(8211), confA(src, dest, 3)), ...
            'fontsize', 28, 'horizontalalignment', 'center', ...
            'verticalalignment', 'bottom');
        t.Rotation  = angle;
        
        edgeLabels{ee}  = t;
    end
    
    box on;
    tightfig;
end

if showSurvival
    allInstanceIndices      = cell2mat(instIndex(colvec(bootIdx)));
    
    convtime            = allData.data.CONVTIME_U(allInstanceIndices);
    censored            = convtime == 999;
    convtime(censored)  = allData.data.CONVTIME_L(...
        allInstanceIndices(censored));
    
    mask                = ismember(stackedBootStates, [8 9]) & ...
        (convtime >= 0);
    
    groupMasks          = bsxfun(@and, ...
        [stackedBootStates == 8 stackedBootStates== 9], mask);
    
    [f8, x8, fl8, fu8]      = ecdf(convtime(groupMasks(:, 1)), ...
        'censoring', censored(groupMasks(:, 1)), ...
        'function', 'survivor');
    [f9, x9, fl9, fu9]      = ecdf(convtime(groupMasks(:, 2)), ...
        'censoring', censored(groupMasks(:, 2)), ...
        'function', 'survivor');
%     [f10, x10, fl10, fu10]  = ecdf(convtime(groupMasks(:, 3)), ...
%         'censoring', censored(groupMasks(:, 3)), 'function', 'survivor');
%     [f11, x11, fl11, fu11]  = ecdf(convtime(groupMasks(:, 4)), ...
%         'censoring', censored(groupMasks(:, 4)), 'function', 'survivor');
    figure; hold on;
    stairs(x8, f8, 'linewidth', 1.5);
    stairs(x9, f9, 'linewidth', 1.5);
%     stairs(x10, f10, 'linewidth', 1.5);
%     stairs(x11, f11, 'linewidth', 1.5);
    legend('Path A', 'Path B');
    ax              = gca;
    ax.FontSize     = 32;
    ax.YTick        = 0:0.1:1;
    ax.YTickLabels  = 0:10:100;
    ax.XTick        = 0:12:120;
    ax.XTickLabel   = 0:1:10;
    ax.YLabel.String    = sprintf('%s \n %s', ...
        '1 - Probability of progression to', ' Alzhimer''s disease (%)');
    ax.XLabel.String    = 'Time (years)';
    box on;
    
    coxMask             = any(groupMasks(:, [1 2]), 2);
    [b, ~, H, stats]    = coxphfit(double(groupMasks(coxMask, [2])), ...
        convtime(coxMask), 'censoring', censored(coxMask));
    
    bootStats           = cell(nBoot, 1);
    
    parfor bb=1:nBoot
        indices         = cell2mat(instIndex(bootIdx(:, bb)));
    
        convtime            = allData.data.CONVTIME_U(indices);
        censored            = convtime == 999;
        convtime(censored)  = allData.data.CONVTIME_L(...
            indices(censored));
    
        states              = cell2mat(allPaths{bb});
        mask                = ismember(states, [8 9]) & ...
            (convtime >= 0);
    
        groupMasks          = bsxfun(@and, ...
            [states == 8 states == 9], mask);
        
        coxMask             = any(groupMasks, 2);
        [~, ~, ~, bootStats{bb}]        = coxphfit(double(groupMasks(...
            coxMask, [2])), convtime(coxMask), 'censoring', ...
            censored(coxMask));
    end
    
    allBeta     = cell2mat(cellfun(@(seq)rowvec(seq.beta), bootStats, ...
        'UniformOutput', false));
    allPval     = cell2mat(cellfun(@(seq)rowvec(seq.p), bootStats, ...
        'UniformOutput', false));
end

if showBackwardElim
    stackedPhi          = cell2mat(rawTrajectories);
    
    catPhiIdx           = 17;
    dummyPhi            = stackedPhi(:, ...
        setdiff(1:size(stackedPhi, 2), catPhiIdx));
    
    phiNames            = allData.phiNames(setdiff(1:size(stackedPhi, 2), ...
        catPhiIdx));
    
    for p=1:length(catPhiIdx)
        dPhi        = dummyvar(stackedPhi(:, catPhiIdx(p)));
        dummyPhi    = cat(2, dummyPhi, dPhi(:, 1:end-1));
        
        uVals       = unique(stackedPhi(:, catPhiIdx(p)));
        uVals       = uVals(~isnan(uVals));
        for u=1:length(uVals)-1
            phiNames    = cat(2, phiNames, sprintf('%s_%d', ...
                allData.phiNames{catPhiIdx(p)}, uVals(u)));
        end
    end
    
    dummyTraj           = mat2cell(dummyPhi, counts);
    ridTraj             = mat2cell(rid, counts);
    
    x                   = [];
    y                   = [];
    pid                 = [];
    groups              = [];
    
    runningCount        = 0;
    
    for bb=1:nBoot
        bootPhi         = dummyTraj(bootIdx(:, bb));
        bootCts         = cellfun(@(seq)size(seq, 1), bootPhi);
        
        bootLabels      = cell2mat(stateMembership{bb}) == 2;
        bootZ           = cell2mat(allPaths{bb});
        idx             = arrayfun(@(ii)ones(bootCts(ii), 1)*...
            (runningCount+ii), colvec(1:length(bootCts)), ...
            'UniformOutput', false);
        idx             = cell2mat(idx);
        
        runningCount    = runningCount + length(bootPhi);
        
        bootPhi         = cell2mat(bootPhi);
        
        mask            = ~any(isnan(bootPhi), 2) & (bootLabels ~= -1) & ...
            (bootZ >= 8);
        
        bootX           = bootPhi(mask, :);
        bootY           = bootLabels(mask, :);
        bootPid         = idx(mask);
        
        [~, firstIdx]   = unique(bootPid);
        bootGrp         = bootY(firstIdx);
        
        x           = cat(1, x, bootX);
        y           = cat(1, y, bootY);
        pid         = cat(1, pid, bootPid);
        groups      = cat(1, groups, bootGrp);
    end
    
    hofolds             = get_folds('lkpo'          , ...
        'groups'        , groups                    , ...
        'rid'           , pid                       , ...
        'patients'      , unique(pid)               , ...
        'nFolds'        , 5                         );
    hofolds.type        = 'lopo';

    cvfolds             = cell(hofolds.NumTestSets, 1);
    for f=1:hofolds.NumTestSets
        trainIdx        = hofolds.training(:, f);
        trainPatMask    = ismember(unique(pid), unique(pid(trainIdx)));

        cvfolds{f}      = get_folds('lkpo'          , ...
            'groups'    , groups(trainPatMask)      , ...
            'rid'       , pid(trainIdx)             , ...
            'nFolds'    , 5                         );
    end
    
    x               = zScale(x);

    [backwardError, backwardModels, backwardHyper, droppedPhi]  = ...
        backwardFeatureElimination(x, y, hofolds, cvfolds, ...
        'lambdaRange', 10.^linspace(-5, 2, 10), 'minFeatures', 1, ...
        'phiNames', phiNames);
    
%     backwardError       = cell(nBoot, 1);
%     droppedPhi          = cell(nBoot, 1);
%     
%     for bb=1:nBoot
%         stackedPhi      = cell2mat(dummyTraj(bootIdx(:, bb)));
%         stackedLabels   = cell2mat(stateMembership{bb});
%         stackedRid      = cell2mat(ridTraj(bootIdx(:, bb)));
%         
%         mask        = ~any(isnan(stackedPhi), 2) & (stackedLabels ~= -1);
%         
%         x           = stackedPhi(mask, :);
%         y           = stackedLabels(mask);
%         pid         = stackedRid(mask);
%         
%         [~, order]      = sort(pid);
%         x               = x(order, :);
%         y               = y(order);
%         pid             = pid(order);
%         
%         [~, firstIdx]       = unique(pid);
%         grp                 = y(firstIdx);
%         
%         hofolds             = get_folds('lopo'          , ...
%             'groups'        , grp                       , ...
%             'rid'           , pid                       , ...
%             'patients'      , unique(pid)               );
%         
%         cvfolds         = cell(hofolds.NumTestSets, 1);
%         for f=1:hofolds.NumTestSets
%             trainIdx        = hofolds.training(:, f);
%             trainPatMask    = ismember(unique(pid), unique(pid(trainIdx)));
%             
%             cvfolds{f}      = get_folds('lkpo'          , ...
%                 'groups'    , grp(trainPatMask)         , ...
%                 'rid'       , pid(trainIdx)             , ...
%                 'nFolds'    , 3                         );
%         end
%         
%         [backwardError{bb}, ~, ~, droppedPhi{bb}]  = ...
%             backwardFeatureElimination(stackedPhi(mask, :), ...
%             stackedLabels(mask), hofolds, cvfolds, ...
%             'lambdaRange', 10.^linspace(-5, 2, 10), 'minFeatures', 1);
%     end
end

if showDistributions
    phiIndices          = 1:19;
    nPhi                = length(phiIndices);
    
    pathAStats          = zeros(nPhi, 3, 3, nBoot);
    pathBStats          = zeros(nPhi, 3, 3, nBoot);
    pathAMean           = zeros(nPhi, 3, nBoot);
    pathBMean           = zeros(nPhi, 3, nBoot);
    pathARank           = zeros(nPhi, 3, nBoot);
    pathBRank           = zeros(nPhi, 3, nBoot);
    
    sA      = [8, 10, 12];
    sB      = [9, 11, 12];
    
    parfor bb=1:nBoot
        sMem            = cell2mat(stateMembership{bb});
        z               = cell2mat(allPaths{bb});
        
        sPhi            = cell2mat(rawTrajectories(bootIdx(:, bb)));
        
        pS      = zeros(nPhi, 3, 3);
        mu      = zeros(nPhi, 3);
        r       = zeros(nPhi, 3);
        
        for ss=1:length(sA)
            pS(:, :, ss)        = quantile(sPhi(sMem == 1 & ...
                z == sA(ss), phiIndices), [0.25, 0.5, 0.75], 1)';
            mu(:, ss)           = nanmean(sPhi(sMem == 1 & z == sA(ss), ...
                phiIndices));
            for pp=1:length(phiIndices)
                r(pp, ss)       = nanmean(tiedrank(sPhi(sMem == 1 & ...
                    z == sA(ss), phiIndices(pp))));
            end
        end
        pathAStats(:, :, :, bb)         = pS;
        pathAMean(:, :, bb)             = mu;
        pathARank(:, :, bb)             = r;
        
        pS      = zeros(nPhi, 3, 3);
        mu      = zeros(nPhi, 3);
        r       = zeros(nPhi, 3);
        for ss=1:length(sB)
            pS(:, :, ss)        = quantile(sPhi(sMem == 2 & ...
                z == sB(ss), phiIndices), [0.25, 0.5, 0.75], 1)';
            mu(:, ss)           = nanmean(sPhi(sMem == 2 & z == sB(ss), ...
                phiIndices));
            for pp=1:length(phiIndices)
                r(pp, ss)       = nanmean(tiedrank(sPhi(sMem == 2 & ...
                    z == sB(ss), phiIndices(pp))));
            end
        end
        pathBStats(:, :, :, bb)         = pS;
        pathBMean(:, :, bb)             = mu;
        pathBRank(:, :, bb)             = r;
    end
    
    medA        = squeeze(pathAStats(:, 2, :, :));
    medB        = squeeze(pathBStats(:, 2, :, :));
    
    phiInterest     = [1, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19];
    
    medVals         = zeros(length(phiInterest), 3, 2);
    iqrVals         = zeros(length(phiInterest), 3, 2);
    sigVals         = zeros(length(phiInterest), 3);
    prctlVals       = zeros(length(phiInterest), 3);
    rankVals        = zeros(length(phiInterest), 3);
    
    for pp=1:length(phiInterest)
        pIdx        = phiInterest(pp);
        for ss=1:3
            medVals(pp, ss, 1)      = nanmedian(colvec(medA(pIdx, ss, :)));
            medVals(pp, ss, 2)      = nanmedian(colvec(medB(pIdx, ss, :)));
            
            iqrVals(pp, ss, 1)      = iqr(colvec(medA(pIdx, ss, :)));
            iqrVals(pp, ss, 2)      = iqr(colvec(medB(pIdx, ss, :)));
            
            sigVals(pp, ss)         = sum((medA(pIdx, ss, :) - ...
                medB(pIdx, ss, :)) > 0)/nBoot;
            
            medDiff         = colvec(medA(pIdx, ss, :) - medB(pIdx, ss, :));
            [~, valIdx]     = min(abs(medDiff));
            ranks           = tiedrank(medDiff);
            
            prctlVals(pp, ss)       = ranks(valIdx)/length(ranks);
            rankVals(pp, ss)        = sum(colvec(pathARank(pIdx, ss, :) - ...
                pathBRank(pIdx, ss, :)) > 0)/nBoot;
            
%             sigVals(pp, ss)     = ranksum(colvec(medA(pIdx, ss, :)), ...
%                 colvec(medB(pIdx, ss, :)));
        end
    end
    
    prctlVals       = min(prctlVals, 1-prctlVals);
    rankVals        = min(rankVals, 1-rankVals);
    sigVals         = min(sigVals, 1-sigVals);
end

if showSlopes
%     allInstanceIndices      = cell2mat(instIndex(colvec(bootIdx)));
%     
%     followup            = mat2cell(allData.data.FOLLOWUP, ...
%         allData.counts);
%     followup            = cellfun(@(seq)seq(1)-seq, followup, ...
%         'uniformoutput', false);
%     followup            = cell2mat(followup);
%     followup            = followup(allInstanceIndices);
%     
%     pathMem             = cell2mat(cat(1, stateMembership{:}));
%     
%     phiIdx              = strcmpi(allData.phiNames, 'TAU');
%     maskA               = ismember(pathMem, [1 2]) & ...
%         ismember(stackedBootStates, [9 10 12 14]) & ...
%         ~isnan(stackedBootPhi(:, phiIdx));
%     maskB               = ismember(pathMem, 3) & ...
%         ismember(stackedBootStates, [11 13 14]) & ...
%         ~isnan(stackedBootPhi(:, phiIdx));
%     
%     xxA         = nan(size(stackedBootStates));
%     xxB         = nan(size(stackedBootStates));
%     
%     xxA(ismember(stackedBootStates, [9 10]))        = 1;
%     xxA(ismember(stackedBootStates, [12]))          = 2;
%     xxA(ismember(stackedBootStates, [14]))          = 3;
%     
%     xxB(ismember(stackedBootStates, [11]))          = 1;
%     xxB(ismember(stackedBootStates, [13]))          = 2;
%     xxB(ismember(stackedBootStates, [14]))          = 3;
%     
%     modelA      = fitlm(xxA(maskA), ...
%         zscore(stackedBootPhi(maskA, phiIdx)));
%     modelB      = fitlm(xxB(maskB), ...
%         zscore(stackedBootPhi(maskB, phiIdx)));
    
    phiIdx          = strcmpi(allData.phiNames, 'TAU');
    modelsA         = cell(nBoot, 1);
    modelsB         = cell(nBoot, 1);
    
    parfor bb=1:nBoot
        states          = cell2mat(allPaths{bb});
        membership      = cell2mat(stateMembership{bb});
        
        phiMat          = cell2mat(rawTrajectories(bootIdx(:, bb)));
        
        maskA           = ismember(membership, 1) & ...
            ismember(states, [8 10 12]) & ~isnan(phiMat(:, phiIdx));
        
        maskB           = ismember(membership, 2) & ...
            ismember(states, [9 11 12]) & ~isnan(phiMat(:, phiIdx));
        
        xxA         = nan(size(states));
        xxB         = nan(size(states));
    
        xxA(ismember(states, 8))        = 1;
        xxA(ismember(states, 10))       = 2;
        xxA(ismember(states, 12))       = 3;
        
        xxB(ismember(states, 9))        = 1;
        xxB(ismember(states, 11))       = 2;
        xxB(ismember(states, 12))       = 3;
        
        modelsA{bb}         = fitlm(xxA(maskA), ...
            phiMat(maskA, phiIdx));
        modelsB{bb}         = fitlm(xxB(maskB), ...
            phiMat(maskB, phiIdx));
    end
    
    slopesA         = cellfun(@(mdl)mdl.Coefficients{'x1', 'Estimate'}, ...
        modelsA);
    slopesB         = cellfun(@(mdl)mdl.Coefficients{'x1', 'Estimate'}, ...
        modelsB);
    pA              = cellfun(@(mdl)mdl.Coefficients{'x1', 'pValue'}, ...
        modelsA);
    pB              = cellfun(@(mdl)mdl.Coefficients{'x1', 'pValue'}, ...
        modelsB);
end

if showTrajectories
    collapseStates      = 1:7;
    expandedStates      = 8:nStates+1;
    
    nMappedStates           = length(expandedStates)+1;
    
    mappedCountMatrix       = [sum(stateCountMatrix(collapseStates, :), 1); ...
        stateCountMatrix(expandedStates, :)];
    mappedTransMatrix       = zeros(nMappedStates, nMappedStates, ...
        size(stateTransMatrix, 3));
    
    for t=1:size(stateTransMatrix, 3)
        mappedTransMatrix(:, :, t)      = [...
            sum(sum(stateTransMatrix(collapseStates, collapseStates, t))) ...
            sum(stateTransMatrix(collapseStates, expandedStates, t), 1); ...
            sum(stateTransMatrix(expandedStates, collapseStates, t), 2) ...
            stateTransMatrix(expandedStates, expandedStates, t)];
    end
    
    subNodes        = [7 1:6];
    tLim            = 10;
    dxCounts        = median(cat(3, allDxFreq{:}), 3);
    
    dxCounts        = [sum(dxCounts(collapseStates, 2:3), 1); ...
        dxCounts(collapseStates(end)+1:end, 2:3)];
    
    colorScaler         = normalize(dxCounts, 2);
    colorScaler         = colorScaler(:, end-1);
    
    plotTrajectories(mappedCountMatrix, mappedTransMatrix, subNodes, ...
        colorScaler, 'T', tLim, 'minTrans', 0.05, 'nodeLabels', ...
        {'13', '1-7', '8', '9', '10', '11', '12'}, ...
        'censColour', ones(1, 3)*0.5, 'arrowColour', [0 0 0]);
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
        'quantileValues'    , [0.01 0.5 0.99]                   , ...
        'nSamples'          , 1000                              , ...
        'edgeMatrices'      , []                                , ...
        'edgeColors'        , []                                , ...
        'plotIndices'       , []                                , ...
        'edgeLines'         , {}                                );
    
    d               = size(baseRawTrajectories{1}, 2);
    
    if sampleEmissions
        allSamples      = stackedBootPhi;
        zSample         = stackedBootStates;
        
        allSamples(:, 1)    = allSamples(:, 1)*100;
        
        obsStats        = quantile(allSamples, quantileValues, 1);
    end
    
    pathMem             = cell2mat(cat(1, stateMembership{:}));
    
    side                = nan(size(zSample));
    group               = nan(size(zSample));
    
    group(ismember(stackedBootStates, [8, 9]))         = 1;
    group(ismember(stackedBootStates, [10, 11]))       = 2;
    group(ismember(stackedBootStates, 12))             = 3;
    
%     group(ismember(stackedBootStates, [9, 10]))         = 1;
%     group(ismember(stackedBootStates, [11]))            = 2;
    
    side(ismember(pathMem, 1))      = 1;
    side(ismember(pathMem, 2))      = 2;

%     side(ismember(stackedBootStates, 9))    = 1;
%     side(ismember(stackedBootStates, 10))   = 2;
%     side(ismember(stackedBootStates, 11))   = 1;

    stateVector         = stateIdx(emissionPath);
    
    mask        = ~isnan(group) & ~isnan(side);

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
    
    MP          = get(0, 'MonitorPositions');
    fig         = figure('units', 'pixels', 'outerposition', ...
        MP(1, :));
    drawnow;
    
    axisWidth       = 0.29*1.5;
    axisHeight      = 0.225*(4/2);
    vertTol         = 0.015*(4/2);
    horzTol         = 0.04*1.5;
    
    handles         = cell(2, 2);
    for p=1:nPlots
        axesVars    = plotIndices{p};
        nAxes       = length(axesVars);
        
        for aa=1:nAxes
            handles{aa, p}      = subplot(nAxes, nPlots, ...
                sub2ind([nPlots nAxes], p, aa));
        end
    end
    
    for p=1:nPlots
        axesVars    = plotIndices{p};
        nAxes       = length(axesVars);

        for aa=1:nAxes
            hSub    = handles{aa, p};
            axes(hSub);

            hold on;

            dimIdx          = axesVars(aa);
            nLines          = length(dimIdx);
            
            violins     = sidedViolin(colvec(allSamples(mask, dimIdx, :)), ...
                group(mask), side(mask), ...
                'showData', false, 'maxWidth', 0.35, ...
                'isContinuous', isContinuous{p}(aa), ...
                'resolution', resolution{p}(aa), 'normType', 'count', ...
                'normWithinSplit', false);

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
            
            leftMargin          = (p-1)*(axisWidth + horzTol) + 0.015;
            bottomMargin        = (nAxes-aa)*(axisHeight + vertTol) + ...
                0.032;
            
            hSub.Position(1)        = leftMargin + (horzTol)/2;
            hSub.Position(2)        = bottomMargin + (vertTol)/2;
            hSub.Position(3)        = axisWidth;
            hSub.Position(4)        = axisHeight;
            
            drawnow;
%             hSub.Position(4)    = hSub.Position(4) + 0.05;
%             hSub.Position(3)    = hSub.Position(3) + 0.04;
            
%             if mod(aa, 2) == 0
%                 hSub.YAxisLocation  = 'right';
%             end
            
            hSub.YLim       = [min(obsStats(:, dimIdx)) ...
                max(obsStats(:, dimIdx))];
            hSub.XTick          = 1:max(group);
            
            if aa == nAxes
                hSub.XTickLabel     = {'8  |  9', '10  |  11', ...
                    '12A  |  12B'};
                hSub.Position(2)    = 0.04;
            else
                hSub.XTickLabel     = [];
            end
            
            dummy           = plot(0, 0, '.', 'Color', [1 1 1]);
            legend(dummy, allData.phiNames(dimIdx), 'Interpreter', 'none');
            hSub.Legend.Location    = 'best';
            
            if dimIdx == 12
                hSub.YLim(2)        = 310;
            end
            
            for vv = 1:size(violins, 1)
                violins{vv, 1}.FaceColor        = hSub.ColorOrder(end-2, :);
                violins{vv, 2}.FaceColor        = hSub.ColorOrder(2, :);
            end
            hold off;
        end
        
        linkaxes(cat(2, handles{:, p}), 'x');
        hSub.XLim       = [0.5 max(group)+0.5];
    end
end