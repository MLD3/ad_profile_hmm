allPi       = cellfun(@(model)colvec(model.pi), bootModels, ...
    'UniformOutput', false);
allA        = cellfun(@(model)model.A, bootModels, 'UniformOutput', false);
allT        = cellfun(@(model)model.emission.T, bootModels, ...
    'UniformOutput', false);

allPi       = cat(2, allPi{:});
piError     = rowvec(std(allPi, 0, 2)/sqrt(nBoot));

allA        = cat(3, allA{:});
aError      = std(allA, 0, 3)/sqrt(nBoot);

allT        = cat(4, allT{:});
tError      = std(allT, 0, 4)/sqrt(nBoot);

model       = hmmModels;
modelLow    = hmmModels;
modelHi     = hmmModels;

modelLow.pi             = modelLow.pi - piError;
modelLow.A              = modelLow.A - aError;
modelLow.emission.T     = modelLow.emission.T - tError;

modelHi.pi              = modelHi.pi + piError;
modelHi.A               = modelHi.A + aError;
modelHi.emission.T      = modelHi.emission.T + tError;

expTransitions              = zeros(nStates, nStates, 3);
expTransitions(:, :, 2)     = bsxfun(@times, colvec(model.pi), model.A);
expTransitions(:, :, 1)     = bsxfun(@times, colvec(modelLow.pi), ...
    modelLow.A);
expTransitions(:, :, 3)     = bsxfun(@times, colvec(modelHi.pi), ...
    modelHi.A);
expTransitions              = expTransitions*size(allData.phi, 1);

paths                   = mat2cell(hmmStates, allData.counts);
empiricalTrans          = countTransitions(paths, nStates);

geneProb                = zeros(nStates, 3, 3);

geneProb(:, :, 2)       = [sum(model.emission.T(:, [1 2 4], 17), 2) ...
    sum(model.emission.T(:, [5 6], 17), 2) ...
    sum(model.emission.T(:, 3, 17), 2)];
geneProb(:, :, 1)       = [sum(modelLow.emission.T(:, [1 2 4], 17), 2) ...
    sum(modelLow.emission.T(:, [5 6], 17), 2) ...
    sum(modelLow.emission.T(:, 3, 17), 2)];
geneProb(:, :, 3)       = [sum(modelHi.emission.T(:, [1 2 4], 17), 2) ...
    sum(modelHi.emission.T(:, [5 6], 17), 2) ...
    sum(modelHi.emission.T(:, 3, 17), 2)];

readProb                = zeros(nStates, 3, 3);

readProb(:, :, 2)       = [sum(model.emission.T(:, 1:4, 19), 2) ...
    sum(model.emission.T(:, 5:6, 19), 2) ...
    sum(model.emission.T(:, 7, 19), 2)];
readProb(:, :, 1)       = [sum(modelLow.emission.T(:, 1:4, 19), 2) ...
    sum(modelLow.emission.T(:, 5:6, 19), 2) ...
    sum(modelLow.emission.T(:, 7, 19), 2)];
readProb(:, :, 3)       = [sum(modelHi.emission.T(:, 1:4, 19), 2) ...
    sum(modelHi.emission.T(:, 5:6, 19), 2) ...
    sum(modelHi.emission.T(:, 7, 19), 2)];

e4Count             = nan(size(allData.phi, 1), 1);
e4Count(ismember(discretePhi(:, 17), [1 2 4]))  = 0;
e4Count(ismember(discretePhi(:, 17), [5 6]))    = 1;
e4Count(ismember(discretePhi(:, 17), [3]))      = 2;

education           = discretize(rawPhi(:, 19), [0 15 18 inf]);

pathMembership      = nan(size(hmmStates));
pathMembership(ismember(hmmStates, [10 12]))      = 1;
pathMembership(ismember(hmmStates, [11 13]))        = 2;

showPathGraph       = true;
showSurvival        = false;
showTrajectories    = true;
showAge             = false;
showReading         = true;
showGenotype        = true;
showDurations       = false;
showEmission        = true;

emStates        = 9:14;

plotIndices     = {[1 6 7 10], [8 9 12 13], [11, 14, 15, 16]};

isContinuous    = {[true, true, true, true], [true, true, false, false], ...
    [true, true, true, true]};
resolution      = {[-1, -1, -1, -1], [-1 -1 3 1], [4 2 2 1]};

emissionParams  = {'emissionPath', emStates, 'quantileValues', ...
    [0.005 0.5 0.995], 'nSamples', 1000, 'plotIndices', plotIndices};

if showReading
    phiIdx          = strcmpi(allData.phiNames, 'EDUCATION');
    reading         = rawPhi(:, phiIdx);
    readingA        = reading(pathMembership == 1 & ~isnan(reading));
    readingB        = reading(pathMembership == 2 & ~isnan(reading));
    
    [~, bootA]      = bootstrp(1000, [], readingA);
    [~, bootB]      = bootstrp(1000, [], readingB);
    readingA        = readingA(bootA);
    readingB        = readingB(bootB);
    
    educationA      = discretize(readingA, [0 15 18 inf]);
    educationB      = discretize(readingB, [0 15 18 inf]);
    
    propA           = histc(educationA, 1:3)/size(educationA, 1);
    propB           = histc(educationB, 1:3)/size(educationB, 1);
    
    ciA             = quantile(propA, [0.025 0.5 0.975], 2);
    ciB             = quantile(propB, [0.025 0.5 0.975], 2);
    
    figure; hold on;
    [hBars, hErrors]    = barwitherr([ciA(:, 2) ciB(:, 2)]', ...
        cat(3, [ciA(:, 2)-ciA(:, 1) ciB(:, 2)-ciB(:, 1)]', ...
        [ciA(:, 3)-ciA(:, 2) ciB(:, 3)-ciB(:, 2)]'));
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
    
    mask                = ismember(hmmStates, 10:13) & ...
        ~isnan(education) & ismember(pathMembership, [1 2]);

    [tbl, ~, pval, lab]     = crosstab(pathMembership(mask), ...
        education(mask));
end

if showGenotype
    e4A             = e4Count(pathMembership == 1 & ~isnan(e4Count));
    e4B             = e4Count(pathMembership == 2 & ~isnan(e4Count));
    
    [~, bootA]      = bootstrp(1000, [], e4A);
    [~, bootB]      = bootstrp(1000, [], e4B);
    e4A             = e4A(bootA);
    e4B             = e4B(bootB);
    
    propA           = histc(e4A, 0:2)/size(e4A, 1);
    propB           = histc(e4B, 0:2)/size(e4B, 1);
    
    ciA             = quantile(propA, [0.025 0.5 0.975], 2);
    ciB             = quantile(propB, [0.025 0.5 0.975], 2);
    
    figure; hold on;
    [hBars, hErrors]    = barwitherr([ciA(:, 2) ciB(:, 2)]', ...
        cat(3, [ciA(:, 2)-ciA(:, 1) ciB(:, 2)-ciB(:, 1)]', ...
        [ciA(:, 3)-ciA(:, 2) ciB(:, 3)-ciB(:, 2)]'));
    ax                  = gca;
    ax.YLabel.String    = 'Probability';
    legend('\epsilon4 Non-Carrier', '\epsilon4 Single Carrier', ...
        '\epsilon4 Double Carrier');
    ax.FontSize         = 32;
    ax.YGrid            = 'on';
    hBars(1).FaceColor  = ax.ColorOrder(5, :);
    hBars(2).FaceColor  = ax.ColorOrder(1, :);
    hBars(3).FaceColor  = ax.ColorOrder(2, :);
    hBars(1).FaceAlpha  = 0.7;
    hBars(2).FaceAlpha  = 0.7;
    hBars(3).FaceAlpha  = 0.7;
    box on;
    ax.XTick            = [1 2];
    ax.XTickLabel       = {'Path A', 'Path B'};
    
    mask                = ismember(hmmStates, 10:13) & ...
        ~isnan(e4Count) & ismember(pathMembership, [1 2]);

    [tbl, ~, pval, lab]     = crosstab(pathMembership(mask), ...
        e4Count(mask));
end

% test if education and genotype are independent risk factors to be in path
% B, or are they actually correlated (they're not)

mask                = ismember(hmmStates, [11 13]) & ~isnan(e4Count) & ...
    ~isnan(education);

[tbl, ~, pval, lab]     = crosstab(education(mask), e4Count(mask));

% test for significant differences in age in stages 12 vs 13

phiIdx              = strcmpi(allData.phiNames, 'AGE');
mask                = ismember(hmmStates, [9 10 11]) & ...
    ~isnan(rawPhi(:, phiIdx));

group               = hmmStates;
group(ismember(hmmStates, [9 10]))      = 1;
group(ismember(hmmStates, 11))          = 2;

pval                = anova1(rawPhi(mask, phiIdx), group(mask), ...
    'off');

% test for significant differences in survival times, given that they start
% stages 9, 10 or 11

if showSurvival
censored            = allData.data.CONVTIME_U == 999;
convtime            = allData.data.CONVTIME_U;
convtime(censored)  = allData.data.CONVTIME_L(censored);

mask                = ismember(hmmStates, [8 9 10 11]) & (convtime >= 0);

groupMasks          = bsxfun(@and, ...
    [hmmStates == 8 hmmStates == 9 hmmStates == 10 hmmStates == 11], mask);

[f8, x8, fl8, fu8]      = ecdf(convtime(groupMasks(:, 1)), ...
    'censoring', censored(groupMasks(:, 1)), ...
    'function', 'survivor');
[f9, x9, fl9, fu9]      = ecdf(convtime(groupMasks(:, 2)), ...
    'censoring', censored(groupMasks(:, 2)), ...
    'function', 'survivor');
[f10, x10, fl10, fu10]  = ecdf(convtime(groupMasks(:, 3)), ...
    'censoring', censored(groupMasks(:, 3)), 'function', 'survivor');
[f11, x11, fl11, fu11]  = ecdf(convtime(groupMasks(:, 4)), ...
    'censoring', censored(groupMasks(:, 4)), 'function', 'survivor');
figure; hold on;
% stairs(x8, f8, 'linewidth', 1.5);
stairs(x9, f9, 'linewidth', 1.5);
stairs(x10, f10, 'linewidth', 1.5);
stairs(x11, f11, 'linewidth', 1.5);
legend('Stage 9', 'Stage 10', 'Stage 11');
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

[b, ~, H, stats]    = coxphfit(double(groupMasks(mask, [1 3])), ...
    convtime(mask), 'censoring', censored(mask));
end

%% plot the trajectories
if showTrajectories
T                       = max(allData.counts);

fullPaths       = cellfun(@(seq)[seq; ...
    ones(T-length(seq), 1)*(nStates+1)], paths, 'uniformoutput', false);
fullStates      = cell2mat(fullPaths);

timeIndex       = arrayfun(@(ii)colvec(1:size(fullPaths{ii}, 1)), ...
    1:length(trajectories), 'UniformOutput', false);

mappedStates        = fullStates;
mappedStates(ismember(fullStates, [1 2 3 4 5 6 7]))     = 1;
mappedPaths         = mat2cell(mappedStates, cellfun(@length, fullPaths));

uMappedStates       = unique(mappedStates);
nMappedStates       = length(uMappedStates);

stateTransMatrix        = zeros(nStates+1, nStates+1, T-1);

stateCountMatrix        = accumarray([mappedStates cell2mat(timeIndex')], 1);

for ii=1:length(fullPaths)
    obs         = mappedPaths{ii};
    
    stateTransMatrix    = stateTransMatrix + ...
        accumarray([obs(1:end-1) obs(2:end) colvec(1:(T-1))], 1, ...
        [nStates+1 nStates+1 T-1], @sum);
end

stateCountMatrix    = stateCountMatrix(uMappedStates, :);
stateTransMatrix    = stateTransMatrix(uMappedStates, uMappedStates, :);

subNodes        = [9 1:8];
tLim            = 10;
mask            = ~isnan(hmmStates) & ~isnan(allData.data.DX);
dxCounts        = accumarray([hmmStates(mask) allData.data.DX(mask)], 1);

dxCounts        = [sum(dxCounts(1:7, 1:2), 1); dxCounts(8:end, 2:3)];

colorScaler         = normalize(dxCounts, 2);
colorScaler         = colorScaler(:, end-1);

plotTrajectories(stateCountMatrix, stateTransMatrix, subNodes, ...
    colorScaler, 'T', tLim, 'minTrans', 0.05, 'nodeLabels', ...
    {'15', '1-7', '8', '9', '10', '11', '12', '13', '14'});
end

if showPathGraph
    subNodes        = 8:14;
    
    subGraph        = empiricalTrans(subNodes, subNodes);
    subGraph(subGraph < 20)             = 0;
    subGraph(1:length(subGraph)+1:end)  = 0;
    
    yScore          = rawPhi(:, strcmpi(allData.phiNames, 'ADAS'));
    yMask           = ~isnan(yScore);
    yLow            = accumarray(hmmStates(yMask), ...
        yScore(yMask), [], @(s)quantile(s, 0.25));
    yMedian         = accumarray(hmmStates(yMask), ...
        yScore(yMask), [], @(s)quantile(s, 0.5));
    yHigh           = accumarray(hmmStates(yMask), ...
        yScore(yMask), [], @(s)quantile(s, 0.75));
    
    xScore          = rawPhi(:, strcmpi(allData.phiNames, 'AMYLOID'));
    xMask           = ~isnan(xScore);% & stackedBootUncStatus;
    xLow            = accumarray(hmmStates(xMask), ...
        xScore(xMask), [], @(s)quantile(s, 0.25));
    xMedian         = accumarray(hmmStates(xMask), ...
        xScore(xMask), [], @(s)quantile(s, 0.5));
    xHigh           = accumarray(hmmStates(xMask), ...
        xScore(xMask), [], @(s)quantile(s, 0.75));
    
    G               = digraph(subGraph);
    
    MP              = get(0, 'MonitorPositions');
    
    fig                 = figure('units', 'pixels', 'position', MP(1, :));
    hold on;
    ax                  = gca;
    ax.Units            = 'inches';
    ax.FontSize         = 32;
    ax.Position         = [1.28 1.28 18.4 9.6];
%     ax.XLim             = [1.24 1.46];
%     ax.XTick            = 1.25:0.05:1.45;
    ax.XLabel.String        = 'Amyloid';
    ax.YLabel.String        = 'ADAS-Cog Score';
    drawnow;
    
    mask            = ~isnan(hmmStates) & ~isnan(allData.data.DX);
    dxCounts        = accumarray([hmmStates(mask) allData.data.DX(mask)], 1);
    
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
        'Marker', 'o', 'showArrows', 'off');
    
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
    
    xScaler         = ax.Position(3)/(ax.XLim(2)-ax.XLim(1));
    yScaler         = ax.Position(4)/(ax.YLim(2)-ax.YLim(1));
    
    arrowColor          = [0.8903 0.3480 0.1983];

    for ee=1:nEdges
        src         = G.Edges.EndNodes(ee, 1);
        dest        = G.Edges.EndNodes(ee, 2);
        
        dY          = graphPlot.YData(dest)-graphPlot.YData(src);
        dX          = graphPlot.XData(dest)-graphPlot.XData(src);
        slope       = (dY*yScaler)/(dX*xScaler);
        
        xPos        = mean([graphPlot.XData(src) graphPlot.XData(dest)]);
        yPos        = mean([graphPlot.YData(src) graphPlot.YData(dest)]);
        angle       = atand(slope);
        
        if src == 3 && dest == 5
            xPos        = xPos - 0.025;
            yPos        = yPos + 1.4;
        end
        
        xSrc    = graphPlot.XData(src);
        ySrc    = graphPlot.YData(src);
        
        ar                  = annotation('arrow');
        ar.Units            = 'inches';
        ar.X        = [ax.Position(1)+(xSrc-ax.XLim(1))*xScaler ...
            ax.Position(1)+(xPos-ax.XLim(1))*xScaler];
        ar.Y        = [ax.Position(2) + (ySrc-ax.YLim(1))*yScaler ...
            ax.Position(2) + ...
            (yPos-ax.YLim(1))*yScaler];
        ar.LineStyle        = 'none';
        ar.HeadWidth        = 28;
        ar.HeadLength       = 28;
        ar.Color            = arrowColor;
        
        if src == 3 && dest == 5
            xPos        = xPos - 0.01;
            yPos        = yPos + 0.8;
        end
        t           = text(xPos, yPos+0.25, num2str(G.Edges.Weight(ee)), ...
            'fontsize', 28, 'horizontalalignment', 'center', ...
            'verticalalignment', 'bottom');
        t.Rotation  = angle;
        
        edgeLabels{ee}  = t;
    end
    
    box on;
end

if showAge
    age         = rawPhi(:, strcmpi(allData.phiNames, 'AGE'));
    group       = nan(size(hmmStates));
    group(ismember(hmmStates, [9 10 11]))           = 1;
    group(ismember(hmmStates, [12 13]))             = 2;
    group(ismember(hmmStates, [14]))                = 3;
    
    side        = nan(size(hmmStates));
    
    isPathA             = cellfun(@(seq)any(ismember(seq, [9 10 12])), ...
        paths);
    isPathB             = cellfun(@(seq)any(ismember(seq, [11 13])), ...
        paths);
    both                = isPathA & isPathB;
    isPathA(both)       = false;
    isPathB(both)       = false;
    
    pathMembership      = arrayfun(@(ii)(double(isPathA(ii)) + ...
        double(isPathB(ii))*2) + zeros(size(paths{ii})), ...
        colvec(1:length(isPathA)), 'UniformOutput', false);
    pathMembership      = cell2mat(pathMembership);
    
    side(pathMembership == 1)       = 1;
    side(pathMembership == 2)       = 2;
    
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
    ax.XTickLabel       = {'9/10  |   11  ', ...
        '12  |  13', '14'};
    box on;
    
    for ii=1:size(ageViolins, 1)
        ageViolins{ii, 1}.FaceColor         = ax.ColorOrder(end-2, :);
        ageViolins{ii, 2}.FaceColor         = ax.ColorOrder(2, :);
    end
    
    legend([ageViolins{ii, 1} ageViolins{ii, 2}], {'Path A', 'Path B'});
end

if showDurations
    startStates     = [9 10 11];
    nPaths          = length(startStates);
    
    nSamples        = 200;
    T               = 50;
    pathTimes       = cell(nBoot, nPaths);
    
    parfor bb=1:nBoot
        model       = bootModels{bb};
        for pp=1:nPaths
            pathTimes{bb, pp}   = getHmmStepCounts(model, startStates(pp), ...
                nStates, 'nSamples', nSamples, 'sampleTillReached', false, ...
                'T', T);
        end
    end
    
    sampledTimes        = cell(nPaths, 1);
    
    for pp=1:nPaths
        sampledTimes{pp}    = cell2mat(pathTimes(:, pp));
    end
    
    allTimes        = cat(2, sampledTimes{:});
    groups          = repmat(1:nPaths, nSamples*nBoot, 1);
    bootIdx         = repmat(colvec(repmat(1:nBoot, nSamples, 1)), 1, nPaths);
    
    frequencies     = accumarray([colvec(allTimes) colvec(groups) ...
        colvec(bootIdx)], 1);
    
    probStats       = quantile(cumsum(frequencies/nSamples, 1), ...
        [0.25 0.5 0.75], 3);
    
    tLim            = 25;
    
    jitter          = repmat(linspace(-0.06, 0.06, 3), tLim, 1);
    
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
    legend('Stage 9 -> Terminal Stage', 'Stage 10 -> Terminal Stage', ...
        'Stage 11 -> Terminal Stage');
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
    
    d               = size(rawTrajectories{1}, 2);
    
    allSamples      = rawPhi;
    zSample         = hmmStates;

%     [~, allSamples]     = getSampleQuantiles(model.emission, ...
%         discretePhi, rawPhi, 1:14, hmmStates, 'sampleObs', true);
%     zSample             = colvec(repmat(1:14, size(allSamples, 1), 1));
%     
%     allSamples      = reshape(allSamples, [], d);
%     allSamples(allSamples(:, 1) > 1.2e-8, 1)    = nan;
        
    obsStats            = quantile(allSamples, quantileValues, 1);

    side                = nan(size(zSample));
    group               = nan(size(zSample));
%     
    group(ismember(zSample, [10, 11]))      = 1;
    group(ismember(zSample, [12, 13]))      = 2;
    group(ismember(zSample, 14))            = 3;

%     group(zSample == 9)         = 1;
%     group(zSample == 10)        = 1;

    isPathA             = cellfun(@(seq)any(ismember(seq, [10 12])), ...
        paths);
    isPathB             = cellfun(@(seq)any(ismember(seq, [11 13])), ...
        paths);
    both                = isPathA & isPathB;
    isPathA(both)       = false;
    isPathB(both)       = false;
    
    pathMembership      = arrayfun(@(ii)(double(isPathA(ii)) + ...
        double(isPathB(ii))*2) + zeros(size(paths{ii})), ...
        colvec(1:length(isPathA)), 'UniformOutput', false);
    pathMembership      = cell2mat(pathMembership);
    
    side(pathMembership == 1)       = 1;
    side(pathMembership == 2)       = 2;
%     side(ismember(zSample, [10 12 14]))     = 1;
%     side(ismember(zSample, [11 13]))        = 2;
    
%     side                    = ones(size(side));
%     side(zSample == 9)      = 1;
%     side(zSample == 10)     = 2;

    stateIdx            = colvec(1:nStates);
    stateVector         = stateIdx(emissionPath);
    
    mask        = ~isnan(colvec(group)) & ~isnan(colvec(side));

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
            
            dimIdx          = axes(aa);
            nLines          = length(dimIdx);
            
            violins         = sidedViolin(...
                colvec(allSamples(mask, dimIdx, :)), ...
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
            
            if aa < nAxes
                hSub.XTickLabel     = [];
            end
            
            hSub.Position(4)    = hSub.Position(4) + 0.05;
            hSub.Position(3)    = hSub.Position(3) + 0.04;
            
            if mod(aa, 2) == 0
                hSub.YAxisLocation  = 'right';
            end
            
            hSub.YLim       = [min(obsStats(:, dimIdx)) ...
                max(obsStats(:, dimIdx))];
            hSub.XTick          = 1:max(group);
%             hSub.XTickLabel     = stateVector;
            hSub.XTickLabel     = {'10  |  11', '12  |  13', '14A  |  14B'};
            
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
        end
        
        linkaxes(cat(2, subPlots{:}), 'x');
        hSub.XLim       = [0.5 max(group)+0.5];
    end
end

startIdx        = cell(size(paths));
endIdx          = cell(size(paths));

for ii=1:length(paths)
    obs         = paths{ii};
    
    startIdx{ii}        = nan(size(obs));
    endIdx{ii}          = nan(size(obs));
    
    if all(ismember([9 12], obs))
        startIdx{ii}(find(obs == 9, 1, 'last'))     = 1;
        endIdx{ii}(find(obs == 9, 1, 'last'))       = 1;
    end
end

startIdx        = cell2mat(startIdx);
endIdx          = cell2mat(endIdx);



