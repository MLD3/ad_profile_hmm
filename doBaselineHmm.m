load profile-hmm-data.mat;

%% pick out the clinical score used to intialize the data; we use ADAS-Cog
scoreIdx            = strcmpi('ADAS', allData.phiNames);
score               = cellfun(@(seq)seq(:, scoreIdx), trajectories, ...
    'UniformOutput', false);
scoreNan            = score;

meanScore           = nanmean(cat(1, score{:}));

for ii=1:length(score)
    patScore        = score{ii};
    patFollow       = followup{ii};
    
    obsScore        = ~isnan(patScore);
    
    if sum(obsScore) >= 2
        score{ii}   = interp1(patFollow(obsScore), patScore(obsScore), ...
            patFollow, 'linear', 'extrap');
    elseif sum(obsScore) == 1
        score{ii}   = ones(size(obsScore))*patScore(obsScore);
    else
        score{ii}   = ones(size(obsScore))*meanScore;
    end
end

%% Discretize the data
rawTrajectories     = trajectories;

rawPhi              = cat(1, rawTrajectories{:});
counts              = cellfun(@(seq)size(seq, 1), rawTrajectories);

nDiscreteValues     = 10;
skipDiscrete        = false(length(allData.phiNames), 1);
skipDiscrete(strcmpi(allData.phiNames, 'ApoE-Profile'))     = true;

%% Assume missingness at random (MAR) and marginalize missed observations
discretePhi         = getProjectedFeatures(rawPhi, ...
    'quantile', 'nQuantiles', nDiscreteValues, ...
    'binaryOutput', false, 'isMissing', isnan(rawPhi), ...
    'removeRedundant', false, 'missingInd', false, ...
    'skipDiscrete', skipDiscrete);

geneIdx             = strcmpi(allData.phiNames, 'ApoE-Profile');
demogIdx            = strcmpi(allData.phiNames, 'AGE') | ...
    strcmpi(allData.phiNames, 'EDUCATION');
clinicalIndex       = find(allData.group == 5);

trajectories        = mat2cell(discretePhi, counts);
series              = cellfun(@transpose, trajectories, ...
    'UniformOutput', false);

trajNoGenetic       = mat2cell(discretePhi(:, ~geneIdx), counts);
seriesNoGenetic     = cellfun(@transpose, trajNoGenetic, ...
    'UniformOutput', false);
trajNoDemog         = mat2cell(discretePhi(:, ~demogIdx), counts);
seriesNoDemog       = cellfun(@transpose, trajNoDemog, ...
    'UniformOutput', false);
trajNoBoth          = mat2cell(discretePhi(:, (~demogIdx) & (~geneIdx)), ...
    counts);
seriesNoBoth        = cellfun(@transpose, trajNoBoth, ...
    'UniformOutput', false);

trajGrid            = {trajectories, trajNoGenetic, trajNoDemog, ...
    trajNoBoth};
seriesGrid          = {series, seriesNoGenetic, seriesNoDemog, ...
    seriesNoBoth};

doBaseline      = true;

if doBaseline
    sc          = scoreNan;
else
    sc          = score;
end
initMethod      = 'convtime';

nSeries             = length(trajectories);
kGrid               = 4:2:30;

[kOpt, phiOpt]      = ndgrid(kGrid, 1:length(trajGrid));
nHmmExp             = numel(kOpt);

hmmModels           = cell(size(kOpt));
loglik              = nan(nSeries, nHmmExp);
dataLogLik          = nan(nSeries, nHmmExp);
loglikZ             = nan(nSeries, nHmmExp);

filteredState       = cell(size(kOpt));
[instTestIdx, ~]    = find(hofInst.test);
[serTestIdx, ~]     = find(hofSeries.test);

parfor exp=1:nHmmExp
    fprintf('Starting HMM experiment %d at %s.\n', exp, datestr(now));
    
    nStates         = kOpt(exp);
    traj            = trajGrid{phiOpt(exp)};
    ser             = seriesGrid{phiOpt(exp)};
    
    serNoClinical   = ser;
    for ii=1:length(serNoClinical)
        serNoClinical{ii}(clinicalIndex, :)     = nan;
    end
    
    alpha           = nan(nInst, nStates);
    pp              = cell(hofSeries.NumTestSets, 1);
    
    ll              = cell(hofSeries.NumTestSets, 1);
    llD             = cell(hofSeries.NumTestSets, 1);
    llZ             = cell(hofSeries.NumTestSets, 1);
    
    for f=1:hofSeries.NumTestSets
        trainIdx    = hofSeries.training(:, f);
        testIdx     = hofSeries.test(:, f);
        
        initModel       = getInitModel(traj(trainIdx), ...
            nStates, 'discrete', 'method', initMethod, ...
            'nDiscreteValues', nDiscreteValues, 'score', ...
            sc(trainIdx), 'nExpansions', 0, 'dx', dx(trainIdx), ...
            'sortModel', false, 'yl', yl(trainIdx), 'yu', yu(trainIdx), ...
            'leftToRight', false, 'nPerturbGroups', 5);
        
        if ~doBaseline
            initFn      = @(~, X, r)getInitModel(cellfun(@transpose, ...
                X, 'uniformoutput', false), nStates, 'discrete', ...
                'method', initMethod, 'nDiscreteValues', ...
                nDiscreteValues, 'score', sc(trainIdx), 'nExpansions', ...
                0, 'dx', dx(trainIdx), 'sortModel', false, 'yl', ...
                yl(trainIdx), 'yu', yu(trainIdx), 'leftToRight', false, ...
                'nPerturbGroups', 5);
            [model, ~]      = hmmFitEm(ser(trainIdx), ...
                initModel.nstates, 'discrete', 'model', initModel, ...
                'initFn', initFn, 'verbose', false, 'nRandomRestarts', 20);
        else
            model           = initModel;
        end
        
        alp                 = cell(sum(testIdx), 1);
        testIndex           = find(testIdx);
        for t=1:sum(testIdx)
            [~, ~, alp{t}]  = hmmInferNodes(model, ...
                serNoClinical{testIndex(t)});
        end
        pp{f}           = cat(2, alp{:})';
        
        trainLL         = hmmLogprob(model, ser(trainIdx));
        testLL          = hmmLogprob(model, ser(testIdx));
        dataLL          = hmmCompleteLogprob(model, ser(testIdx));
        testLLZ         = getLoglikZscore(trainLL, testLL, ...
            trainIdx, testIdx, counts);
        
        ll{f}           = testLL;
        llZ{f}          = testLLZ;
        llD{f}          = dataLL;
    end
    
    loglik(serTestIdx, exp)         = cat(1, ll{:});
    loglikZ(serTestIdx, exp)        = cat(1, llZ{:});
    dataLogLik(serTestIdx, exp)     = cat(1, llD{:});
    
    initModel       = getInitModel(traj, nStates, 'discrete', ...
        'method', initMethod, 'nDiscreteValues', nDiscreteValues, ...
        'score', sc, 'nExpansions', 0, 'dx', dx, 'sortModel', false, ...
        'yl', yl, 'yu', yu, 'leftToRight', false, 'nPerturbGroups', 5);
    
    if ~doBaseline
        initFn      = @(~, X, r)getInitModel(cellfun(@transpose, ...
            X, 'uniformoutput', false), nStates, 'discrete', 'method', ...
            initMethod, 'nDiscreteValues', nDiscreteValues, ...
            'score', sc, 'nExpansions', 0, 'dx', dx, 'sortModel', false, ...
            'yl', yl, 'yu', yu, 'leftToRight', false, 'nPerturbGroups', 5);
        
        model       = hmmFitEm(ser, initModel.nstates, 'discrete', ...
            'model', initModel, 'initFn', initFn, 'verbose', false, ...
            'nRandomRestarts', 20);
    else
        model       = initModel;
    end
    
    alpha(instTestIdx, :)   = cat(1, pp{:});   
    filteredState{exp}      = alpha;
    hmmModels{exp}          = model;
    
    fprintf('Done HMM experiment %d at %s!\n\n', exp, datestr(now));
end

%% Setup the classification experiments

[hmmExpIdx, durationIdx]        = ndgrid(1:nHmmExp, 1:nDurations);

nClfExp             = numel(hmmExpIdx);

auroc               = cell(size(hmmExpIdx));
classifiers         = cell(size(hmmExpIdx));
clfPredictions      = cell(size(hmmExpIdx));
clfHypers           = cell(size(hmmExpIdx));

costRange           = 10.^linspace(-5, 1, 7);

for exp=1:nClfExp
    fprintf('Start classifier experiment %d at %s.\n\n', exp, datestr(now));
    
    durOpt      = durationIdx(exp);
    alpha       = filteredState{hmmExpIdx(exp)};
    
    mask        = ~censored(:, durOpt);
    lab         = labels(:, durOpt);
    
    hof         = hofClf{durOpt};
    cvf         = cvfClf{durOpt};
    cvfHandle   = cvfHandleClf{durOpt};
    
    assert(all(ismember(unique(lab(mask)), [-1 1])), ...
        'Problem in filtering out censored patients!');
    assert(~any(colvec(isnan(alpha(mask, :)))), ...
        'Problem in filtering HMM states!');
    
    [~, auroc{exp}, clfPredictions{exp}, classifiers{exp}, ...
        clfHypers{exp}]    = ...
        hyperparameterSearch({alpha(mask, :), lab(mask)}, hof, cvf, ...
        @getSvmModel, ...
        {'classifier', 'logreg', 'kernel', 'linear', ...
        'weightMethod', 'label'}, {'metric', 'auroc'}, ...
        {'cost'}, {costRange}, 'searchMethod', 'exhaustive-search', ...
        'convergenceProp', 0, 'exhaustiveArgs', ...
        {'hyperBase', [10], 'hyperRes', [1], 'hyperPowLims', {[-5 4]}, ...
        'budget', 120, 'nExtensions', 0}, ...
        'outputDebug', false, 'nRandSamples', 1, 'errorFunc', ...
        @(seq)(1-seq), 'cvfoldHandle', cvfHandle, ...
        'evalFunc', @combineClassifierError, 'cvtype', 'not-lopo', ...
        'hotype', hof.type);
    
    fprintf('Done classifier experiment %d at %s!\n\n\n', exp, datestr(now));
    if mod(exp, 60) == 0
        save('inter.mat');
    end
end

%% Setup regression experiments

testError           = cell(nHmmExp, 1);
regressors          = cell(nHmmExp, 1);
regPredictions      = cell(nHmmExp, 1);
regHypers           = cell(nHmmExp, 1);

for exp=1:nHmmExp
    fprintf('Start regression experiment %d at %s.\n\n', exp, datestr(now));
    
    alpha           = filteredState{exp};

    assert(~any(colvec(isnan(alpha(regressionMask, :)))), ...
        'Problem in filtering HMM states!');
    
    [~, testError{exp}, regPredictions{exp}, regressors{exp}, ...
        regHypers{exp}]    = ...
        hyperparameterSearch({alpha(regressionMask, :), target, ...
        target}, hofReg, cvfReg, @getSvcrModel, ...
        {'classifier', 'svr', 'kernel', 'linear'}, ...
        {'errorType', 'rmse-ci'}, ...
        {'costL'}, {costRange}, 'searchMethod', 'exhaustive-search', ...
        'convergenceProp', 0, 'exhaustiveArgs', ...
        {'hyperBase', [10], 'hyperRes', [1], 'hyperPowLims', {[-5 4]}, ...
        'budget', 120, 'nExtensions', 0}, ...
        'outputDebug', false, 'nRandSamples', 1, 'errorFunc', ...
        @(seq)seq.all(1), 'cvfoldHandle', cvfHandleReg, ...
        'evalFunc', @combineSvcrError, 'cvtype', 'not-lopo', ...
        'evalArgs', {'errorType', 'rmse-ci'}, 'hotype', hofReg.type);
    
    fprintf('Done experiment %d at %s!\n\n\n', exp, datestr(now));
    if mod(exp, 10) == 0
        save('inter.mat');
    end
end

for ii=1:length(regressors)
    regressors{ii}      = cellfun(@compressModel, regressors{ii}, ...
        'Uniformoutput', false);
end

save('profile-hmm-convtime-application-no-clinical');

clear; clc;
load profile-hmm-data.mat;

%% pick out the clinical score used to intialize the data; we use ADAS-Cog
scoreIdx            = strcmpi('ADAS', allData.phiNames);
score               = cellfun(@(seq)seq(:, scoreIdx), trajectories, ...
    'UniformOutput', false);
scoreNan            = score;

meanScore           = nanmean(cat(1, score{:}));

for ii=1:length(score)
    patScore        = score{ii};
    patFollow       = followup{ii};
    
    obsScore        = ~isnan(patScore);
    
    if sum(obsScore) >= 2
        score{ii}   = interp1(patFollow(obsScore), patScore(obsScore), ...
            patFollow, 'linear', 'extrap');
    elseif sum(obsScore) == 1
        score{ii}   = ones(size(obsScore))*patScore(obsScore);
    else
        score{ii}   = ones(size(obsScore))*meanScore;
    end
end

%% Discretize the data
rawTrajectories     = trajectories;

rawPhi              = cat(1, rawTrajectories{:});
counts              = cellfun(@(seq)size(seq, 1), rawTrajectories);

nDiscreteValues     = 10;
skipDiscrete        = false(length(allData.phiNames), 1);
skipDiscrete(strcmpi(allData.phiNames, 'ApoE-Profile'))     = true;

%% Assume missingness at random (MAR) and marginalize missed observations
discretePhi         = getProjectedFeatures(rawPhi, ...
    'quantile', 'nQuantiles', nDiscreteValues, ...
    'binaryOutput', false, 'isMissing', isnan(rawPhi), ...
    'removeRedundant', false, 'missingInd', false, ...
    'skipDiscrete', skipDiscrete);

geneIdx             = strcmpi(allData.phiNames, 'ApoE-Profile');
demogIdx            = strcmpi(allData.phiNames, 'AGE') | ...
    strcmpi(allData.phiNames, 'EDUCATION');
clinicalIndex       = find(allData.group == 5);

trajectories        = mat2cell(discretePhi, counts);
series              = cellfun(@transpose, trajectories, ...
    'UniformOutput', false);

trajNoGenetic       = mat2cell(discretePhi(:, ~geneIdx), counts);
seriesNoGenetic     = cellfun(@transpose, trajNoGenetic, ...
    'UniformOutput', false);
trajNoDemog         = mat2cell(discretePhi(:, ~demogIdx), counts);
seriesNoDemog       = cellfun(@transpose, trajNoDemog, ...
    'UniformOutput', false);
trajNoBoth          = mat2cell(discretePhi(:, (~demogIdx) & (~geneIdx)), ...
    counts);
seriesNoBoth        = cellfun(@transpose, trajNoBoth, ...
    'UniformOutput', false);

trajGrid            = {trajectories, trajNoGenetic, trajNoDemog, ...
    trajNoBoth};
seriesGrid          = {series, seriesNoGenetic, seriesNoDemog, ...
    seriesNoBoth};

doBaseline      = true;

if doBaseline
    sc          = scoreNan;
else
    sc          = score;
end
initMethod      = 'score';

nSeries             = length(trajectories);
kGrid               = 4:2:30;

[kOpt, phiOpt]      = ndgrid(kGrid, 1:length(trajGrid));
nHmmExp             = numel(kOpt);

hmmModels           = cell(size(kOpt));
loglik              = nan(nSeries, nHmmExp);
dataLogLik          = nan(nSeries, nHmmExp);
loglikZ             = nan(nSeries, nHmmExp);

filteredState       = cell(size(kOpt));
[instTestIdx, ~]    = find(hofInst.test);
[serTestIdx, ~]     = find(hofSeries.test);

parfor exp=1:nHmmExp
    fprintf('Starting HMM experiment %d at %s.\n', exp, datestr(now));
    
    nStates         = kOpt(exp);
    traj            = trajGrid{phiOpt(exp)};
    ser             = seriesGrid{phiOpt(exp)};
    
    serNoClinical   = ser;
    for ii=1:length(serNoClinical)
        serNoClinical{ii}(clinicalIndex, :)     = nan;
    end
    
    alpha           = nan(nInst, nStates);
    pp              = cell(hofSeries.NumTestSets, 1);
    
    ll              = cell(hofSeries.NumTestSets, 1);
    llD             = cell(hofSeries.NumTestSets, 1);
    llZ             = cell(hofSeries.NumTestSets, 1);
    
    for f=1:hofSeries.NumTestSets
        trainIdx    = hofSeries.training(:, f);
        testIdx     = hofSeries.test(:, f);
        
        initModel       = getInitModel(traj(trainIdx), ...
            nStates, 'discrete', 'method', initMethod, ...
            'nDiscreteValues', nDiscreteValues, 'score', ...
            sc(trainIdx), 'nExpansions', 0, 'dx', dx(trainIdx), ...
            'sortModel', false, 'yl', yl(trainIdx), 'yu', yu(trainIdx), ...
            'leftToRight', false, 'nPerturbGroups', 5);
        
        if ~doBaseline
            initFn      = @(~, X, r)getInitModel(cellfun(@transpose, ...
                X, 'uniformoutput', false), nStates, 'discrete', ...
                'method', initMethod, 'nDiscreteValues', ...
                nDiscreteValues, 'score', sc(trainIdx), 'nExpansions', ...
                0, 'dx', dx(trainIdx), 'sortModel', false, 'yl', ...
                yl(trainIdx), 'yu', yu(trainIdx), 'leftToRight', false, ...
                'nPerturbGroups', 5);
            [model, ~]      = hmmFitEm(ser(trainIdx), ...
                initModel.nstates, 'discrete', 'model', initModel, ...
                'initFn', initFn, 'verbose', false, 'nRandomRestarts', 20);
        else
            model           = initModel;
        end
        
        alp                 = cell(sum(testIdx), 1);
        testIndex           = find(testIdx);
        for t=1:sum(testIdx)
            [~, ~, alp{t}]  = hmmInferNodes(model, ...
                serNoClinical{testIndex(t)});
        end
        pp{f}           = cat(2, alp{:})';
        
        trainLL         = hmmLogprob(model, ser(trainIdx));
        testLL          = hmmLogprob(model, ser(testIdx));
        dataLL          = hmmCompleteLogprob(model, ser(testIdx));
        testLLZ         = getLoglikZscore(trainLL, testLL, ...
            trainIdx, testIdx, counts);
        
        ll{f}           = testLL;
        llZ{f}          = testLLZ;
        llD{f}          = dataLL;
    end
    
    loglik(serTestIdx, exp)         = cat(1, ll{:});
    loglikZ(serTestIdx, exp)        = cat(1, llZ{:});
    dataLogLik(serTestIdx, exp)     = cat(1, llD{:});
    
    initModel       = getInitModel(traj, nStates, 'discrete', ...
        'method', initMethod, 'nDiscreteValues', nDiscreteValues, ...
        'score', sc, 'nExpansions', 0, 'dx', dx, 'sortModel', false, ...
        'yl', yl, 'yu', yu, 'leftToRight', false, 'nPerturbGroups', 5);
    
    if ~doBaseline
        initFn      = @(~, X, r)getInitModel(cellfun(@transpose, ...
            X, 'uniformoutput', false), nStates, 'discrete', 'method', ...
            initMethod, 'nDiscreteValues', nDiscreteValues, ...
            'score', sc, 'nExpansions', 0, 'dx', dx, 'sortModel', false, ...
            'yl', yl, 'yu', yu, 'leftToRight', false, 'nPerturbGroups', 5);
        
        model       = hmmFitEm(ser, initModel.nstates, 'discrete', ...
            'model', initModel, 'initFn', initFn, 'verbose', false, ...
            'nRandomRestarts', 20);
    else
        model       = initModel;
    end
    
    alpha(instTestIdx, :)   = cat(1, pp{:});   
    filteredState{exp}      = alpha;
    hmmModels{exp}          = model;
    
    fprintf('Done HMM experiment %d at %s!\n\n', exp, datestr(now));
end

%% Setup the classification experiments

[hmmExpIdx, durationIdx]        = ndgrid(1:nHmmExp, 1:nDurations);

nClfExp             = numel(hmmExpIdx);

auroc               = cell(size(hmmExpIdx));
classifiers         = cell(size(hmmExpIdx));
clfPredictions      = cell(size(hmmExpIdx));
clfHypers           = cell(size(hmmExpIdx));

costRange           = 10.^linspace(-5, 1, 7);

for exp=1:nClfExp
    fprintf('Start classifier experiment %d at %s.\n\n', exp, datestr(now));
    
    durOpt      = durationIdx(exp);
    alpha       = filteredState{hmmExpIdx(exp)};
    
    mask        = ~censored(:, durOpt);
    lab         = labels(:, durOpt);
    
    hof         = hofClf{durOpt};
    cvf         = cvfClf{durOpt};
    cvfHandle   = cvfHandleClf{durOpt};
    
    assert(all(ismember(unique(lab(mask)), [-1 1])), ...
        'Problem in filtering out censored patients!');
    assert(~any(colvec(isnan(alpha(mask, :)))), ...
        'Problem in filtering HMM states!');
    
    [~, auroc{exp}, clfPredictions{exp}, classifiers{exp}, ...
        clfHypers{exp}]    = ...
        hyperparameterSearch({alpha(mask, :), lab(mask)}, hof, cvf, ...
        @getSvmModel, ...
        {'classifier', 'logreg', 'kernel', 'linear', ...
        'weightMethod', 'label'}, {'metric', 'auroc'}, ...
        {'cost'}, {costRange}, 'searchMethod', 'exhaustive-search', ...
        'convergenceProp', 0, 'exhaustiveArgs', ...
        {'hyperBase', [10], 'hyperRes', [1], 'hyperPowLims', {[-5 4]}, ...
        'budget', 120, 'nExtensions', 0}, ...
        'outputDebug', false, 'nRandSamples', 1, 'errorFunc', ...
        @(seq)(1-seq), 'cvfoldHandle', cvfHandle, ...
        'evalFunc', @combineClassifierError, 'cvtype', 'not-lopo', ...
        'hotype', hof.type);
    
    fprintf('Done classifier experiment %d at %s!\n\n\n', exp, datestr(now));
    if mod(exp, 60) == 0
        save('inter.mat');
    end
end

%% Setup regression experiments

testError           = cell(nHmmExp, 1);
regressors          = cell(nHmmExp, 1);
regPredictions      = cell(nHmmExp, 1);
regHypers           = cell(nHmmExp, 1);

for exp=1:nHmmExp
    fprintf('Start regression experiment %d at %s.\n\n', exp, datestr(now));
    
    alpha           = filteredState{exp};

    assert(~any(colvec(isnan(alpha(regressionMask, :)))), ...
        'Problem in filtering HMM states!');
    
    [~, testError{exp}, regPredictions{exp}, regressors{exp}, ...
        regHypers{exp}]    = ...
        hyperparameterSearch({alpha(regressionMask, :), target, ...
        target}, hofReg, cvfReg, @getSvcrModel, ...
        {'classifier', 'svr', 'kernel', 'linear'}, ...
        {'errorType', 'rmse-ci'}, ...
        {'costL'}, {costRange}, 'searchMethod', 'exhaustive-search', ...
        'convergenceProp', 0, 'exhaustiveArgs', ...
        {'hyperBase', [10], 'hyperRes', [1], 'hyperPowLims', {[-5 4]}, ...
        'budget', 120, 'nExtensions', 0}, ...
        'outputDebug', false, 'nRandSamples', 1, 'errorFunc', ...
        @(seq)seq.all(1), 'cvfoldHandle', cvfHandleReg, ...
        'evalFunc', @combineSvcrError, 'cvtype', 'not-lopo', ...
        'evalArgs', {'errorType', 'rmse-ci'}, 'hotype', hofReg.type);
    
    fprintf('Done experiment %d at %s!\n\n\n', exp, datestr(now));
    if mod(exp, 10) == 0
        save('inter.mat');
    end
end

for ii=1:length(regressors)
    regressors{ii}      = cellfun(@compressModel, regressors{ii}, ...
        'Uniformoutput', false);
end

save('profile-hmm-score-application-no-clinical');