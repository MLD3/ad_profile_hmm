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
[discretePhi, qLims]    = getProjectedFeatures(rawPhi, ...
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

trajGrid            = rowvec({trajectories});
seriesGrid          = rowvec({series});

clfSeriesGrid       = [seriesGrid; seriesGrid];
for ii=1:size(clfSeriesGrid, 2)
    for jj=1:length(clfSeriesGrid{2, ii})
        clfSeriesGrid{2, ii}{jj}(clinicalIndex, :)     = nan;
    end
end
nClfSeriesOpt       = 2;

doBaseline      = false;

if doBaseline
    sc          = scoreNan;
else
    sc          = score;
end
initMethod      = 'score-perturbed';

nSeries             = length(trajectories);
kGrid               = 4:2:30;
nRestarts           = 100;

[kOpt, phiOpt]      = ndgrid(rowvec(repmat(kGrid, nRestarts, 1)), ...
    1:length(trajGrid));
nHmmExp             = numel(kOpt);

% How many HMM experiments there
hmmModels           = cell(size(kOpt));
loglik              = nan(nSeries, nHmmExp);
dataLogLik          = nan(nSeries, nHmmExp);
loglikZ             = nan(nSeries, nHmmExp);

% Given a HMM model, we can infer the posterior over states using different
% combinations of states, i.e. use or do not use clinical features in this
% case
filteredState       = cell([nHmmExp nClfSeriesOpt]);
[instTestIdx, ~]    = find(hofInst.test);
[serTestIdx, ~]     = find(hofSeries.test);

parfor exp=1:nHmmExp
    fprintf('Starting HMM experiment %d at %s.\n', exp, datestr(now));
    
    nStates         = kOpt(exp);
    traj            = trajGrid{phiOpt(exp)};
    ser             = seriesGrid{phiOpt(exp)};
    
    pp              = cell(hofSeries.NumTestSets, nClfSeriesOpt);
    
    ll              = cell(hofSeries.NumTestSets, 1);
    llD             = cell(hofSeries.NumTestSets, 1);
    llZ             = cell(hofSeries.NumTestSets, 1);
    
    % train the HMM and infer the posterior for held-out test data
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
                'initFn', initFn, 'verbose', false, ...
                'nRandomRestarts', 1);
        else
            model           = initModel;
        end
        
        % infer the posterior over states for each option in clfSeriesGrid
        alp                 = cell(sum(testIdx), nClfSeriesOpt);
        testIndex           = find(testIdx);
        for cc=1:nClfSeriesOpt
            for t=1:sum(testIdx)
                [~, ~, alp{t, cc}]  = hmmInferNodes(model, ...
                    clfSeriesGrid{cc, phiOpt(exp)}{testIndex(t)});
            end
            pp{f, cc}   = cat(2, alp{:, cc})';
        end
        
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

    % Train the final HMM using ALL the data
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
            'nRandomRestarts', 1);
    else
        model       = initModel;
    end
    
    fs              = cell(1, nClfSeriesOpt);
    for cc=1:nClfSeriesOpt
        alpha                   = nan(nInst, nStates);
        alpha(instTestIdx, :)   = cat(1, pp{:, cc});   
        fs{cc}                  = alpha;
    end
    
    filteredState(exp, :)       = fs;
    
    hmmModels{exp}          = model;
    
    fprintf('Done HMM experiment %d at %s!\n\n', exp, datestr(now));
end

save('inter.mat');

%% Setup the classification experiments

[hmmExpOpt, durationOpt, phiOpt]    = ndgrid(1:nHmmExp, 1:nDurations, ...
    1:nClfSeriesOpt);

nClfExp             = numel(hmmExpOpt);

auroc               = cell(size(hmmExpOpt));
classifiers         = cell(size(hmmExpOpt));
clfPredictions      = cell(size(hmmExpOpt));
clfHypers           = cell(size(hmmExpOpt));

costRange           = 10.^linspace(-5, 1, 7);

for exp=1:nClfExp
    fprintf('Start classifier experiment %d at %s.\n\n', exp, datestr(now));
    
    durIdx      = durationOpt(exp);
    phiIdx      = phiOpt(exp);
    alpha       = filteredState{hmmExpOpt(exp), phiIdx};
    
    mask        = ~censored(:, durIdx);
    lab         = labels(:, durIdx);
    
    hof         = hofClf{durIdx};
    cvf         = cvfClf{durIdx};
    cvfHandle   = cvfHandleClf{durIdx};
    
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
    if mod(exp, 500) == 0
        save('inter.mat');
    end
end
% 
% %% Setup regression experiments
% 
% testError           = cell(nHmmExp, 1);
% regressors          = cell(nHmmExp, 1);
% regPredictions      = cell(nHmmExp, 1);
% regHypers           = cell(nHmmExp, 1);
% 
% for exp=1:nHmmExp
%     fprintf('Start regression experiment %d at %s.\n\n', exp, datestr(now));
%     
%     % Must exclude clinical scores for these experiments
%     alpha           = filteredState{exp, 2};
% 
%     assert(~any(colvec(isnan(alpha(regressionMask, :)))), ...
%         'Problem in filtering HMM states!');
%     
%     [~, testError{exp}, regPredictions{exp}, regressors{exp}, ...
%         regHypers{exp}]    = ...
%         hyperparameterSearch({alpha(regressionMask, :), log(target+1), ...
%         log(target+1)}, hofReg, cvfReg, @getSvcrModel, ...
%         {'classifier', 'svr', 'kernel', 'linear'}, ...
%         {'errorType', 'rmse-ci'}, ...
%         {'costL'}, {costRange}, 'searchMethod', 'exhaustive-search', ...
%         'convergenceProp', 0, 'exhaustiveArgs', ...
%         {'hyperBase', [10], 'hyperRes', [1], 'hyperPowLims', {[-5 4]}, ...
%         'budget', 120, 'nExtensions', 0}, ...
%         'outputDebug', false, 'nRandSamples', 1, 'errorFunc', ...
%         @(seq)seq.all(1), 'cvfoldHandle', cvfHandleReg, ...
%         'evalFunc', @combineSvcrError, 'cvtype', 'not-lopo', ...
%         'evalArgs', {'errorType', 'rmse-ci'}, 'hotype', hofReg.type);
%     
%     fprintf('Done experiment %d at %s!\n\n\n', exp, datestr(now));
%     if mod(exp, 200) == 0
%         save('inter.mat');
%     end
% end
% 
% % Compress regression model structures because they're quite large and
% % wasteful
% for ii=1:length(regressors)
%     regressors{ii}      = cellfun(@compressModel, regressors{ii}, ...
%         'Uniformoutput', false);
% end

save('profile-hmm-application');