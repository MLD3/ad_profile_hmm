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

trajectories        = mat2cell(discretePhi, counts);
series              = cellfun(@transpose, trajectories, ...
    'UniformOutput', false);

initMethod          = 'score-perturbed';
nRestarts           = 100;
nBoot               = 1000;

nSeries             = length(trajectories);
kGrid               = 12;

[~, bootSamples]    = bootstrp(nBoot, [], colvec(1:nSeries));

[kOpt, bootOpt]     = ndgrid(kGrid, 1:nBoot);
nHmmExp             = numel(kOpt);

hmmModels           = cell(size(kOpt));
hmmStates           = cell(size(kOpt));
loglik              = nan(nSeries, nHmmExp);
dataLoglik          = nan(nSeries, nHmmExp);

bootIdx             = bootSamples(:, bootOpt);

parfor exp=1:nHmmExp
    fprintf('Starting HMM experiment %d at %s.\n', exp, datestr(now));
    
    nStates         = kOpt(exp);
    
    traj            = trajectories(bootIdx(:, exp));
    ser             = series(bootIdx(:, exp));
    sc              = score(bootIdx(:, exp));
    
    initModel       = getInitModel(traj, nStates, 'discrete', ...
        'method', initMethod, 'nDiscreteValues', nDiscreteValues, ...
        'score', sc, 'nExpansions', 0, 'dx', dx, 'sortModel', false, ...
        'yl', yl, 'yu', yu, 'leftToRight', false, 'nPerturbGroups', 5);    
    
    initFn      = @(~, X, r)getInitModel(cellfun(@transpose, ...
        X, 'uniformoutput', false), nStates, 'discrete', 'method', ...
        initMethod, 'nDiscreteValues', nDiscreteValues, ...
        'score', sc, 'nExpansions', 0, 'dx', dx, 'sortModel', false, ...
        'yl', yl, 'yu', yu, 'leftToRight', false, 'nPerturbGroups', 5);
    
    model       = hmmFitEm(ser, initModel.nstates, 'discrete', ...
        'model', initModel, 'initFn', initFn, 'verbose', false, ...
        'nRandomRestarts', nRestarts);
    
    sortedModel             = sortHmmModel(model, 'sortMethod', ...
        'adas-amyloid', 'qLims', qLims);
    hmmModels{exp}          = sortedModel;
    loglik(:, exp)          = hmmLogprob(model, ser);
    dataLogLik(:, exp)      = hmmCompleteLogprob(model, ser);
    
    hmmStates{exp}          = cell2mat(cellfun(@(seq)...
        hmmMap(sortedModel, seq)', ser, 'uniformoutput', false));
    
    fprintf('Done HMM experiment %d at %s!\n\n', exp, datestr(now));
end

save('profile-hmm-bootstrap');