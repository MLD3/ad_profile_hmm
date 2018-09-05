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
nRestarts           = 1000;
nBoot               = 1000;

nSeries             = length(trajectories);
nStates             = 14;

randModels          = cell(nRestarts, 1);
randLoglik          = nan(nSeries, nRestarts);

parfor rr=1:nRestarts
    fprintf('Starting HMM experiment %d at %s.\n', rr, datestr(now));
    
    initModel       = getInitModel(trajectories, nStates, 'discrete', ...
        'method', initMethod, 'nDiscreteValues', nDiscreteValues, ...
        'score', score, 'nExpansions', 0, 'dx', dx, 'sortModel', false, ...
        'yl', yl, 'yu', yu, 'leftToRight', false, 'nPerturbGroups', 5);    
    
    initFn      = @(~, X, r)getInitModel(cellfun(@transpose, ...
        X, 'uniformoutput', false), nStates, 'discrete', 'method', ...
        initMethod, 'nDiscreteValues', nDiscreteValues, ...
        'score', score, 'nExpansions', 0, 'dx', dx, 'sortModel', false, ...
        'yl', yl, 'yu', yu, 'leftToRight', false, 'nPerturbGroups', 5);
    
    model       = hmmFitEm(series, initModel.nstates, 'discrete', ...
        'model', initModel, 'initFn', initFn, 'verbose', false, ...
        'nRandomRestarts', 1);
    
    sortedModel             = sortHmmModel(model, 'sortMethod', ...
        'adas-amyloid', 'qLims', qLims);
    randModels{rr}          = sortedModel;
    randLoglik(:, rr)       = hmmLogprob(model, series);
    
    fprintf('Done HMM experiment %d at %s!\n', rr, datestr(now));
end

[~, bestIdx]        = max(sum(randLoglik, 1));

hmmModels           = randModels{bestIdx};
hmmStates           = cell2mat(cellfun(@(seq)...
    hmmMap(hmmModels, seq)', series, 'uniformoutput', false));

[bootModels, bootLoglik]    = getBootstrapHmmModels(...
    hmmModels, allData.counts, 'nBoot', nBoot);

save('profile-hmm-bootstrap-ci');