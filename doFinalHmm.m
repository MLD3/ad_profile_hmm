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

trajGrid            = {trajectories};
seriesGrid          = {series};

initMethod          = 'score-perturbed';
nRestarts           = 100;

nSeries             = length(trajectories);
kGrid               = 8:2:30;

[kOpt, phiOpt]      = ndgrid(colvec(repmat(kGrid, nRestarts, 1)), ...
    1:length(trajGrid));
nHmmExp             = numel(kOpt);

hmmModels           = cell(size(kOpt));
loglik              = nan(nSeries, nHmmExp);
dataLoglik          = nan(nSeries, nHmmExp);

parfor exp=1:nHmmExp
    fprintf('Starting HMM experiment %d at %s.\n', exp, datestr(now));
    
    nStates         = kOpt(exp);
    traj            = trajGrid{phiOpt(exp)};
    ser             = seriesGrid{phiOpt(exp)};
    
    initModel       = getInitModel(traj, nStates, 'discrete', ...
        'method', initMethod, 'nDiscreteValues', nDiscreteValues, ...
        'score', score, 'nExpansions', 0, 'dx', dx, 'sortModel', false, ...
        'yl', yl, 'yu', yu, 'leftToRight', false, 'nPerturbGroups', 5);    
    
    initFn      = @(~, X, r)getInitModel(cellfun(@transpose, ...
        X, 'uniformoutput', false), nStates, 'discrete', 'method', ...
        initMethod, 'nDiscreteValues', nDiscreteValues, ...
        'score', score, 'nExpansions', 0, 'dx', dx, 'sortModel', false, ...
        'yl', yl, 'yu', yu, 'leftToRight', false, 'nPerturbGroups', 5);
    
    model       = hmmFitEm(ser, initModel.nstates, 'discrete', ...
        'model', initModel, 'initFn', initFn, 'verbose', false, ...
        'nRandomRestarts', 1);
    
    hmmModels{exp}          = sortHmmModel(model, 'qLims', qLims, ...
        'sortMethod', 'adas-amyloid');
    loglik(:, exp)          = hmmLogprob(model, ser);
    dataLogLik(:, exp)      = hmmCompleteLogprob(model, ser);
    
    fprintf('Done HMM experiment %d at %s!\n\n', exp, datestr(now));
end

save('profile-hmm-final');