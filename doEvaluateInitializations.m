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

kGrid               = 4:2:30;
nRandomRestarts     = 20;
initGrid            = {'score-perturbed', 'random', 'convtime', 'score'};
% initGrid            = cat(1, initGrid{:});
ltrGrid             = [false, true];

nSeries             = length(trajectories);

[kOpt, initOpt, trajOpt, ltrOpt]    = ndgrid(kGrid, 1:length(initGrid), ...
    1:length(trajGrid), ltrGrid);
initOpt             = initGrid(initOpt);
% seriesOpt           = seriesGrid(trajOpt);
% trajOpt             = trajGrid(trajOpt);
nHmmExp             = numel(kOpt);

hmmModels           = cell(size(kOpt));
loglik              = nan(nSeries, nHmmExp);

parfor exp=1:nHmmExp
    fprintf('Starting HMM experiment %d at %s.\n', exp, datestr(now));
    
    nStates         = kOpt(exp);
    initMethod      = initOpt{exp};
    leftToRight     = ltrOpt(exp);
    traj            = trajGrid{trajOpt(exp)};
    ser             = seriesGrid{trajOpt(exp)};
    
    initFn          = @(~, X, r)getInitModel(cellfun(@transpose, X, ...
        'UniformOutput', false), nStates, 'discrete', 'method', ...
        initMethod, 'nDiscreteValues', nDiscreteValues, ...
        'score', score, 'nExpansions', 0, 'dx', dx, 'sortModel', false, ...
        'yl', yl, 'yu', yu, 'nPerturbGroups', 5, 'leftToRight', leftToRight);
    
    initModel       = getInitModel(traj, nStates, 'discrete', ...
        'method', initMethod, 'nDiscreteValues', nDiscreteValues, ...
        'score', score, 'nExpansions', 0, 'dx', dx, 'sortModel', false, ...
        'yl', yl, 'yu', yu, 'nPerturbGroups', 5, ...
        'leftToRight', leftToRight);
    [model, llHist]     = hmmFitEm(ser, initModel.nstates, 'discrete', ...
        'model', initModel, 'initFn', initFn, 'verbose', false, ...
        'nRandomRestarts', nRandomRestarts);
    
    hmmModels{exp}          = model;
    loglik(:, exp)          = hmmLogprob(model, ser);
    
    fprintf('Done HMM experiment %d at %s!\n\n', exp, datestr(now));
end

save('profile-hmm-initialization-comparison-2.mat');