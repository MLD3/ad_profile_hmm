function [bootModels, bootLoglik] = getBootstrapHmmModels(model, counts, ...
    varargin)

[   nBoot               ] = process_options(varargin, ...
    'nBoot'             , 1000                      );

bootModels      = cell(nBoot, 1);
bootLoglik      = zeros(nBoot, 1);

parfor bb=1:nBoot
    seriesSample    = arrayfun(@(seq)generateDataSample(model, seq), ...
        counts(randperm(length(counts))), 'uniformoutput', false);
    
    bootModels{bb}  = hmmFitEm(seriesSample, model.nstates, 'discrete', ...
        'model', model, 'initFn', [], 'verbose', false, ...
        'nRandomRestarts', 1);
    bootLoglik(bb)  = sum(hmmLogprob(bootModels{bb}, seriesSample));
end

end

function obsSeq = generateDataSample(model, n)

stateSeq        = [];
obsSeq          = [];

while length(stateSeq) < n
    if isempty(stateSeq)
        nextState       = find(mnrnd(1, model.pi));
    else
        nextState       = find(mnrnd(1, model.A(stateSeq(end), :)));
    end
    
    stateSeq        = cat(2, stateSeq, nextState);
    
    obsDraws        = mnrnd(1, squeeze(model.emission.T(nextState, :, :))');
    [obsVals, ~]    = find(obsDraws');

    obsSeq          = cat(2, obsSeq, obsVals);
end

end