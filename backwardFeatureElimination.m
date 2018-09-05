function [backwardError, backwardModels, backwardHyper, droppedFeatures] = ...
    backwardFeatureElimination(phi, labels, hofolds, cvfolds, varargin)

[   minFeatures                                         , ...
    lambdaRange                                         , ...
    weights                                             , ...
    clfType                                             , ...
    cvfoldHandle                                        , ...
    phiNames            ] = process_options(varargin    , ...
    'minFeatures'       , 1                             , ...
    'lambdaRange'       , 10.^linspace(-5, 2, 15)       , ...
    'weights'           , ones(size(labels))            , ...
    'clfType'           , 'l2-logreg'                   , ...
    'cvfoldHandle'      , []                            , ...
    'phiNames'          , {}                            );

nPhi        = size(phi, 2);

currPhi             = phi;
fullPhiIdx          = 1:nPhi;

if isempty(phiNames)
    phiNames        = cellfun(@num2str, num2cell(1:nPhi), ...
        'UniformOutput', false);
end

droppedFeatures     = [];

% Train the baseline model
[~, testError, ~, models, hyper]        = elnetSearch(...
    currPhi, labels, hofolds, cvfolds, 'lambdaRange', lambdaRange, ...
    'weights', weights, 'clfType', clfType, ...
    'cvfoldHandle', cvfoldHandle, 'verbose', 0, 'nRandSamples', 1);

backwardError       = {1-testError};
backwardModels      = {models};
backwardHyper       = {hyper};

while size(currPhi, 2) > minFeatures
    currD           = size(currPhi, 2);
    phiIdx          = 1:currD;
    
    fprintf('AUROC after dropping %d feature(s): %.4f\n', nPhi-currD, ...
        max(backwardError{end}));
    
    backwardError       = cat(1, backwardError, zeros(currD, 1));
    backwardModels      = cat(1, backwardModels, {cell(currD, 1)});
    backwardHyper       = cat(1, backwardModels, {cell(currD, 1)});

    % Drop each feature one at a time
    for d=1:currD
        x           = currPhi(:, setdiff(phiIdx, d));
        [~, testError, ~, models, hyper]        = elnetSearch(...
            x, labels, hofolds, cvfolds, 'lambdaRange', lambdaRange, ...
            'weights', weights, 'clfType', clfType, ...
            'cvfoldHandle', cvfoldHandle, 'verbose', 0, ...
            'convergenceProp', 0, 'nRandSamples', 1);
        
        backwardError{end}(d)       = mean(1-testError);
        backwardModels{end}{d}      = models;
        backwardHyper{end}{d}       = hyper;
    end
    
    % Drop the feature that reduces (affects) AUROC the least
    [~, dropIdx]        = max(backwardError{end});
    
    droppedFeatures     = cat(1, droppedFeatures, fullPhiIdx(dropIdx));
    fullPhiIdx          = setdiff(fullPhiIdx, fullPhiIdx(dropIdx));
    
    currPhi             = currPhi(:, setdiff(1:currD, dropIdx));
    
    fprintf('Dropping feature %s. Num features=%d.\n', ...
        phiNames{dropIdx}, size(currPhi, 2));
    phiNames            = setdiff(phiNames, phiNames{dropIdx});
end

end