function [logp, pZ] = convtimeModelLogProb(model, data, yl, yu)
% logp(i) = log p(X{i} | model), X{i} is 1*T
% if X is a single sequence, we compute logp = log p( {X} | model)

% This file is from pmtk3.googlecode.com

data    = cellwrap(data);

assert(~isempty(yl) && ~isempty(yu), ...
    'Must provide convtime for convtime-model!');
    
stackedYl       = cat(1, yl{:});
stackedYu       = cat(1, yu{:});
stackedPhi      = cat(1, data{:});
counts          = cellfun(@length, yl);
[n, d]          = size(stackedPhi);

censored        = (stackedYl == -999) | (stackedYu == 999);
    
meanConvtime    = mean([stackedYl stackedYu], 2);
convtimeBins    = model.bins;
nStates         = model.nstates;
cpd             = model.cpd;

% flip the states; a high state number should mean a terminal disease
% stage, but not flipping it would assign a high state number of
% someone that's actually furthest away from disease!
stackedStates   = nStates - discretize(meanConvtime, convtimeBins) + 1;
logP            = zeros(n, 1);
pZ              = zeros(n, nStates);

softev          = mkProfileSoftEvidence(cpd, stackedPhi')';

% hard assign state dist. probability to uncensored instances
index           = sub2ind(size(softev), find(~censored), ...
    stackedStates(~censored));
logP(~censored)     = softev(index);
pZ(index)           = 1;

[pzM, llMarg]   = mixDiscreteInferLatent(model, stackedPhi);
logP(censored)      = llMarg(censored);
pZ(censored, :)     = pzM(censored, :);

ll          = mat2cell(logP, counts);

logp        = cellfun(@sum, ll);

end