function cpd = condDiscreteProdEsCpdCreate(T, varargin)
%% Create a conditional discrete product distribution for the ES-HMM
% This differs from a tabularCPD in that it supports vector valued discrete
% observations, (it also supports scalar value observations). 
% These are assumed to be conditionally independent.
%
% T is of size nstates-by-nObsStates-nChildren
% so T(k,j,i) = p(i'th child = j | parent = k)
%
%% Optional inputs
% 'prior' - a struct with the the field 'alpha', which must be
% either a scalar or a matrix the same size as T. 
%%

% This file is from pmtk3.googlecode.com

[prior, parentState]    = process_options(varargin, ...
    'prior', [], 'parentState', []);
if isempty(prior)
   prior.alpha  = 1.1; % add 0.1 as pseudo counts (implicitly replicated) 
end

if isempty(parentState)
    parentState         = 1:size(T, 1);
end

[nstates, nObsStates, d]    = size(T); 
cpd             = structure(T, nstates, nObsStates, d, prior, parentState);
cpd.cpdType     = 'condDiscreteProd';
cpd.fitFn       = @condDiscreteProdEsCpdFit;
cpd.fitFnEss    = @condDiscreteProdEsCpdFitEss;
cpd.essFn       = @condDiscreteProdEsCpdComputeEss;
cpd.logPriorFn  = @(m)sum(log(m.T(:) + eps).*(m.prior.alpha(:)-1));
cpd.rndInitFn   = @rndInit;
end

function cpd = rndInit(cpd)
%% Randomly initialize
cpd.T       = normalize(rand(size(CPT.T), 2)); 
end

function cpd = condDiscreteProdEsCpdFit(cpd, Z, Y)
%% Fit  given fully observed data
% (MAP estimate with Dirichlet prior)
% Z(i) is the state of the parent Z in case i.
% Y(i, :) is the ith 1-by-d observation of the children
%%
nstates         = cpd.nstates;
nObsStates      = cpd.nObsStates; 
T               = cpd.T;

Z               = cpd.parentState(Z);

if isempty(cpd.prior)
    alpha   = 1;
else
    alpha   = cpd.prior.alpha;
end

nParents        = max(cpd.parentState);
T               = T(1:nParents, :, :);

for k = 1:nParents
    % nObsStates x d matrix
    T(k, :, :)  = reshape(histc(Y(Z==k, :), 1:nObsStates, 1), ...
        [1 nObsStates size(Y, 2)]);
end

T               = T(cpd.parentState, :, :);

T               = normalize(T + alpha - 1, 2);
cpd.T           = T;

end

function ess = condDiscreteProdEsCpdComputeEss(cpd, data, weights, B)
%% Compute the expected sufficient statistics for a condDiscreteProd CPD
% data     -  nobs-by-d
% weights  -  nobs-by-nstates; the marginal probability of the parent    
% B        -  ignored, but required by the interface, 
%             (since mixture emissions, e.g. condMixGaussTied, use it). 
%%
[nstates, nObsStates, d]    = size(cpd.T);

% counts(k, c, d) = p(x_d = c | Z = k)
counts                      = zeros(nstates, nObsStates, d);
if d < nObsStates*nstates
    for j = 1:d
        counts(:, :, j)     = weights'*bsxfun(@eq, data(:, j), 1:nObsStates);
    end
else
    for c = 1:nObsStates
        for k = 1:nstates
            counts(k, c, :)     = sum(bsxfun(@times, (data==c), ...
                weights(:, k)));
        end
    end
end

parentState     = cpd.parentState;
nParents        = max(parentState);

for ii=1:nParents
    counts(parentState == ii, :, :)     = repmat(sum(...
        counts(parentState == ii, :, :), 1), sum(parentState == ii), 1, 1);
end

ess.counts      = counts;

% this is redundant, we can simply normalize and get the same result since
% this is a discrete probability distribution that must add to 1
ess.wsum        = sum(weights, 1);
ess.wsum        = accumarray(parentState, colvec(ess.wsum), [], @sum);
ess.wsum        = rowvec(ess.wsum(parentState));

end

function cpd = condDiscreteProdEsCpdFitEss(cpd, ess)
%% Fit a condDiscreteProdCpd given the expected sufficient statistics
prior       = cpd.prior;
if isempty(prior)
    alpha   = 1;
else
    alpha   = prior.alpha;
end

cpd.T       = normalize(ess.counts + alpha-1, 2);

end