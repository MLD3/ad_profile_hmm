function CPD = condMixedCpdCreate(mu, Sigma, T, distType, varargin)
%% Create a mixed conditional distribution for use in a graphical model
% Some dimensions of this distribution are normally distribution and others
% are discrete. The Gaussian components are jointly Gaussian, the discrete
% components are independent though.
% distType      -- d dimensional vector that is 'gaussian' or 'discrete'
% 
% mu is a matrix of size dn-by-nstates
% Sigma is of size dn-by-dn-by-nstates
% T is a matrix of nstates-by-nobsstates-dd
% dn = number of gaussian dimensions
% dd = number of discrete dimensions
%
% 'gaussPrior' is a Gauss-inverseWishart distribution, namely, a struct with
% fields  mu, Sigma, dof, k
% Set 'prior' to 'none' to do mle. 
% 'discretePrior' is a struct with the field alpha which must either be a
% scalar or a matrix the same size as T
%%

% This file is from pmtk3.googlecode.com

[gaussPrior, discretePrior]     = process_options(varargin, ...
    'gaussPrior', [], 'discretePrior', []);

dGauss          = sum(strcmpi(distType, 'gauss'));
dDiscrete       = sum(strcmpi(distType, 'discrete'));
d               = dGauss + dDiscrete;

if isempty(gaussPrior) 
      gaussPrior.mu         = zeros(1, d);
      gaussPrior.Sigma      = 0.1*eye(d);
      gaussPrior.k          = 0.01;
      gaussPrior.dof        = d + 1; 
end

if isempty(discretePrior)
    % add 0.1 as pseudo counts (implicitly replicated) 
    discretePrior.alpha     = 1.1; 
end

if isvector(Sigma)
   Sigma        = permute(rowvec(Sigma), [1 3 2]);  
end

nstates         = size(Sigma, 3); 
if size(mu, 2) ~= nstates && size(mu, 1) == nstates
    mu          = mu';
end

assert(size(mu, 2) == dGauss, ...
    'Number of gaussian dimensions does not match size of mu!');
assert(size(T, 3) == dDiscrete, ...
    'Number of discrete dimensions does not match size of T!');

CPD             = structure(mu, Sigma, nstates, gaussPrior, T, ...
    discretePrior, distType, dGauss, dDiscrete); 

CPD.cpdType     = 'condMixes'; 
CPD.fitFn       = @condMixedCpdFit; 
CPD.fitFnEss    = @condMixedCpdFitEss;
CPD.essFn       = @condMixedCpdComputeEss;
CPD.logPriorFn  = @logPriorFn;
CPD.rndInitFn   = @rndInit;

end

function cpd = rndInit(cpd)
%% randomly initialize
dGauss          = cpd.dGauss;
nstates         = cpd.nstates;

cpd.mu          = randn(dGauss, nstates);
Sigma           = zeros(dGauss, dGauss, nstates);
for i=1:nstates
    Sigma(:, :, i)      = randpd(d) + 2*eye(d); 
end

cpd.Sigma       = Sigma;

cpd.T           = normalize(rand(size(CPT.T), 2));

end

function logp = logPriorFn(cpd)
%% TODO: Fix to work for both the Gaussian and the discrete distribution!!

% calculate the gaussian logprior
logp        = 0;
prior       = cpd.gaussPrior; 
if ~isempty(prior) && isstruct(prior)
    nstates     = cpd.nstates; 
    mu          = cpd.mu;
    Sigma       = cpd.Sigma; 
    for k = 1:nstates
    logp        = logp + gaussInvWishartLogprob(...
        prior, mu(:, k), Sigma(:, :, k));
    end
end

end

function ess = condMixedCpdComputeEss(cpd, data, weights, B)
%% Compute the expected sufficient statistics for a condGaussCpd
% data is nobs-by-d
% weights is nobs-by-nstates; the marginal probability of the discrete
% parent for each observation. 
% B is ignored, but required by the interface, (since mixture emissions use
% it). 
%%
d               = cpd.d; 
nstates         = cpd.nstates; 
wsum            = sum(weights, 1);
xbar            = bsxfun(@rdivide, data'*weights, wsum); %d-by-nstates
XX              = zeros(d, d, nstates);
for j=1:nstates
    Xc          = bsxfun(@minus, data, xbar(:, j)');
    XX(:, :, j) = bsxfun(@times, Xc, weights(:, j))'*Xc;
end
ess = structure(xbar, XX, wsum); 
end

function cpd = condGaussCpdFitEss(cpd, ess)
%% Fit a condGaussCpd given expected sufficient statistics
% ess is a struct containing wsum, XX, and xbar
% cpd is a condGaussCpd as created by e.g condGaussCpdCreate
%
%%
wsum    = ess.wsum;
XX      = ess.XX;
xbar    = ess.xbar;
d       = cpd.d;
nstates = cpd.nstates;
prior   = cpd.prior;
if ~isstruct(prior) || isempty(prior) % do mle
    cpd.mu    = reshape(xbar, d, nstates);
    cpd.Sigma = bsxfun(@rdivide, XX, reshape(wsum, [1 1 nstates]));
else % do map
    kappa0 = prior.k;
    m0     = prior.mu(:);
    nu0    = prior.dof;
    S0     = prior.Sigma;
    mu     = zeros(d, nstates);
    Sigma  = zeros(d, d, nstates);
    for k = 1:nstates
        xbark          = xbar(:, k);
        XXk            = XX(:, :, k);
        wk             = wsum(k);
        mn             = (wk*xbark + kappa0*m0)./(wk + kappa0);
        a              = (kappa0*wk)./(kappa0 + wk);
        b              = nu0 + wk + d + 2;
        Sprior         = (xbark-m0)*(xbark-m0)';
        Sigma(:, :, k) = (S0 + XXk + a*Sprior)./b;
        mu(:, k)       = mn;
    end
    cpd.mu    = mu;
    cpd.Sigma = Sigma;
end
end

function cpd = condMixedCpdFit(cpd, Z, Y)
%% Fit a conditional Gaussian CPD
% Z(i) is the state of the parent Z in observation i.
% Y(i, :) is the ith 1-by-d observation of the child corresponding to Z(i)
% 
% By default we lightly regularize the parameters so we are doing map
% estimation, not mle. The Gauss-invWishart prior is set by 
% condGaussCpdCreate. 
%
% cpd.mu is a matrix of size dGauss-by-nstates
% cpd.Sigma is of size dGauss-by-dGauss-by-nstates
% cpd.T is a matrix of nstates-by-nobsstates-by-dDiscrete

dGauss      = cpd.dGauss; 
Z           = colvec(Z); 
nstates     = cpd.nstates; 

%% do Gaussian first; support only MLE for now, implement MAP for nans later
prior       = cpd.gaussPrior; 
mu          = zeros(dGauss, nstates);
Sigma       = zeros(dGauss, dGauss, nstates);

for k=1:nStates
    mask            = Z == k;
    mu(:, k)        = nanmean(Y(mask, :), 1);
    Sigma(:, :, k)  = nancov(Y(mask, :));
end

cpd.mu          = mu;
cpd.Sigma       = Sigma;

%% do Discrete dimensions now
T               = cpd.T;
dDiscrete       = cpd.dDiscrete;
nObsStates      = size(T, 2);
prior           = cpd.discretePrior;

if isempty(CPD.prior)
    alpha   = 1;
else
    alpha   = prior.alpha;
end
for k = 1:nstates
    % each matrix is size nObStates x dDiscrete
    T(k, :, :)      = reshape(normalize(histc(Y(Z==k, :), 1:nObsStates) + ...
        alpha - 1, 1), [1 nObsStates dDiscrete]);
end

end
