function CPD = condGaussMissingCpdCreate(mu, Sigma, varargin)
%% Create a conditional Gaussian distribution for use in a graphical model
%
% mu is a matrix of size d-by-nstates
% Sigma is of size d-by-d-by-nstates
%
% 'prior' is a Gauss-inverseWishart distribution, namely, a struct with
% fields  mu, Sigma, dof, k
% Set 'prior' to 'none' to do mle. 
%%

% This file is from pmtk3.googlecode.com

prior = process_options(varargin, 'prior', []); 
d = size(Sigma, 1); 
if isempty(prior) 
      prior.mu    = zeros(1, d);
      prior.Sigma = 0.1*eye(d);
      prior.k     = 0.01;
      prior.dof   = d + 1; 
end

if isvector(Sigma)
   Sigma = permute(rowvec(Sigma), [1 3 2]);  
end

nstates = size(Sigma, 3); 

if size(mu, 2) ~= nstates && size(mu, 1) == nstates
    mu = mu';
end

CPD = structure(mu, Sigma, nstates, d, prior); 
CPD.cpdType    = 'condMissingGauss'; 
CPD.fitFn      = @condGaussMissingCpdFit; 
CPD.fitFnEss   = @condGaussMissingCpdFitEss;
CPD.essFn      = @condGaussMissingCpdComputeEss;

end

function ess = condGaussMissingCpdComputeEss(cpd, data, weights, B)
%% Compute the expected sufficient statistics for a condGaussCpd
% data is nobs-by-d
% weights is nobs-by-nstates; the marginal probability of the discrete
% parent for each observation. 
% B is ignored, but required by the interface, (since mixture emissions use
% it). 
%%
d           = cpd.d; 
nstates     = cpd.nstates; 
wsum        = sum(weights, 1);

xbar        = zeros(d, nstates);
XX          = zeros(d, d, nstates);

for ii=1:d
    nonMissing      = ~isnan(data(:, ii));
    xbar(ii, :)     = bsxfun(@rdivide, data(nonMissing, ii)'*...
        weights(nonMissing, :), sum(weights(nonMissing, :)));
end

for k=1:nstates
    for ii=1:d
        nonMissing  = ~isnan(data(:, ii));
        
        Xc              = data(nonMissing, ii) - xbar(ii, k);
        XX(ii, ii, k)   = (Xc.*weights(nonMissing, k))'*Xc;
        
        XX(ii, ii, k)   = XX(ii, ii, k)/sum(weights(nonMissing, k));
    end
end

ess = structure(xbar, XX, wsum); 
end

function cpd = condGaussMissingCpdFitEss(cpd, ess)
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
    cpd.mu      = reshape(xbar, d, nstates);
    cpd.Sigma   = XX;
%     cpd.Sigma = bsxfun(@rdivide, XX, reshape(wsum, [1 1 nstates]));
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

function cpd = condGaussMissingCpdFit(cpd, Z, Y)
%% Fit a conditional Gaussian CPD
% Z(i) is the state of the parent Z in observation i.
% Y(i, :) is the ith 1-by-d observation of the child corresponding to Z(i)
% 
% By default we lightly regularize the parameters so we are doing map
% estimation, not mle. The Gauss-invWishart prior is set by 
% condGaussCpdCreate. 
%
%  cpd.mu is a matrix of size d-by-nstates
%  cpd.Sigma is of size d-by-d-by-nstates
%%
d = cpd.d; 
Z = colvec(Z); 
prior = cpd.prior; 
nstates = cpd.nstates; 
if ~isstruct(prior) || isempty(prior) % do mle
    cpd.mu    = partitionedMean(Y, Z, nstates)';
    cpd.Sigma = partitionedCov(Y, Z,  nstates); 
else  % map
    mu = zeros(d, nstates);
    Sigma = zeros(d, d, nstates);
    for s = 1:nstates
        m              = gaussFit(Y(Z == s, :), prior);
        mu(:, s)       = m.mu(:);
        Sigma(:, :, s) = m.Sigma;
    end
    cpd.mu = mu;
    cpd.Sigma = Sigma; 
end
end
