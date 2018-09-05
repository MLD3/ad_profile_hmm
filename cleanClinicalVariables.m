function allData = cleanClinicalVariables(allData)
%% Function to clean clinical variables
% Most clinical variables are discrete, which presents a problem when
% trying to model them using normal distributions. A quick hack to solve
% this problem is to simply add white noise with a small standard deviation
% to the data in order to make it continuous

phiNames                = allData.phiNames;
phi                     = allData.phiUnscaled;

% The empirically observed resolution of the clinical scores are:
% ADAS      = 0.33
% TRAILS    = 1
% FAQ       = 1
% RAVLT     = 1
% MMSE      = 1
% CDR       = 0.5

scoreNames          = {'ADAS', 'TRAILS', 'FAQ', 'RAVLT', 'MMSE', 'CDR'};
sigma               = [0.075 0.1 0.1 0.1 0.1 0.1];

for ii=1:length(scoreNames)
    phiIdx          = strcmpi(phiNames, scoreNames{ii});
    
    obsIdx          = ~isnan(phi(:, phiIdx));
    phi(obsIdx, phiIdx)     = phi(obsIdx, phiIdx) + ...
        normrnd(0, sigma(ii), sum(obsIdx), 1);
end

allData.phiUnscaled     = phi;

end