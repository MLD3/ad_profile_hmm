function [phi, phiSum, summaryNames] = getMriSummary(data, varargin)

[   scale       ,                               ...
    imputeNan   ,                               ...
    age         ,                               ...
    normIcv     ,                               ...
    counts      ,                               ...
    phiNames    ] = process_options(varargin,   ...
    'scale'     ,   true                    ,   ...
    'imputeNan' ,   true                    ,   ...
    'age'       ,   true                    ,   ...
    'normIcv'   ,   true                    ,   ...
    'counts'    ,   []                      ,   ...
    'phiNames'  ,   {}                      );

% hippocampus, entorhinal, fusiform, midtemp, ventricles, wholebrain

summaryNames    = {'HIPPOCAMPUS', 'ENTORHINAL', 'FUSIFORM', 'MIDTEMP', ...
    'VENTRICLES'}; %, 'WHOLEBRAIN'};

if ~istable(data)
    assert(~isempty(phiNames), ...
        'Must provide features names to retrieve summary');
end

phiCodes        = {{'ST29SV', 'ST88SV'}, ...
    {'ST24CV', 'ST83CV'}, ...
    {'ST26CV', 'ST85CV'}, ...
    {'ST40CV', 'ST99CV'}, ...
    {'ST30SV', 'ST37SV', 'ST89SV', 'ST96SV'}};, ...
%     {'ST128SV', 'ST17SV', 'ST18SV', 'ST61SV', 'ST16SV', 'ST53SV', ...
%     'ST42SV', 'ST29SV', 'ST12SV', 'ST11SV', 'ST65SV', 'ST76SV', ...
%     'ST77SV', 'ST120SV', 'ST75SV', 'ST112SV', 'ST101SV', 'ST88SV', ...
%     'ST71SV', 'ST70SV', 'ST124SV', 'ST147SV', 'ST148SV', 'ST150SV', ...
%     'ST151SV'}};

phiSum      = zeros(size(data, 1), length(phiCodes));

for regIdx=1:length(phiCodes)
    region      = phiCodes{regIdx};
    for voxel=region
        if istable(data)
            col     = data{:, voxel};
        elseif ismatrix(data)
            col     = data(:, strcmpi(phiNames, voxel));
        end
        if imputeNan
            col(isnan(col))     = nanmean(col);
        end
        phiSum(:, regIdx)   = phiSum(:, regIdx) + col;
    end
end

if normIcv
    if istable(data)
        icv         = data{:, 'ICV'};
    elseif ismatrix(data)
        if any(strcmpi(phiNames, 'ICV'))
            icv     = data(:, strcmpi(phiNames, 'ICV'));
        else
            warning('No ICV feature found! Skipping normalization...');
            normIcv         = false;
        end
    end
    if normIcv
        if imputeNan
            icv(isnan(icv))         = nanmean(icv);
        end
        phiSum      = bsxfun(@rdivide, phiSum, icv);
    end
end

% phiSum(phiSum > 6)      = 6;
% phiSum(phiSum < -6)     = -6;

if age
    summaryNames    = cat(2, summaryNames, 'AGE');
    if istable(data)
        phiSum      = [phiSum data.AGE];
    elseif ismatrix(data)
        phiSum      = [phiSum data(:, strcmpi(phiNames, 'AGE'))];
    end
end

if scale
    phiSum      = zScale(phiSum);
end

if ~isempty(counts)
    phi         = mat2cell(phiSum, counts);
else
    phi         = phiSum;
end

end