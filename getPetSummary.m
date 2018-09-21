function [phi, phiSum, regionNames]  = getPetSummary(data, varargin)

[   scale                                   , ...
    imputeNan                               , ...
    phiNames                                , ...
    counts                                  , ...
    age         ] = process_options(varargin, ...
    'scale'     ,   true                    , ...
    'imputeNan' ,   true                    , ...
    'phiNames'  ,   {}                      , ...
    'counts'    ,   []                      , ...
    'age'       ,   true                    );

if ~istable(data)
    assert(~isempty(phiNames), ...
        'Must provide features names to retrieve summary');
end

regionNames     = {'Angular_Left', 'Angular_Right', ...
    'CingulumPost_Bilateral', 'Temporal_Left', 'Temporal_Right'};
phiTypes        = {'MEAN'};

[region, typ]   = ndgrid(1:length(regionNames), 1:length(phiTypes));
phiCodes        = arrayfun(@(ii)strcat(regionNames{region(ii)}, '_', ...
    phiTypes{typ(ii)}), 1:numel(region), 'UniformOutput', false);

if istable(data)
    phiSum      = table2array(data(:, phiCodes));
elseif ismatrix(data)
    include     = cellfun(@(seq)any(strcmpi(phiCodes, seq)), phiNames);
    phiSum      = data(:, include);
end

if age
    regionNames     = cat(2, regionNames, 'AGE');
    if istable(data)
        phiSum      = [phiSum data.AGE];
    elseif ismatrix(data)
        phiSum      = [phiSum data(:, strcmpi(phiNames, 'AGE'))];
    end
end

if imputeNan
    % Replace NaN values with the mean of the column
    for ii=1:size(phiSum, 2)
        col                         = phiSum(:, ii);
        phiSum(isnan(col), ii)      = mean(col(~isnan(col)));
    end
end

phiSum          = mean(phiSum, 2);
regionNames     = {'FDG_SUVR'};

if scale
    phiSum      = zScale(phiSum);
end

% phiSum(phiSum > 6)      = 6;
% phiSum(phiSum < -6)     = -6;

if ~isempty(counts)
    phi         = mat2cell(phiSum, counts);
else
    phi         = phiSum;
end

end