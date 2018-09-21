function [phiGrouped, data, taskId, phiNames, phiType, uniqueGroups] = ...
    getMultiModalData(varargin)

[   minVisits                                       , ...
    groups                                          , ...
    scale               ] = process_options(varargin, ...
    'minVisits'         , 1                         , ...
    'groups'            , {{'MRI', 'DEMOGRAPHIC', 'GENETIC'}, ...
    {'FDG'}, {'CSF'}, {'AMYLOID'}}                  , ...
    'scale'             , true                      );

if exist('../data/', 'dir') == 7
    data_path   = sprintf('../data/all_data.csv');
    dict_path   = sprintf('../data/all_data_map.csv');
elseif exist('/phobos/alzheimers', 'dir') == 7
    data_path   = sprintf('/phobos/alzheimers/adni/all_data.csv');
    dict_path   = sprintf('/phobos/alzheimers/adni/all_data_map.csv');
elseif exist('/z/data', 'dir') == 7    
    data_path   = sprintf('/z/data/alzheimers/adni/all_data.csv');
    dict_path   = sprintf('/z/data/alzheimers/adni/all_data_map.csv');
else
    fprintf('Multi-modal data file not found!\n');
end

data                = readtable(data_path);
dataDict            = readtable(dict_path, 'ReadVariableNames', false);
dataDict            = dataDict.Var2;
[isSorted, data]    = isDataSorted(data);

if ~isSorted
    fprintf('RID not sorted... Sorting table!\n');
end

% Make diagnosis/gender numeric
dx                              = zeros(size(data.DX));
dx(strcmpi(data.DX, 'NL'))      = 1;
dx(strcmpi(data.DX, 'MCI'))     = 2;
dx(strcmpi(data.DX, 'AD'))      = 3;
data.DX                         = dx;

data            = removeNoise(data, minVisits);
columns         = data.Properties.VariableNames;

n               = size(data, 1);
isinGroup       = false(n, length(groups));
phiGrouped      = cell(length(groups), 1);
phiNames        = cell(length(groups), 1);
phiType         = dataDict;

isNumeric       = colvec(arrayfun(@(ii)isnumeric(data{:, ii}), ...
    1:length(dataDict)));

for g=1:length(groups)
    modalities      = groups{g};
    groupMask       = strcmpi(data{:, cellfun(@(seq)cat(2, 'HAS_', seq), ...
        modalities, 'UniformOutput', false)}, 'True');
    isinGroup(:, g)     = all(groupMask, 2);
    
    phiLocation     = cellfun(@(seq)strcmpi(dataDict, seq), modalities, ...
        'UniformOutput', false);
    phiLocation     = any(cat(2, phiLocation{:}), 2);
    
    phiGrouped{g}   = data{:, phiLocation & isNumeric};
    phiNames{g}     = columns(phiLocation & isNumeric);
end

if scale
    phiGrouped      = cellfun(@zScale, phiGrouped, 'UniformOutput', false);
end

[uniqueGroups, ~, taskId]   = unique(isinGroup, 'rows');

end