function [phi, phiSum, summaryNames]  = getGeneticSummary(data, varargin)
%% Get genetic profile
% 1         = 2,2
% 2         = 3,3
% 3         = 4,4
% 4         = 

[   phiNames                                , ...
    counts      ] = process_options(varargin, ...
    'phiNames'  ,   {}                      , ...
    'counts'    ,   []                      );

if ~istable(data)
    assert(~isempty(phiNames), ...
        'Must provide features names to retrieve summary');
end

geneNames       = {'CAT_APGEN1_2.0', 'CAT_APGEN1_3.0', ...
    'CAT_APGEN1_4.0', 'CAT_APGEN2_2.0', 'CAT_APGEN2_3.0', ...
    'CAT_APGEN2_4.0'};
summaryNames    = {'ApoE-Profile'};

if istable(data)
%     apgen           = data{:, geneNames};
    phi             = getGeneProfile(data);
else
    apgen           = [];
    for subgene = rowvec(geneNames)
        apgen       = cat(2, data(:,strcmpi(phiNames, subgene)));
    end
    phi             = getGeneProfile(apgen);
end
% [apgen1, ii]    = find(apgen(:, 1:3)');
% [apgen2, jj]    = find(apgen(:, 4:6)');
% 
% alleles         = nan(size(apgen, 1), 2);
% alleles(ii, 1)  = apgen1+1;
% alleles(jj, 2)  = apgen2+1;
% 
% phi             = nan(size(data, 1), 1);
% 
% for ii=1:length(phi)
%     if isequal(alleles(ii, :), [2 2])
%         phi(ii)     = 1;
%     elseif isequal(alleles(ii, :), [3 3])
%         phi(ii)     = 2;
%     elseif isequal(alleles(ii, :), [4 4])
%         phi(ii)     = 3;
%     elseif isequal(alleles(ii, :), [2 3])
%         phi(ii)     = 4;
%     elseif isequal(alleles(ii, :), [3 4])
%         phi(ii)     = 5;
%     elseif isequal(alleles(ii, :), [2 4])
%         phi(ii)     = 6;
%     end
% end

phiSum      = phi;

if ~isempty(counts)
    phi         = mat2cell(phiSum, counts);
else
    phi         = phiSum;
end

end