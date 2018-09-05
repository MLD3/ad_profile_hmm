load profile-hmm-data.mat;

phiIdx          = strcmpi(allData.phiNames, 'ADAS');
nonMissing      = ~isMissing(:, phiIdx);

data            = allData.phiUnscaled(nonMissing, phiIdx);
data(data < 0)  = median(data);
dx              = allData.data.DX(nonMissing);

wnData          = data + normrnd(0, 0.1, length(data), 1);

% peak            = max(data); 
% trough          = min(data); 
% res             = (peak-trough)/(N-1);
% lin_array       = trough:res:peak;

% sortedData      = sort(data);

% transform to uniform distributed data
% uniform_data    =interp1(sortedData, lin_array, data);  

% u1_data         = uniform_data/(max(uniform_data)-min(uniform_data));  
% u1_data         = u1_data - min(u1_data);

% gauss_data      = norminv(u1_data, 0, 1);

order           = tiedrank(wdata);
unifProb        = order/(length(order) + 1);
transformed     = norminv(unifProb, 0, 1);

figure; 
subplot(221);
histogram(wdata);
grid on; box on;
subplot(222);
histogram(transformed);
grid on; box on;
subplot(223);
boxplot(wdata, allData.data.DX(nonMissing));
grid on;box on;
subplot(224);
