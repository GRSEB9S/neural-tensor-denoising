function [ summary ] = tensorDenoiseWrapper( Y, r1, method, dataset )
% wrapper for tensorDenoiseGridSearchCV.m
% iterate over trial count and method

options = [];
options.gridStep = 2;
options.minRank = [20 2 2];
options.maxRank = [80 20 20];
options.resample = 1;
options.verbose = 0;

options.r1 = r1;
options.method = method;

summary = tensorDenoiseGridSearchCV(Y, options);

%% write summary to file
L = length(summary.err);
summary.r1 = repmat(r1,L,1);
summary.minrank = repmat(summary.minrank,L,1);

table_out = struct2table(rmfield(summary,{'options'}));
writetable(table_out, ['dat-files/dat-' num2str(r1, '%04d') '-' num2str(method) '-' num2str(dataset)]);

end