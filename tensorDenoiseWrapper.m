function [ summary ] = tensorDenoiseWrapper( Y, r1, method, dataset, rng_iter, datpath )
% wrapper for tensorDenoiseGridSearchCV.m
% iterate over trial count and method
if nargin < 6
  datpath = 'dat-files/';
end

options = [];
options.gridStep = 3;
options.minRank = [20 2 2];
options.maxRank = [80 20 20];
options.resample = 1;
options.verbose = 1;

options.r1 = r1;
options.method = method;
options.rng_iter = rng_iter;

summary = tensorDenoiseGridSearchCV(Y, options);

%%
table_out = summary;

% provide metadata (UserData? ok)
%table_out.Properties.UserData.r1 = r1;
%table_out.Properties.UserData.method = method;
%table_out.Properties.UserData.rng_iter = rng_iter;
%table_out.Properties.UserData.minrank = summary.minrank;
%table_out.Properties.UserData.minind = summary.minind;

%table_out.Properties.UserData.Ynorm = norm(Y(:)).^2; % handle nans???

filename = [datpath 'dat-' num2str(rng_iter, '%03d') '-' num2str(r1, '%03d') '-' num2str(method) '-' num2str(dataset)];
save(filename, 'table_out'); 

%writetable(table_out, [datpath 'dat-' num2str(r1, '%04d') '-' num2str(method) '-' num2str(dataset)]);

end