%% specify
datindex = 1;
dataset = 2;

%% load
files_r = dir('dat-file*');
datpath = [files_r(datindex).name '/'];

files = dir(datpath);
files = files(~[files.isdir]);
n = length(files);
load([datpath files(1).name]);

load([datpath 'plotfiles']);

if summary.options.simrun
  rnk = summary.options.rnks(dataset,:);
  disp(rnk)
end

%% plot
h = figure; 
ax = subplot(2,1,1);
hold all
xaxis = 2:size(errPlot{dataset},1)+1;
plot(xaxis,mean(errPlot{dataset},3),'linewidth',2);
ax.YLim(1) = 0;
legend('tensor','neuron','time','condition','average');
xlabel('trials (per neuron per condition)');
ylabel('relative error');
title('main analysis | relative error vs trial')
h.Color = 'w';

%h = figure; 
ax = subplot(2,1,2);
hold all
plot(xaxis,mean(minRankPlot{dataset}{1},3),'linewidth',2)
ax.ColorOrderIndex = 1;
plot(xaxis,mean(minRankPlot{dataset}{2}(:,1,:),3),'--')
plot(xaxis,mean(minRankPlot{dataset}{3}(:,2,:),3),'--')
plot(xaxis,mean(minRankPlot{dataset}{4}(:,3,:),3),'--')

ax.ColorOrderIndex = 1;
if summary.options.simrun
  plot([xaxis(1) xaxis(end)], [rnk(1) rnk(1)])
  plot([xaxis(1) xaxis(end)], [rnk(2) rnk(2)])
  plot([xaxis(1) xaxis(end)], [rnk(3) rnk(3)]) 
end

legend('neuron','time','condition');
xlabel('trials (per neuron per condition)');
ylabel('rank');
title('main analysis | optimal rank vs trial')
h.Color = 'w';

