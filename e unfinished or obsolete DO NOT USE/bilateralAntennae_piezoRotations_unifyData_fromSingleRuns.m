%% load and average all data

clear
flyNum     = 263;
cellNum    = 1;
% basename = 'fly233_cell03'; % temp


datafiles = uipickfiles('FilterSpec', sprintf('D:/Dropbox (HMS)/p2/fly%3d_PP/data_*', flyNum));

for i = 1:length(datafiles)
    data(i) = load(datafiles{i});
end

meanFR = cat(1, data.meanFR);
meanfiltVm = cat(1, data.meanfiltVm);

meanFR = mean(meanFR);
meanfiltVm = mean(meanfiltVm);

% metadata = data(1).metadata;
aL = data(1).aL;
aR = data(1).aR;
intensities = data(1).intensities;
directions = data(1).directions;
% trialIndicesFullBlock = data(1).trialIndicesFullBlock;

%% replot single scatters - new
% for i = 2:4
% figure; hold on; axis image
% xlabel(' Left antenna (ipsi) - angular displ')
% ylabel(' Right antenna (contra) - angular displ')
% scatter(aL, aR, 120, data(i).meanfiltVm-data(i).meanfiltVm(1), 'filled', 'MarkerEdgeColor', 'k')
% title(sprintf('%s - average Vm', data(1).metadata.experimentHandle.basename), 'Interpreter', 'none')
% colormap(bluewhitered(256)), cb = colorbar;
% cb.TickLabels = num2str((cb.Ticks+data(i).meanfiltVm(1))', '%2.1f');
% cb.Label.String = '(mV)';
% export_fig(sprintf('scatter_Vm_%s.pdf', data(i).metadata.experimentHandle.basename))
% 
% 
% figure; hold on; axis image
% xlabel(' Left antenna (ipsi) - angular displ')
% ylabel(' Right antenna (contra) - angular displ')
% scatter(aL, aR, 120, data(i).meanFR-data(i).meanFR(1), 'filled', 'MarkerEdgeColor', 'k')
% title(sprintf('%s - firing rate', data(i).metadata.experimentHandle.basename), 'Interpreter', 'none')
% colormap(bluewhitered(256)), cb = colorbar;
% cb.TickLabels = num2str((cb.Ticks+data(i).meanFR(1))', '%2.1f');
% cb.Label.String = '(Hz)';
% export_fig(sprintf('scatter_FR_%s.pdf', data(i).metadata.experimentHandle.basename))
% end
% 
for i = 1:length(data)


figure; hold on; axis image
xlabel(' Left antenna (ipsi) - angular displ')
ylabel(' Right antenna (contra) - angular displ')
scatter(aL, aR, 120, data(i).meanfiltVm-data(i).meanfiltVm(1), 'filled', 'MarkerEdgeColor', 'k')
title(sprintf('Vm (averaged within block)\n%s\n%s', data(i).metadata.experimentHandle.basename, data(i).metadata.key.ID), 'Interpreter', 'none')
colormap(bluewhitered(256)), cb = colorbar;
cb.TickLabels = num2str((cb.Ticks+data(i).meanfiltVm(1))', '%2.1f');
cb.Label.String = '(mV)';
export_fig(sprintf('scatter_Vm_%s.pdf', data(i).metadata.key.ID))


figure; hold on; axis image
xlabel(' Left antenna (ipsi) - angular displ')
ylabel(' Right antenna (contra) - angular displ')
scatter(aL, aR, 120, data(i).meanFR-data(i).meanFR(1), 'filled', 'MarkerEdgeColor', 'k')
title(sprintf('firing rate (averaged within block)\n%s\n%s', data(i).metadata.experimentHandle.basename, data(i).metadata.key.ID), 'Interpreter', 'none')
colormap(bluewhitered(256)), cb = colorbar;
cb.TickLabels = num2str((cb.Ticks+data(i).meanFR(1))', '%2.1f');
cb.Label.String = '(Hz)';
export_fig(sprintf('scatter_FR_%s.pdf', data(i).metadata.key.ID))

end

%% consider only single antenna movements - split in 5x blocks and plot over time
figure
subplot(2,2,1); hold on
title(sprintf('L (ispi) antenna only\n7 sequential blocks of 5 reps'), 'FontSize', 12 )
ylabel('Vm (mV)')
ax(1) = gca;


ixL = find(aR==0);
ixL = ixL([2:4,1,5:7]); 
aL_ixL = aL(ixL);   % xaxis - intervals are not linear, plus values might not be sorted? - but they are

filtVmZero = [];
filtVmSt = [];
for bl = 1:length(data)
    filtVmZero = cat(1, filtVmZero, data(bl).filtVm{1});
    filtVmSt = cat(1, filtVmSt, cat(2,data(bl).filtVm{setdiff(ixL,1)}) );
end
filtVmZero = reshape(filtVmZero, [], 7);
filtVmZero = squeeze(mean(filtVmZero,1));
filtVmZero = filtVmZero';

filtVmSt = permute(filtVmSt, [1,3,2]);
filtVmSt = reshape(filtVmSt, 5,7,6);
filtVmSt = squeeze(mean(filtVmSt,1));
allBlocksL = cat(2, filtVmSt(:,1:3), filtVmZero, filtVmSt(:,4:6)); %time is within rows, each column is a stim

plot(aL_ixL, allBlocksL');


subplot(2,2,2); hold on
title(sprintf('R (contra) antenna only\n7 sequential blocks of 5 reps'), 'FontSize', 12 )
ax(2) = gca;
ixR = find(aL==0);
ixR = ixR([2:4,1,5:7]); 
aR_ixR = aR(ixR);   % xaxis - intervals are not linear, plus values might not be sorted? - but they are

filtVmZero = [];
filtVmSt = [];
for bl = 1:length(data)
    filtVmZero = cat(1, filtVmZero, data(bl).filtVm{1});
    filtVmSt = cat(1, filtVmSt, cat(2,data(bl).filtVm{setdiff(ixR,1)}) );
end
filtVmZero = reshape(filtVmZero, [], 7);
filtVmZero = squeeze(mean(filtVmZero,1));
filtVmZero = filtVmZero';

filtVmSt = permute(filtVmSt, [1,3,2]);
filtVmSt = reshape(filtVmSt, 5,7,6);
filtVmSt = squeeze(mean(filtVmSt,1));
allBlocksR = cat(2, filtVmSt(:,1:3), filtVmZero, filtVmSt(:,4:6)); %time is within rows, each column is a stim
plot(aR_ixR, allBlocksR');

legend(num2str((1:7)'))




subplot(2,2,3); hold on
plot(aL_ixL, mean(allBlocksL), '-k')
title(sprintf('L (ipsi) antenna only\naverage'), 'FontSize', 12 )
xlabel('antennal rotation (deg)')
ylabel('FR (Hz)')
ax(3) = gca;


subplot(2,2,4); hold on
plot(aR_ixR, mean(allBlocksR), '-k')
xlabel('antennal rotation (deg)')
ylabel('FR (Hz)')
title(sprintf('R (contra) antenna only\naverage'), 'FontSize', 12 )
ax(4) = gca;

linkaxes(ax, 'y')


%%
export_fig(sprintf('AVG_runs1to4_singleAntennae_Vm_%s.pdf', data(1).metadata.experimentHandle.basename))



%% consider only single antenna movements - split in 5x blocks and plot over time --------FR
figure
subplot(2,2,1); hold on
title(sprintf('L (ispi) antenna only\n7 sequential blocks of 5 reps'), 'FontSize', 12 )
ylabel('FR (Hz)')
ax(1) = gca;


ixL = find(aR==0);
ixL = ixL([2:4,1,5:7]); 
aL_ixL = aL(ixL);   % xaxis - intervals are not linear, plus values might not be sorted? - but they are

filtVmZero = [];
filtVmSt = [];
for bl = 1:length(data)
    filtVmZero = cat(1, filtVmZero, data(bl).FR{1});
    filtVmSt = cat(1, filtVmSt, cat(2,data(bl).FR{setdiff(ixL,1)}) );
end
filtVmZero = reshape(filtVmZero, [], 7);
filtVmZero = squeeze(mean(filtVmZero,1));
filtVmZero = filtVmZero';

filtVmSt = permute(filtVmSt, [1,3,2]);
filtVmSt = reshape(filtVmSt, 5,7,6);
filtVmSt = squeeze(mean(filtVmSt,1));
allBlocksL = cat(2, filtVmSt(:,1:3), filtVmZero, filtVmSt(:,4:6)); %time is within rows, each column is a stim

plot(aL_ixL, allBlocksL');


subplot(2,2,2); hold on
title(sprintf('R (contra) antenna only\n7 sequential blocks of 5 reps'), 'FontSize', 12 )
ax(2) = gca;
ixR = find(aL==0);
ixR = ixR([2:4,1,5:7]); 
aR_ixR = aR(ixR);   % xaxis - intervals are not linear, plus values might not be sorted? - but they are

filtVmZero = [];
filtVmSt = [];
for bl = 1:length(data)
    filtVmZero = cat(1, filtVmZero, data(bl).FR{1});
    filtVmSt = cat(1, filtVmSt, cat(2,data(bl).FR{setdiff(ixR,1)}) );
end
filtVmZero = reshape(filtVmZero, [], 7);
filtVmZero = squeeze(mean(filtVmZero,1));
filtVmZero = filtVmZero';

filtVmSt = permute(filtVmSt, [1,3,2]);
filtVmSt = reshape(filtVmSt, 5,7,6);
filtVmSt = squeeze(mean(filtVmSt,1));
allBlocksR = cat(2, filtVmSt(:,1:3), filtVmZero, filtVmSt(:,4:6)); %time is within rows, each column is a stim

plot(aR_ixR, allBlocksR');

legend(num2str((1:7)'))



subplot(2,2,3); hold on
plot(aL_ixL, mean(allBlocksL), '-k')
title(sprintf('L (ipsi) antenna only\naverage'), 'FontSize', 12 )
xlabel('antennal rotation (deg)')
ylabel('FR (Hz)')
ax(3) = gca;


subplot(2,2,4); hold on
plot(aR_ixR, mean(allBlocksR), '-k')
xlabel('antennal rotation (deg)')
ylabel('FR (Hz)')
title(sprintf('R (contra) antenna only\naverage'), 'FontSize', 12 )
ax(4) = gca;

linkaxes(ax, 'y')


%%
export_fig(sprintf('AVG_runs1to4_singleAntennae_FR_%s.pdf', data(1).metadata.experimentHandle.basename))

%% revert this to general scatter!!

for i = 1:length(data)


figure; hold on; axis image
xlabel(' Left antenna (ipsi) - angular displ')
ylabel(' Right antenna (contra) - angular displ')
scatter(aL, aR, 120, data(i).meanfiltVm-data(i).meanfiltVm(1), 'filled', 'MarkerEdgeColor', 'k')
title(sprintf('Vm (averaged within block)\n%s\n%s', data(i).metadata.experimentHandle.basename, data(i).metadata.key.ID), 'Interpreter', 'none')
colormap(bluewhitered(256)), cb = colorbar;
cb.TickLabels = num2str((cb.Ticks+data(i).meanfiltVm(1))', '%2.1f');
cb.Label.String = '(mV)';
export_fig(sprintf('scatter_Vm_%s.pdf', data(i).metadata.key.ID))


figure; hold on; axis image
xlabel(' Left antenna (ipsi) - angular displ')
ylabel(' Right antenna (contra) - angular displ')
scatter(aL, aR, 120, data(i).meanFR-data(i).meanFR(1), 'filled', 'MarkerEdgeColor', 'k')
title(sprintf('firing rate (averaged within block)\n%s\n%s', data(i).metadata.experimentHandle.basename, data(i).metadata.key.ID), 'Interpreter', 'none')
colormap(bluewhitered(256)), cb = colorbar;
cb.TickLabels = num2str((cb.Ticks+data(i).meanFR(1))', '%2.1f');
cb.Label.String = '(Hz)';
export_fig(sprintf('scatter_FR_%s.pdf', data(i).metadata.key.ID))

end

%% interpolate scatter
% bring aL and aR in mash grid format
V = meanfiltVm'-meanfiltVm(1);
F = scatteredInterpolant(aL,aR,V);

[Xq,Yq] = meshgrid(min(aL):0.25:max(aL));
Vq = F(Xq, Yq);

figure; 
xlabel('L (deg)')
ylabel('R (deg)')
surf(Xq,Yq,Vq)

az = 0;
el = 90;
view(az, el);

ylim([min(aL), max(aL)])
xlim([min(aL), max(aL)])
axis square

colormap(bluewhitered(256)), cb = colorbar;
cb.TickLabels = num2str((cb.Ticks+meanfiltVm(1))', '%2.1f');
cb.Label.String = '(mV)';

title(data(1).metadata.experimentHandle.basename, 'Interpreter', 'none')

export_fig(fullfile(data(1).metadata.experimentHandle.flyfolder,'interpolated_scatter_Vm.jpg'))
savefig(fullfile(data(1).metadata.experimentHandle.flyfolder,'interpolated_scatter_Vm.fig'))


% Vq = interp2(aL, aR, meanfiltVm, Xq,Yq); % this does not work for unevenly sampled data (or non-stricly monotonically sampled)



%% make average direction tuning curves
I30 = intensities == 30; %patch
I60 = intensities == 60;
I120 = intensities == 120;
f = figure; hold on
plot(directions(I30), meanfiltVm(I30)-meanfiltVm(1)); 
plot(directions(I60), meanfiltVm(I60)-meanfiltVm(1));
plot(directions(I120), meanfiltVm(I120)-meanfiltVm(1));
legend('30', '60', '120')
legend boxoff
legend('Location', 'NorthWest')
xlabel('estimated wind direction')
ylabel('mean filtered Vm change with respect to rest')
title({data(1).metadata.experimentHandle.basename}, 'Interpreter', 'none')
f.WindowStyle = 'normal';
f.Position = [1922 516 1019 364];
savefig(sprintf('AVG_dirTuning_filteredVoltageRestSubtr_%s', data(1).metadata.experimentHandle.basename))
export_fig(sprintf('AVG_dirTuning_filteredVoltageRestSubtr_%s.pdf', data(1).metadata.experimentHandle.basename))


f = figure; hold on
plot(directions(I30), meanFR(I30)-meanFR(1)); 
plot(directions(I60), meanFR(I60)-meanFR(1));
plot(directions(I120), meanFR(I120)-meanFR(1));
legend('30', '60', '120')
legend boxoff
legend('Location', 'NorthWest')
xlabel('estimated wind direction')
ylabel('mean firing rate change with respect to rest')
title({data(1).metadata.experimentHandle.basename}, 'Interpreter', 'none')
f.WindowStyle = 'normal';
f.Position = [1922 516 1019 364];

savefig(sprintf('dirTuning_firingRateRestSubtr_%s', data(1).metadata.experimentHandle.basename))
export_fig(sprintf('dirTuning_firingRateRestSubtr_%s.pdf', data(1).metadata.experimentHandle.basename))


