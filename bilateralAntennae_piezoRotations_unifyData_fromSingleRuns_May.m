%% load and average all data
clear
flyNum     = 283;
cellNum    = 1;
% basename = 'fly233_cell03'; % temp


datafiles = uipickfiles('FilterSpec', sprintf('/Users/galileo/Dropbox (HMS)/p2/fly%3d_PP/data_*', flyNum));
disp('loading...')
for i = 1:length(datafiles)
    disp(i)
    data(i) = load(datafiles{i});
end
disp('done.')
cd(fileparts(datafiles{1}))

meanFR = cat(1, data.meanFR);
meanfiltVm = cat(1, data.meanfiltVm);

meanFR = mean(meanFR,1);
meanfiltVm = mean(meanfiltVm,1);

% metadata = data(1).metadata;
aL = data(1).aL;
aR = data(1).aR;
intensities = data(1).intensities;
directions = data(1).directions;
% trialIndicesFullBlock = data(1).trialIndicesFullBlock;

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
close

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

savefig(sprintf('AVG_dirTuning_firingRateRestSubtr_%s', data(1).metadata.experimentHandle.basename))
export_fig(sprintf('AVG_dirTuning_firingRateRestSubtr_%s.pdf', data(1).metadata.experimentHandle.basename))
close




%% show single antenna contribution
ixL = find(aR==0);
ixL = ixL([2:4,1,5:7]); 
aL_ixL = aL(ixL);  

ixR = find(aL==0);
ixR = ixR([2:4,1,5:7]); 
aR_ixR = aR(ixR);   % xaxis - intervals are not linear, plus values might not be sorted? - but they are


% want to include errorbars?
filtVmZero = [];
filtVmStL = []; %move L
filtVmStR = []; %move R
FRZero = [];
FRStL = [];
FRStR = [];

for bl = 1:length(data)
    filtVmZero = cat(1, filtVmZero, data(bl).filtVm{1});
    filtVmStL = cat(1, filtVmStL, cat(2,data(bl).filtVm{setdiff(ixL,1)}) );
    filtVmStR = cat(1, filtVmStR, cat(2,data(bl).filtVm{setdiff(ixR,1)}) );
    
    FRZero = cat(1, FRZero, data(bl).FR{1});
    FRStL = cat(1, FRStL, cat(2,data(bl).FR{setdiff(ixL,1)}) );
    FRStR = cat(1, FRStR, cat(2,data(bl).FR{setdiff(ixR,1)}) );
end
if length(directions) == 84
    filtVmZero = mean(reshape(filtVmZero, 7, []));
    filtVmZero = filtVmZero(:);

    FRZero = mean(reshape(FRZero, 7, []));
    FRZero = FRZero(:);
else
    filtVmZero = mean(reshape(filtVmZero, 4, []));
    filtVmZero = filtVmZero(:);

    FRZero = mean(reshape(FRZero, 4, []));
    FRZero = FRZero(:);
end
    
FRsingleT_L = cat(2, FRStL(:,1:3), FRZero, FRStL(:,4:6));
FRsingleT_R = cat(2, FRStR(:,1:3), FRZero, FRStR(:,4:6));
filtVmsingleT_L = cat(2, filtVmStL(:,1:3), filtVmZero, filtVmStL(:,4:6));
filtVmsingleT_R = cat(2, filtVmStR(:,1:3), filtVmZero, filtVmStR(:,4:6));
xs = repmat(aL_ixL', size(FRsingleT_L,1),1);


figure; hold on
h1 = plot(aL_ixL, meanFR(ixL)', '-<');
h2 = plot(aR_ixR, meanFR(ixR)', '->');
% plot(aL_ixL, FRsingleT_L+0.1*rand(size(FRsingleT_L)), 'Color', [0.3 0.3 0.3])
% plot(aR_ixR, FRsingleT_R+0.1*rand(size(FRsingleT_R)), '--k')
plot([aL_ixL(1), aL_ixL(end)], [meanFR(ixL(4)), meanFR(ixL(4))], ':k')
legend([h1, h2],{'L moving, R at rest', 'R moving, L at rest'}, 'Location', 'best')
title(flyNum)
export_fig(sprintf('/Users/galileo/Dropbox (HMS)/OkuboPatellaWilson ms/pp_bits/AVG_singleAntennaeContribution_fly%d.pdf', flyNum))


%% plot it
% figure('Position', [1           1        1190         827]);
% 
% subplot(2,3,1); hold on
% ylabel(sprintf('FR (Hz), N = %d', size(FRsingleT_L,1)))
% notBoxPlot(FRsingleT_L(:), xs(:))
% ax(1) = gca;
% title('L moving, R at rest')
% 
% subplot(2,3,2); hold on
% ax(2) = gca;
% h1 = plot(aL_ixL, meanFR(ixL)', '-<');
% h2 = plot(aR_ixR, meanFR(ixR)', '->');
% % plot(aL_ixL, FRsingleT_L+0.1*rand(size(FRsingleT_L)), 'Color', [0.3 0.3 0.3])
% % plot(aR_ixR, FRsingleT_R+0.1*rand(size(FRsingleT_R)), '--k')
% plot([aL_ixL(1), aL_ixL(end)], [meanFR(ixL(4)), meanFR(ixL(4))], ':k')
% legend([h1, h2],{'L moving, R at rest', 'R moving, L at rest'}, 'Location', 'best')
% % legend boxoff
% 
% subplot(2,3,3); hold on
% notBoxPlot(FRsingleT_R(:), xs(:))
% ax(3) = gca;
% title('R moving, L at rest')
% 
% 
% subplot(2,3,4); hold on
% ylabel('filtered membrane voltage (mV)')
% notBoxPlot(filtVmsingleT_L(:), xs(:))
% ax(4) = gca;
% title('L moving, R at rest')
% 
% 
% subplot(2,3,5); hold on
% ax(5) = gca;
% h2 = plot(aL_ixL, meanfiltVm(ixL)', '-<');
% h3 = plot(aR_ixR, meanfiltVm(ixR)', '->');
% % plot(aL_ixL, filtVmsingleT_L, 'Color', [0.3 0.3 0.3])
% % plot(aR_ixR, filtVmsingleT_R, '--k')
% plot([aL_ixL(1), aL_ixL(end)], [meanfiltVm(ixL(4)), meanfiltVm(ixL(4))], ':k')
% legend([h2, h3],{'L moving, R at rest', 'R moving, L at rest'}, 'Location', 'best')
% % legend boxoff
% 
% subplot(2,3,6); hold on
% notBoxPlot(filtVmsingleT_R(:), xs(:))
% ax(6) = gca;
% title('R moving, L at rest')
% 
% 
% axYLim = cat(1, ax(1:3).YLim);
% axYLim = max(axYLim);
% YLIM_FR = [0, axYLim(2)];
% 
% axYLim = cat(1, ax(4:6).YLim);
% YLIM_filtVm = [min(axYLim(:,1)), max(axYLim(:,2))];
% 
% 
% 
% 
% ax(1).YLim = YLIM_FR;
% linkaxes([ax(1),ax(2),ax(3) ], 'y')
% ax(4).YLim = YLIM_filtVm;
% linkaxes([ax(4),ax(5),ax(6) ], 'y')
% for ia = 1:6
%     ax(ia).XTick = round(aR_ixR*10)/10;
% end
% export_fig(sprintf('AVG_singleAntennaContribution_%s.pdf', data(1).metadata.experimentHandle.basename))




%% make pooled rasterplots
figFolder = fullfile(fileparts(datafiles{1}), sprintf('singleStims_scatterPSTH_%s',data(1).metadata.experimentHandle.basename));
mkdir(figFolder);

clear pooled
for s = 1:length(meanFR)
    prepost2SpikeTimes = [];
    for i = 1:length(datafiles)
        prepost2SpikeTimes = cat(1, prepost2SpikeTimes, data(i).perist(s).prepost2SpikeTimes);
    end
    pooled(s).prepost2SpikeTimes = prepost2SpikeTimes;
end
for s = 1:length(meanFR)
    
    figure('Position', [1141         583         288         150]);
    [xPoints, yPoints] = plotSpikeRaster(pooled(s).prepost2SpikeTimes, 'PlotType', 'vertline');
    hold on
    ax = gca;
    ax.XTick = data(1).T_prepost2(data(1).prepost2_stimDelimitations);
    ax.XAxis.Visible = 'off';
    ax.XGrid = 'on';
    %primo trial in alto, ultimo in basso
    
    fileName = sprintf('AVG_scatterPSTH_spikes_stim%02d_%s.pdf',s, data(1).metadata.experimentHandle.basename);
    export_fig(fullfile(figFolder, fileName))
    close
end


%% general color-coded scatter
figure; hold on; axis image
xlabel(' Left antenna (ipsi) - angular displ')
ylabel(' Right antenna (contra) - angular displ')
scatter(aL, aR, 180, meanfiltVm-meanfiltVm(1), 'filled', 'MarkerEdgeColor', 'k')
title(sprintf('Vm (N = %d)\n%s', size(FRsingleT_L,1), data(1).metadata.experimentHandle.basename), 'Interpreter', 'none')
colormap(bluewhitered(256)), cb = colorbar;
cb.TickLabels = num2str((cb.Ticks+meanfiltVm(1))', '%2.1f');
cb.Label.String = '(mV)';
xlim([-15 15])
ylim([-15 15])
set(gca, 'TickDir', 'out')
export_fig(sprintf('AVG_scatter_Vm_%s.pdf', data(1).metadata.experimentHandle.basename))


figure; hold on; axis image
xlabel(' Left antenna (ipsi) - angular displ')
ylabel(' Right antenna (contra) - angular displ')
scatter(aL, aR, 120, meanFR-meanFR(1), 'filled', 'MarkerEdgeColor', 'k')
title(sprintf('firing rate (N = %d)\n%s', size(FRsingleT_L,1), data(1).metadata.experimentHandle.basename), 'Interpreter', 'none')
colormap(bluewhitered(256)), cb = colorbar;
cb.TickLabels = num2str((cb.Ticks+meanFR(1))', '%2.1f');
cb.Label.String = '(Hz)';
export_fig(sprintf('AVG_scatter_FR_%s.pdf', data(1).metadata.experimentHandle.basename))




%% replot single scatters - only if needed!
% 
% for i = 1:length(data)
% 
% figure; hold on; axis image
% xlabel(' Left antenna (ipsi) - angular displ')
% ylabel(' Right antenna (contra) - angular displ')
% scatter(aL, aR, 120, data(i).meanfiltVm-data(i).meanfiltVm(1), 'filled', 'MarkerEdgeColor', 'k')
% title(sprintf('Vm (averaged within block)\n%s\n%s', data(i).metadata.experimentHandle.basename, data(i).metadata.key.ID), 'Interpreter', 'none')
% colormap(bluewhitered(256)), cb = colorbar;
% cb.TickLabels = num2str((cb.Ticks+data(i).meanfiltVm(1))', '%2.1f');
% cb.Label.String = '(mV)';
% export_fig(sprintf('scatter_Vm_%s.pdf', data(i).metadata.key.ID))
% 
% 
% figure; hold on; axis image
% xlabel(' Left antenna (ipsi) - angular displ')
% ylabel(' Right antenna (contra) - angular displ')
% scatter(aL, aR, 120, data(i).meanFR-data(i).meanFR(1), 'filled', 'MarkerEdgeColor', 'k')
% title(sprintf('firing rate (averaged within block)\n%s\n%s', data(i).metadata.experimentHandle.basename, data(i).metadata.key.ID), 'Interpreter', 'none')
% colormap(bluewhitered(256)), cb = colorbar;
% cb.TickLabels = num2str((cb.Ticks+data(i).meanFR(1))', '%2.1f');
% cb.Label.String = '(Hz)';
% export_fig(sprintf('scatter_FR_%s.pdf', data(i).metadata.key.ID))
% 
% end

%% consider only single antenna movements - split in 5x blocks and plot over time - this was made ad hoc for 35 reps - not general - needs work
% figure
% subplot(2,2,1); hold on
% title(sprintf('L (ispi) antenna only\n7 sequential blocks of 5 reps'), 'FontSize', 12 )
% ylabel('Vm (mV)')
% ax(1) = gca;
% 
% 
% ixL = find(aR==0);
% ixL = ixL([2:4,1,5:7]); 
% aL_ixL = aL(ixL);   % xaxis - intervals are not linear, plus values might not be sorted? - but they are
% 
% filtVmZero = [];
% filtVmSt = [];
% for bl = 1:length(data)
%     filtVmZero = cat(1, filtVmZero, data(bl).filtVm{1});
%     filtVmSt = cat(1, filtVmSt, cat(2,data(bl).filtVm{setdiff(ixL,1)}) );
% end
% filtVmZero = reshape(filtVmZero, [], 7);
% filtVmZero = squeeze(mean(filtVmZero,1));
% filtVmZero = filtVmZero';
% 
% filtVmSt = permute(filtVmSt, [1,3,2]);
% filtVmSt = reshape(filtVmSt, 5,7,6);
% filtVmSt = squeeze(mean(filtVmSt,1));
% allBlocksL = cat(2, filtVmSt(:,1:3), filtVmZero, filtVmSt(:,4:6)); %time is within rows, each column is a stim
% 
% plot(aL_ixL, allBlocksL');
% 
% 
% subplot(2,2,2); hold on
% title(sprintf('R (contra) antenna only\n7 sequential blocks of 5 reps'), 'FontSize', 12 )
% ax(2) = gca;
% ixR = find(aL==0);
% ixR = ixR([2:4,1,5:7]); 
% aR_ixR = aR(ixR);   % xaxis - intervals are not linear, plus values might not be sorted? - but they are
% 
% filtVmZero = [];
% filtVmSt = [];
% for bl = 1:length(data)
%     filtVmZero = cat(1, filtVmZero, data(bl).filtVm{1});
%     filtVmSt = cat(1, filtVmSt, cat(2,data(bl).filtVm{setdiff(ixR,1)}) );
% end
% filtVmZero = reshape(filtVmZero, [], 7);
% filtVmZero = squeeze(mean(filtVmZero,1));
% filtVmZero = filtVmZero';
% 
% filtVmSt = permute(filtVmSt, [1,3,2]);
% filtVmSt = reshape(filtVmSt, 5,7,6);
% filtVmSt = squeeze(mean(filtVmSt,1));
% allBlocksR = cat(2, filtVmSt(:,1:3), filtVmZero, filtVmSt(:,4:6)); %time is within rows, each column is a stim
% plot(aR_ixR, allBlocksR');
% 
% legend(num2str((1:7)'))
% 
% 
% 
% 
% subplot(2,2,3); hold on
% plot(aL_ixL, mean(allBlocksL), '-k')
% title(sprintf('L (ipsi) antenna only\naverage'), 'FontSize', 12 )
% xlabel('antennal rotation (deg)')
% ylabel('FR (Hz)')
% ax(3) = gca;
% 
% 
% subplot(2,2,4); hold on
% plot(aR_ixR, mean(allBlocksR), '-k')
% xlabel('antennal rotation (deg)')
% ylabel('FR (Hz)')
% title(sprintf('R (contra) antenna only\naverage'), 'FontSize', 12 )
% ax(4) = gca;
% 
% linkaxes(ax, 'y')
% 
% 
% %%
% export_fig(sprintf('AVG_runs1to4_singleAntennae_Vm_%s.pdf', data(1).metadata.experimentHandle.basename))



%% consider only single antenna movements - split in 5x blocks and plot over time --------FR
% figure
% subplot(2,2,1); hold on
% title(sprintf('L (ispi) antenna only\n7 sequential blocks of 5 reps'), 'FontSize', 12 )
% ylabel('FR (Hz)')
% ax(1) = gca;
% 
% 
% ixL = find(aR==0);
% ixL = ixL([2:4,1,5:7]); 
% aL_ixL = aL(ixL);   % xaxis - intervals are not linear, plus values might not be sorted? - but they are
% 
% filtVmZero = [];
% filtVmSt = [];
% for bl = 1:length(data)
%     filtVmZero = cat(1, filtVmZero, data(bl).FR{1});
%     filtVmSt = cat(1, filtVmSt, cat(2,data(bl).FR{setdiff(ixL,1)}) );
% end
% filtVmZero = reshape(filtVmZero, [], 7);
% filtVmZero = squeeze(mean(filtVmZero,1));
% filtVmZero = filtVmZero';
% 
% filtVmSt = permute(filtVmSt, [1,3,2]);
% filtVmSt = reshape(filtVmSt, 5,7,6);
% filtVmSt = squeeze(mean(filtVmSt,1));
% allBlocksL = cat(2, filtVmSt(:,1:3), filtVmZero, filtVmSt(:,4:6)); %time is within rows, each column is a stim
% 
% plot(aL_ixL, allBlocksL');
% 
% 
% subplot(2,2,2); hold on
% title(sprintf('R (contra) antenna only\n7 sequential blocks of 5 reps'), 'FontSize', 12 )
% ax(2) = gca;
% ixR = find(aL==0);
% ixR = ixR([2:4,1,5:7]); 
% aR_ixR = aR(ixR);   % xaxis - intervals are not linear, plus values might not be sorted? - but they are
% 
% filtVmZero = [];
% filtVmSt = [];
% for bl = 1:length(data)
%     filtVmZero = cat(1, filtVmZero, data(bl).FR{1});
%     filtVmSt = cat(1, filtVmSt, cat(2,data(bl).FR{setdiff(ixR,1)}) );
% end
% filtVmZero = reshape(filtVmZero, [], 7);
% filtVmZero = squeeze(mean(filtVmZero,1));
% filtVmZero = filtVmZero';
% 
% filtVmSt = permute(filtVmSt, [1,3,2]);
% filtVmSt = reshape(filtVmSt, 5,7,6);
% filtVmSt = squeeze(mean(filtVmSt,1));
% allBlocksR = cat(2, filtVmSt(:,1:3), filtVmZero, filtVmSt(:,4:6)); %time is within rows, each column is a stim
% 
% plot(aR_ixR, allBlocksR');
% 
% legend(num2str((1:7)'))
% 
% 
% 
% subplot(2,2,3); hold on
% plot(aL_ixL, mean(allBlocksL), '-k')
% title(sprintf('L (ipsi) antenna only\naverage'), 'FontSize', 12 )
% xlabel('antennal rotation (deg)')
% ylabel('FR (Hz)')
% ax(3) = gca;
% 
% 
% subplot(2,2,4); hold on
% plot(aR_ixR, mean(allBlocksR), '-k')
% xlabel('antennal rotation (deg)')
% ylabel('FR (Hz)')
% title(sprintf('R (contra) antenna only\naverage'), 'FontSize', 12 )
% ax(4) = gca;
% 
% linkaxes(ax, 'y')


%%
% export_fig(sprintf('AVG_runs1to4_singleAntennae_FR_%s.pdf', data(1).metadata.experimentHandle.basename))


