% load data for a givn experiment
flyNum     = 248;
cellNum    = 1;
% basename = 'fly233_cell03'; % temp


SPIKEDET_STDMULTIPL = 2;
YLIM = [-50, 20];

[metadataname, cellfolder] = uigetfile(sprintf('D:/Dropbox (HMS)/p2/fly%3d_PP/fly%3d_cell%02d/metadata*', flyNum, flyNum, cellNum));
metadata = load(fullfile(cellfolder, metadataname));
flyfolder = metadata.experimentHandle.flyfolder;
basename = sprintf('fly%d_cell%02d',flyNum,cellNum);


ephysdataname = dir(fullfile(cellfolder, sprintf('ephysData_%s_*.bin', metadata.key.ID)));
ephysdataname = ephysdataname.name;


HWsettings = sensor_settings_PP;

fid = fopen(fullfile(cellfolder, ephysdataname), 'r');
data = fread(fid, 'double' );
fclose(fid);

data = reshape(data, 10, []);        % use 10 after fly 246 (june 1 2018)
timestamps = (data(1,:))';
data = data(2:end,:);
% figure; hold on
% plot(timestamps, data(HWsettings.DATAch.PIEZO_SENSOR_LEFT,:))
% plot(timestamps, data(HWsettings.DATAch.PIEZO_SENSOR_RIGHT,:))

cd(flyfolder)

%% convert signals and plot voltage

[ current , voltage, scaled ] = get_scaled_voltage_and_current_PP( data' );
% figure;
% plot(timestamps, scaled)

% 
% %%
% XLIM = [100, 150];
% 
% figure;
% plot(timestamps, voltage)
% title(sprintf('%s - %s', basename, metadata.key.ID), 'Interpreter', 'none')
% ylabel('Voltage (mV)')
% xlabel('time (s)')
% ylim([-65, 10])
% xlim(XLIM)
% box off
% set(gca, 'TickDir', 'out')
% axV = gca;
% 
% figure;
% plot(timestamps, current)
% axC = gca;
% ylabel('Current (pA)')
% xlabel('time (s)')
% ylim([-100, 100])
% xlim(XLIM)
% box off
% set(gca, 'TickDir', 'out')
% 
% % linkaxes([axV, axC], 'x')

%% check onset position
%was this experiment stopped early?

onsetsUsed = metadata.joint.allOnsetsPositions ( metadata.joint.allOnsetsPositions < length(timestamps) );
trialIndixes = metadata.joint.trialIndices(1:length(onsetsUsed));
repsPerBlock = (length(metadata.dec.joint2L_decoder_i) + metadata.userinput.factor_jointZero - 1) * metadata.userinput.N_withinPseudoBlock;
Nfullblocks = floor(length(onsetsUsed) / repsPerBlock);
onsetsUsedFullBlock = onsetsUsed(1: Nfullblocks*repsPerBlock);
trialIndicesFullBlock = metadata.joint.trialIndices(1:length(onsetsUsedFullBlock));

    
if length(metadata.joint.allOnsetsPositions) > length(onsetsUsed)
    STOPPEDeARLY = 1;
elseif timestamps(end) < timestamps(onsetsUsed(end)) + metadata.userinput.singleStimDur + 3
    STOPPEDeARLY = 1;
else
    STOPPEDeARLY = 0;
end


% figure; hold on
% plot(timestamps, data(HWsettings.DATAch.PIEZO_SENSOR_RIGHT,:))
% plot(timestamps(onsetsUsed), 5*ones(size(onsetsUsed)), 'xr')

% OK!

%% sample stim traces
maxAbsDisp = max(metadata.dec.i2displ_decoder);
maxTailW = find(metadata.dec.joint2L_decoder_displ == -maxAbsDisp  &  metadata.dec.joint2R_decoder_displ == -maxAbsDisp);
maxContraW = find(metadata.dec.joint2L_decoder_displ == maxAbsDisp  &  metadata.dec.joint2R_decoder_displ == -maxAbsDisp);
maxIpsiW = find(metadata.dec.joint2L_decoder_displ == -maxAbsDisp  &  metadata.dec.joint2R_decoder_displ == maxAbsDisp);
maxHeadW = find(metadata.dec.joint2L_decoder_displ == maxAbsDisp  &  metadata.dec.joint2R_decoder_displ == maxAbsDisp);


VOLTAGE_SR = metadata.userinput.DAQ_fs;
NCONDITS = length(metadata.dec.i2displ_decoder);





STIMS_USE = [maxTailW, maxContraW, metadata.dec.jointZero_index, maxIpsiW, maxHeadW];
STIM_LABS = {'max tail wind', 'max contra wind', 'none', 'max ipsi wind', 'max head wind'};
secs_baseline = metadata.userinput.singleStimDur*1.5;
secs_postStimOnset = metadata.userinput.singleStimDur + metadata.userinput.singleStimDur*1.5;

tic
for stim = 1:NCONDITS^2
    alltrials = find(trialIndicesFullBlock==stim);
    for i = 1:sum(trialIndicesFullBlock==stim)
        trialN = alltrials(i);
        onset = onsetsUsedFullBlock(trialN);
        firstI = onset - metadata.userinput.DAQ_fs*secs_baseline;
        lastI = onset + metadata.userinput.DAQ_fs*( secs_postStimOnset);
        if lastI > length(voltage)
            perist(stim).voltage(i,:) = nan(1, length(firstI:lastI));
            perist(stim).voltage(i,1: length(voltage)-firstI+1) = voltage(firstI:length(voltage));
        else
            perist(stim).voltage(i,:) = voltage(firstI:lastI);
        end

        if i == 1 && stim == 1
            T = (timestamps(1:length(perist(stim).voltage)) - secs_baseline);
        end
        
        perist(stim).medFiltTrace(i,:) = medfilt1(perist(stim).voltage(i,:), 0.08*VOLTAGE_SR, 'truncate');
        perist(stim).delta_Vm(i,:) = perist(stim).voltage(i,:) - perist(stim).medFiltTrace(i,:);
        
        if i == 1 && stim == 1
            
            SPIKE_THRESHOLD = std(perist(stim).delta_Vm(i,:))*SPIKEDET_STDMULTIPL;
            
            figure; hold on
            plot(T, perist(stim).voltage(i,:))
            plot(T, perist(stim).medFiltTrace(i,:))
            xlabel('peristimulus time (sec)')
            ylabel('mV')
            figure, hold on
            plot(T, perist(stim).delta_Vm(i,:))
            plot([T(1), T(end)], [SPIKE_THRESHOLD, SPIKE_THRESHOLD], '--r')
            xlabel('peristimulus time (sec)')
            ylabel('mV')
        end    
        % psth
        volts_thresh = perist(stim).delta_Vm(i,:);
        volts_thresh(find(volts_thresh < SPIKE_THRESHOLD)) = 0;
        [~, locs] = findpeaks(volts_thresh, 'MinPeakDistance',50); %40 should be 1 ms
        spikes = zeros(1, length(volts_thresh));
        spikes(locs) = 1;
        if i == 1 && stim == 1
            plot(T(locs), SPIKE_THRESHOLD*ones(size(locs)), 'or')
        end
        perist(stim).spikes(i,:) = spikes;
    end
    perist(stim).T = T;
end
toc

save(sprintf('%s_%s_periStimulusData.mat', basename, metadata.key.ID ), 'perist')

%% for each stimulus, plot 5 single trials (sampled over time), mean Vm, psth - SAVE ALL


Ntr2plot = 5;
trNumbers = round(linspace(1, sum(trialIndicesFullBlock==1), Ntr2plot));
for stim = 1:NCONDITS^2
    figure('Color', [1,1,1])
    for i_tr = 1 : Ntr2plot
        trN = trNumbers(i_tr);
    	subplot(Ntr2plot+2, 1, i_tr)
        plot(T, perist(stim).voltage(trN,:))
        ylim(YLIM)
        axis off
        hold on
        plot([0,1,1,0,0], [YLIM(1) YLIM(1) YLIM(2) YLIM(2) YLIM(1)], 'Color', [0.67 0.67 0.67])
        if i_tr == 1
            plot([0 0], [YLIM(2) YLIM(2) - 10], '-k')
            plot([0,0.1], [YLIM(2) YLIM(2)], '-k')
        end   
    end
    
    subplot(Ntr2plot+2, 1, Ntr2plot+1)
    meanVM = mean(perist(stim).medFiltTrace);
    plot(T, meanVM, 'r')
    ylim(YLIM)
    axis off
    hold on
    plot([0,1,1,0,0], [YLIM(1) YLIM(1) YLIM(2) YLIM(2) YLIM(1)], 'Color', [0.67 0.67 0.67])
    
    subplot(Ntr2plot+2, 1, Ntr2plot+2)
    plotSpikeRaster(logical(perist(stim).spikes),'PlotType','vertline', 'AxHandle',gca);
    axis off
    
    export_fig(sprintf('%s_%s_stim%02d.pdf', basename, metadata.key.ID, stim ))
    close
end

%% psth and mean firing rate 
psth_dt = 0.1; %don't change. it will screw hard numbers (4:5) few lines later
BIN_SIZE = psth_dt*VOLTAGE_SR; %PSTH
time_psth = (0:psth_dt:length(perist(1).T)/VOLTAGE_SR) + perist(1).T(1);

meanFR = [];
for stim = 1:NCONDITS^2
% stim = 7;
    psth_discr = mean(perist(stim).spikes);
    psth_discr = psth_discr';
    psth_discr = psth_discr(1 : BIN_SIZE*floor(length(psth_discr)/BIN_SIZE));
    psth_discr = reshape(psth_discr, BIN_SIZE, []);
    psth = sum(psth_discr)/psth_dt;
%     figure; bar(time_psth(1:end-1), psth)
    meanFR(stim) = mean(psth(17:25)); %excluding firts 100 ms [0.1,1]
end

%% matrix FR
matrixFR = nan(NCONDITS);
for R = 1:NCONDITS % by rows
    for L = 1:NCONDITS % by cols
        stimN = find(metadata.dec.joint2L_decoder_i==L & metadata.dec.joint2R_decoder_i==R);
        matrixFR(R,L) = meanFR(stimN);
    end
end
zeroFR = matrixFR(ceil(NCONDITS/2), ceil(NCONDITS/2));
matrixFR_BS = matrixFR - zeroFR; %baseline subtracted

figure;
imagesc(matrixFR_BS)
ax = gca;
colormap(bluewhitered(256)), cb = colorbar;
cb.TickLabels = num2str((cb.Ticks+zeroFR)', '%2.1f')
cb.Label.String = {'firing rate (Hz)'; 'averaged during stimulus ON'};
axis square
ax.XAxisLocation = 'top';
ax.XTick = 1:NCONDITS;
ax.YTick = 1:NCONDITS;

ax.XTickLabel = metadata.dec.i2displ_decoder;
ax.YTickLabel = metadata.dec.i2displ_decoder;

L = ax.YTickLabel;
LY = cat(2, repmat('     ', NCONDITS,1), L);
LY(1,1:4) = 'pull';
LY(end,1:4) = 'push';
ax.YTickLabel = LY;
ax.YTickLabel(1,:) = cat(2, 'pull ',  L(1,:));

xlabel('IPSI displacement (um)')
ylabel('CONTRA displacement (um)')

export_fig(sprintf('%s_%s_matrix_FR.pdf', basename, metadata.key.ID ))


%% matrix filtered Vm
T_stON = perist(1).T > 0.1 & perist(1).T<1;
matrixfiltVm = nan(NCONDITS);
for R = 1:NCONDITS % by rows
    for L = 1:NCONDITS % by cols
        stimN = find(metadata.dec.joint2L_decoder_i==L & metadata.dec.joint2R_decoder_i==R);
        matrixfiltVm(R,L) = mean(mean(perist(stimN).medFiltTrace(:,T_stON)));
    end
end
zerofiltVm = matrixfiltVm(ceil(NCONDITS/2), ceil(NCONDITS/2));
matrixfiltVm_BS = matrixfiltVm - zerofiltVm; %baseline subtracted


figure;
imagesc(matrixfiltVm_BS)
ax = gca;
colormap(bluewhitered(256)), cb = colorbar;
cb.Label.String = 'average mV change from [0,0] displacement';

axis square
ax.XAxisLocation = 'top';
ax.XTick = 1:NCONDITS;
ax.YTick = 1:NCONDITS;
ax.XTickLabel = metadata.dec.i2displ_decoder;
ax.YTickLabel = metadata.dec.i2displ_decoder;
L = ax.YTickLabel;
LY = cat(2, repmat('     ', NCONDITS,1), L);
LY(1,1:4) = 'pull';
LY(end,1:4) = 'push';
ax.YTickLabel = LY;
ax.YTickLabel(1,:) = cat(2, 'pull ',  L(1,:));
xlabel('IPSI displacement (um)')
ylabel('CONTRA displacement (um)')

export_fig(sprintf('%s_%s_matrix_filteredVm.pdf', basename, metadata.key.ID ))

% %% matrix -old way
% NCONDITS = length(metadata.dec.i2displ_decoder);
% matrix = nan(NCONDITS);
% 
% % pre-build - stim[onset,offset] window
% for stim = 1:NCONDITS^2
%     alltrials = find(trialIndicesFullBlock==stim);
%     clear periV
%     for i = 1:sum(trialIndicesFullBlock==stim)
%         trialN = alltrials(i);
%         onset = onsetsUsedFullBlock(trialN);
%         firstI = onset;
%         lastI = onset + metadata.userinput.DAQ_fs*(metadata.userinput.singleStimDur);
%         periV(i,:) = voltage(firstI:lastI);
%     end
%     mPeriV_joint(stim,:) = mean(periV);
% end
% 
% % build matrix
% for R = 1:NCONDITS % by rows
%     for L = 1:NCONDITS % by cols
%         stimN = find(metadata.dec.joint2L_decoder_i==L & metadata.dec.joint2R_decoder_i==R);
%         matrix(R,L) = mean(mPeriV_joint(stimN,:));
%     end
% end
% matrixBS = matrix - matrix(ceil(NCONDITS/2), ceil(NCONDITS/2)); %baseline subtracted
% 
% 
% figure;
% imagesc(matrixBS)
% ax = gca;
% colormap(bluewhitered(256)), cb = colorbar;
% cb.Label.String = 'average mV change from [0,0] displacement';
% axis square
% ax.XAxisLocation = 'top';
% ax.XTick = 1:NCONDITS;
% ax.YTick = 1:NCONDITS;
% 
% ax.XTickLabel = metadata.dec.i2displ_decoder;
% 
% ax.YTickLabel = metadata.dec.i2displ_decoder;
% L = ax.YTickLabel;
% LY = cat(2, repmat('     ', NCONDITS,1), L);
% LY(1,1:4) = 'pull';
% LY(end,1:4) = 'push';
% ax.YTickLabel = LY;
% ax.YTickLabel(1,:) = cat(2, 'pull ',  L(1,:));
% 
% xlabel('IPSI displacement (um)')
% ylabel('CONTRA displacement (um)')
% export_fig(sprintf('%s_%s_matrix.pdf', basename, metadata.key.ID ))
%% show antennae
figure;
fr = 1;
subplot(1,3,1)
imshowpair(allFr(2).frames(:,:,fr), allFr(2).frames(:,:,2))
subplot(1,3,2)
imshowpair(allFr(1).frames(:,:,fr), allFr(1).frames(:,:,2))
subplot(1,3,3)
imshowpair(allFr(3).frames(:,:,fr), allFr(3).frames(:,:,2))


%% botttom camera, piezo reaching sequentially

flyNum = 235;
ID = '20180525T184936'; % 20180609T171444


folder = sprintf('D:\\Dropbox (HMS)\\p2\\fly%d_PP',flyNum);
cd(folder)
load(sprintf('piezosDispl_fly%d_%s.mat',flyNum, ID))

for cam = 1:3
    for fr = 1:9
        
        frameim = allFr(cam).frames(:,:,fr);
        frameim = imadjust(frameim);

        imwrite(frameim, sprintf('cam%d_fr%02d_fly%d_%s.tiff', cam, fr, flyNum, ID))
    end
end



%% read in command
fidC = fopen(fullfile(cellfolder, metadata.key.logname), 'r');
command = fread(fidC, 'single');
fclose(fidC);
command = reshape(command, 2, []);


% %% left piezo sensor and command
% figure; hold on
% DC_left = mean(data(HWsettings.DATAch.PIEZO_SENSOR_LEFT,1:20e4));
% plot(timestamps, data(HWsettings.DATAch.PIEZO_SENSOR_LEFT,:))
% plot(timestamps, DC_left + command(1,1:length(timestamps)))
% xlim([20 60])
% title(sprintf('left piezo sensor and command - %s', basename), 'Interpreter', 'none')
% 
% %% right piezo sensor and command
% figure; hold on
% DC_right = mean(data(HWsettings.DATAch.PIEZO_SENSOR_RIGHT,1:20e4));
% plot(timestamps, data(HWsettings.DATAch.PIEZO_SENSOR_RIGHT,:))
% plot(timestamps, DC_right + command(2,1:length(timestamps)))
% xlim([20 60])
% title(sprintf('right piezo sensor and command - %s', basename), 'Interpreter', 'none')
% 
% %%
%  tic
% 
%     fprintf('\nRECORDING: elapsed time (s): %8.2f\t',toc)
%     % start the loop
%     i = 0;
%     while i < 100       % Check if the loop has to be stopped
%         pause(0.05)
%         i = i+1;
%         fprintf('%c',repmat(8,9,1)) ;   % clear up previous time
%         fprintf('%8.2f\t',toc) ;        % display elapsed time
%     end
    