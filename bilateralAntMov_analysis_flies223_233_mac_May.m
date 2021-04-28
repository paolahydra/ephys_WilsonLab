% load data for a givn experiment
clear
flyNum     = 223;
cellNum    = 1;
basename = 'fly223_cell01'; %'fly233_cell01_run1'

Nfullblocks2keep = [];
NfirstBlocksToExclude = [];
% basename = 'fly233_cell03'; % temp


armL = 300;
armR = 300;


% [metadataname, cellfolder] = uigetfile(sprintf('/Volumes/CurBioP&W18/project2_2piezo_ephys/p2/fly%3d_PP/fly%3d_cell%02d/metadata*', flyNum, flyNum, cellNum));
[metadataname, cellfolder] = uigetfile(sprintf('/Users/galileo/Dropbox (HMS)/p2/fly%3d_PP/fly%3d_cell%02d/metadata*', flyNum, flyNum, cellNum));

metadata = load(fullfile(cellfolder, metadataname));
flyfolder = metadata.experimentHandle.flyfolder;
basename = metadata.experimentHandle.basename;
flyNum = metadata.experimentHandle.flyNum;      % in case you picked a different one
cellNum = metadata.experimentHandle.cellNum;    % in case you picked a different one

ephysdataname = dir(fullfile(cellfolder, sprintf('ephysData_%s_*.bin', metadata.key.ID)));
ephysdataname = ephysdataname(1).name;



HWsettings = sensor_settings_PP_olderFlies;

fid = fopen(fullfile(cellfolder, ephysdataname), 'r');
data = fread(fid, 'double' );
fclose(fid);

data = reshape(data, 9, []);       % 9 was used before 246   % use 10 after fly 246 (june 1 2018)
timestamps = (data(1,:))';
data = data(2:end,:);
% figure; hold on
% plot(timestamps, data(HWsettings.DATAch.PIEZO_SENSOR_LEFT,:))
% plot(timestamps, data(HWsettings.DATAch.PIEZO_SENSOR_RIGHT,:))

% cd(flyfolder)
ff = strsplit(flyfolder, '\');
fmac = strsplit('Users/galileo/Dropbox (HMS)', '/');
flyfolder = fullfile(fmac{:}, ff{3:4});
flyfolder = cat(2, '/', flyfolder);
cd(flyfolder)

%% convert signals and plot voltage

[ current , voltage, scaled ] = get_scaled_voltage_and_current_PP( data' );
% figure;
% plot(timestamps, scaled)
sensorLfull = data(HWsettings.DATAch.PIEZO_SENSOR_LEFT,:);
sensorRfull = data(HWsettings.DATAch.PIEZO_SENSOR_RIGHT,:);
clear data


%% check onset position
%was this experiment stopped early?
clearvars -except armL armR voltage sensorLfull sensorRfull flyfolder timestamps cellNum flyNum basename metadata metadataname cellfolder Nfullblocks2keep NfirstBlocksToExclude

onsetsUsed = metadata.joint.allOnsetsPositions ( metadata.joint.allOnsetsPositions < length(timestamps) );
trialIndixes = metadata.joint.trialIndices(1:length(onsetsUsed));
repsPerBlock = (length(metadata.dec.joint2L_decoder_i) + metadata.userinput.factor_jointZero - 1) * metadata.userinput.N_withinPseudoBlock;
Nfullblocks = floor(length(onsetsUsed) / repsPerBlock)


if ~isempty(Nfullblocks2keep)
    if Nfullblocks >= Nfullblocks2keep
        Nfullblocks =  Nfullblocks2keep;
    end
end

if ~isempty(NfirstBlocksToExclude)  %introduced on May 19 - 3:30 pm
    indexFirstOnset = NfirstBlocksToExclude*repsPerBlock + 1;
    fprintf('keeping %d full blocks\n', Nfullblocks-NfirstBlocksToExclude)
else
    indexFirstOnset = 1;
    fprintf('keeping %d full blocks\n', Nfullblocks)
end



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
% plot(timestamps, sensorRfull)
% plot(timestamps(onsetsUsed), 5*ones(size(onsetsUsed)), 'xr')

% OK!

%% sample stim traces
maxAbsDisp = max(metadata.dec.i2displ_decoder);
maxTailW = find(metadata.dec.joint2L_decoder_displ == -maxAbsDisp  &  metadata.dec.joint2R_decoder_displ == -maxAbsDisp);
maxContraW = find(metadata.dec.joint2L_decoder_displ == maxAbsDisp  &  metadata.dec.joint2R_decoder_displ == -maxAbsDisp);
maxIpsiW = find(metadata.dec.joint2L_decoder_displ == -maxAbsDisp  &  metadata.dec.joint2R_decoder_displ == maxAbsDisp);
maxHeadW = find(metadata.dec.joint2L_decoder_displ == maxAbsDisp  &  metadata.dec.joint2R_decoder_displ == maxAbsDisp);

STIMS_USE = [maxTailW, maxContraW, metadata.dec.jointZero_index, maxIpsiW, maxHeadW];
STIM_LABS = {'max tail wind', 'max contra wind', 'none', 'max ipsi wind', 'max head wind'};
secs_baseline = metadata.userinput.singleStimDur*1.5;
secs_postStimOnset = metadata.userinput.singleStimDur + metadata.userinput.singleStimDur*1.5;


NCONDITS = length(metadata.dec.i2displ_decoder);


clear perist
VOLTAGE_SR = metadata.userinput.DAQ_fs;

SPIKEDET_STDMULTIPL = 1.25;

delayBetweenSpikePeaksMs = 10; %5; %this really also depends on the threshold...
delayBetweenSpikePeaksPOINTS = delayBetweenSpikePeaksMs * VOLTAGE_SR/1e3;
% nTrialTypes = length(unique(metadata.dec.trialTypes));
nTrialTypes = length(unique(trialIndicesFullBlock));  %patch
durHalfRamp_sec = 0; %there was no ramping in these old experiments

% starting baseline and spike threshold
lastBaseline_I = onsetsUsedFullBlock(1) - durHalfRamp_sec*VOLTAGE_SR;
baselineStart.voltage = voltage(1:lastBaseline_I);
baselineStart.medFiltTrace = medfilt1(baselineStart.voltage, 0.08*VOLTAGE_SR, 'truncate');
baselineStart.delta_Vm = baselineStart.voltage - baselineStart.medFiltTrace; 
baselineStart.timestamps = timestamps(1:lastBaseline_I);

SPIKE_THRESHOLD = std(baselineStart.delta_Vm)*SPIKEDET_STDMULTIPL;
figure; hold on; plot(baselineStart.timestamps, baselineStart.delta_Vm)
plot([baselineStart.timestamps(1), baselineStart.timestamps(end)], [SPIKE_THRESHOLD, SPIKE_THRESHOLD], '--r')
xlabel('time (sec)')
ylabel('mV')
title('high-pass filtered and centered voltage: starting baseline')


%%
volts_thresh = baselineStart.delta_Vm;
volts_thresh(volts_thresh < SPIKE_THRESHOLD) = 0;
[~, baselineStart.locs] = findpeaks(volts_thresh, 'MinPeakDistance',delayBetweenSpikePeaksPOINTS); 
baselineStart.spikes = zeros(1, length(volts_thresh));
baselineStart.spikes(baselineStart.locs) = 1;


for stim = 1:NCONDITS^2
    alltrials = find(trialIndicesFullBlock==stim);
    for i = 1:sum(trialIndicesFullBlock==stim)
        trialN = alltrials(i);
        onset = onsetsUsedFullBlock(trialN);
        firstI_ON = onset + VOLTAGE_SR * durHalfRamp_sec; %stim on, ramp excluded - onset
        lastI_ON = onset + VOLTAGE_SR * (metadata.userinput.singleStimDur - durHalfRamp_sec); %stim on, ramp excluded - offset
        
        perist(stim).voltage(i,:) = voltage(firstI_ON:lastI_ON);
        perist(stim).medFiltTrace(i,:) = medfilt1(perist(stim).voltage(i,:), 0.08*VOLTAGE_SR, 'truncate');
        perist(stim).delta_Vm(i,:) = perist(stim).voltage(i,:) - perist(stim).medFiltTrace(i,:);  %high-pass filtered for spike detection        
        
        
        pre2_ON = onset - VOLTAGE_SR*2*metadata.userinput.singleStimDur;
        post2_ON = onset + VOLTAGE_SR*3*metadata.userinput.singleStimDur;
        
        perist(stim).prepost2voltage(i,:) = nan(length(pre2_ON:post2_ON),1);
        try
            perist(stim).prepost2voltage(i,:) = voltage(pre2_ON:post2_ON);
        catch
            temp = voltage(pre2_ON:end);
            perist(stim).prepost2voltage(i,1:length(temp)) = voltage(pre2_ON:end);
            clear temp
            warning('skipped some points in a prepost2 stim record')
            stim
            i
        end
        prepost2medFiltTrace = medfilt1(perist(stim).prepost2voltage(i,:), 0.08*VOLTAGE_SR, 'truncate');
        prepost2delta_Vm = perist(stim).prepost2voltage(i,:) - prepost2medFiltTrace;  %high-pass filtered for spike detection
            
        if i == 1 && stim == 1
            T_ON = (timestamps(1:length(perist(stim).voltage)));
            prepost2_stimDelimitations = cumsum(repmat(VOLTAGE_SR * [ durHalfRamp_sec, metadata.userinput.singleStimDur-2*durHalfRamp_sec, durHalfRamp_sec], 1, 5 ));
            prepost2_stimDelimitations(3:3:end) = [];
            T_prepost2 = (timestamps(1:length(perist(stim).prepost2voltage)));
        end
        
        volts_thresh = perist(stim).delta_Vm(i,:);
        volts_thresh(volts_thresh < SPIKE_THRESHOLD) = 0;
        [~, locs] = findpeaks(volts_thresh, 'MinPeakDistance',delayBetweenSpikePeaksPOINTS); 
        spikes = zeros(1, length(volts_thresh));
        spikes(locs) = 1;
        perist(stim).locs{i} = locs;
        perist(stim).spikes(i,:) = spikes;
        perist(stim).spikeCount(i) = length(locs);
        perist(stim).sensorL(i,:) = sensorLfull(firstI_ON:lastI_ON);
        perist(stim).sensorR(i,:) = sensorRfull(firstI_ON:lastI_ON);
        
        
        
        volts_thresh = prepost2delta_Vm;
        volts_thresh(volts_thresh < SPIKE_THRESHOLD) = 0;
        [~, locs] = findpeaks(volts_thresh, 'MinPeakDistance',delayBetweenSpikePeaksPOINTS); 
        perist(stim).prepost2SpikeTimes{i,1} = T_prepost2(locs);
        clear prepost2medFiltTrace prepost2delta_Vm
    end
        
% %  old:        
        
%         trialN = alltrials(i);
%         onset = onsetsUsedFullBlock(trialN);
%         firstI = onset - metadata.userinput.DAQ_fs*secs_baseline;
%         lastI = onset + metadata.userinput.DAQ_fs*( secs_postStimOnset);
%         
%         if lastI > length(voltage)
%             perist(stim).voltage(i,:) = nan(1, length(firstI:lastI));
%             perist(stim).voltage(i,1: length(voltage)-firstI+1) = voltage(firstI:length(voltage));
%         else
%             perist(stim).voltage(i,:) = voltage(firstI:lastI);
%         end
% 
%         if i == 1 && stim == 1
%             T = (timestamps(1:length(perist(stim).voltage)) - secs_baseline);
%         end
%         
%         perist(stim).medFiltTrace(i,:) = medfilt1(perist(stim).voltage(i,:), 0.08*VOLTAGE_SR, 'truncate');
%         perist(stim).delta_Vm(i,:) = perist(stim).voltage(i,:) - perist(stim).medFiltTrace(i,:);
%            
%         % psth
%         volts_thresh = perist(stim).delta_Vm(i,:);
%         volts_thresh(volts_thresh < SPIKE_THRESHOLD) = 0;
%         [~, locs] = findpeaks(volts_thresh, 'MinPeakDistance',delayBetweenSpikePeaksPOINTS);
%         spikes = zeros(1, length(volts_thresh));
%         spikes(locs) = 1;
%         perist(stim).locs{i} = locs;
%         perist(stim).spikes(i,:) = spikes;
%         perist(stim).spikeCount(i) = length(locs);
%         perist(stim).sensorL(i,:) = sensorLfull(onset : onset + metadata.userinput.DAQ_fs*metadata.userinput.singleStimDur);
%         perist(stim).sensorR(i,:) = sensorRfull(onset : onset + metadata.userinput.DAQ_fs*metadata.userinput.singleStimDur);
%     end
%     perist(stim).T = T;
end

DCoffset_L = mean(sensorLfull(1:4*VOLTAGE_SR));
DCoffset_R = mean(sensorRfull(1:4*VOLTAGE_SR));
clear sensorL sensorR
for stim = 1:nTrialTypes
    sensorL(stim) = mean(mean(perist(stim).sensorL)) - DCoffset_L;
    sensorR(stim) = mean(mean(perist(stim).sensorR)) - DCoffset_R;
end

%decode angular displancement from arm lengths and sensor data
piezo90_ratio_umperVolt = 9;
dispL_dec = sensorL*piezo90_ratio_umperVolt;
dispR_dec = -sensorR*piezo90_ratio_umperVolt;


angL_dec = atand(dispL_dec/armL);
angR_dec = atand(dispR_dec/armR);

% displacementsL = armL * tand(angL);
% displacementsR = armR * tand(angR);
% piezo90_ratio_umperVolt = 9;

aL = atand(metadata.dec.joint2L_decoder_displ/armL);
aR = atand(metadata.dec.joint2R_decoder_displ/armR);


figure; hold on
plot(angL_dec); 
plot(aL)
legend('decoded', 'command')
title('L angular displacements')
ylabel('deg')

figure; hold on
plot(angR_dec); 
plot(aR)
legend('decoded', 'command')
title('R angular displacements')
ylabel('deg')

save(sprintf('%s_%s_periStimulusData.mat', basename, metadata.key.ID ), 'perist')

%% try new 
% plot all example traces (early, middle, late) per sample stimulus
Ntr2plot = 5;

figure('WindowStyle', 'normal', 'Position', [1921 41 1024 1163])

[ha, pos] = tight_subplot(length(STIMS_USE), Ntr2plot, [0.025, 0.05]);
YLIM = [-40, 10];  %[-50, 10];

ha_idx = 1;
for stim = STIMS_USE
    Ntr2plot = 5;
    trNumbers = round(linspace(1, sum(trialIndicesFullBlock==1), Ntr2plot));
    for i_tr = 1 : Ntr2plot
        trN = trNumbers(i_tr);
        axes(ha(ha_idx));
        plot(T_ON, perist(stim).voltage(trN,:))
        hold on
        plot(T_ON(perist(stim).locs{trN}), -28*ones(size(perist(stim).locs{trN})), 'xg')
        xlim([0, metadata.userinput.singleStimDur])
        if stim == STIMS_USE(1)
            if i_tr == 1
                title(sprintf('%s %s\nblock %d', metadata.key.ID(1:8), metadata.key.ID(10:end), i_tr))
            else
                title(sprintf('%d',i_tr))
            end
        end
        
        ylim(YLIM)
        ha(ha_idx).Box = 'Off';
        ha(ha_idx).YGrid = 'on';
        if ~ismember(ha_idx, 1:Ntr2plot:length(ha))
            ha(ha_idx).YAxis.Visible = 'off';
        else
            lab = STIM_LABS( ha_idx == 1:Ntr2plot:length(ha));
            ylabel(sprintf('%s (mV)', lab{1}))
        end
        if ~ismember(ha_idx, length(ha)-Ntr2plot+1:length(ha))
            ha(ha_idx).XAxis.Visible = 'off';
        else
            xlabel(sprintf('time (s)'))
        end
        ha_idx = ha_idx+1;
    end
end

savefig(sprintf('sampleTrials_%s_voltageTrace_%s', basename, metadata.key.ID))
export_fig(sprintf('sampleTrials_%s_voltageTrace_%s.pdf', basename, metadata.key.ID))

%% new
for stim = 1:nTrialTypes
    filtVm{stim} = mean(perist(stim).medFiltTrace, 2);
    FR{stim} = sum(perist(stim).spikes,2);
    
    meanfiltVm(stim) =  mean(mean(perist(stim).medFiltTrace));
    meanFR(stim) = mean(sum(perist(stim).spikes,2));
end


%% matrix FR - ok
matrixFR = nan(NCONDITS);
for R = 1:NCONDITS % by rows
    for L = 1:NCONDITS % by cols
        stimN = find(metadata.dec.joint2L_decoder_i==L & metadata.dec.joint2R_decoder_i==R);
        matrixFR(R,L) = meanFR(stimN);
    end
end
zeroFR = matrixFR(ceil(NCONDITS/2), ceil(NCONDITS/2));
matrixFR_BS = matrixFR - zeroFR; %baseline subtracted


%%
figure;
imagesc(matrixFR_BS)
ax = gca;
colormap(bluewhitered(256)), cb = colorbar;
cb.TickLabels = num2str((cb.Ticks+zeroFR)', '%2.1f')
cb.Label.String = {'firing rate (Hz)'; 'averaged during stimulus ON'};
axis square
% ax.XAxisLocation = 'top';
ax.YDir = 'normal';
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

export_fig(sprintf('%s_%s_matrix_FR_mac.pdf', basename, metadata.key.ID ))
export_fig(sprintf('%s_%s_matrix_FR_mac.jpg', basename, metadata.key.ID ))



%% 1 estimate each arm length before starting analysis, and input it in first block

%% to do: 2 save data with same format as regular analysis so that it can be reaccessed and combined across runs
directions = [];
intensities = [];
save(fullfile(flyfolder, sprintf('data_%s_%s.mat',basename, metadata.key.ID)), ...
    'meanfiltVm', 'meanFR', 'filtVm', 'FR', 'metadata', 'perist', 'directions', 'intensities', 'aL', 'aR', 'trialIndicesFullBlock', 'prepost2_stimDelimitations', 'T_prepost2' )

%% to do: 3 plot data points in angular displacement space

figure; hold on; axis image
xlabel(' Left antenna (ipsi) - angular displ')
ylabel(' Right antenna (contra) - angular displ')
scatter(aL, aR, 200, meanfiltVm-meanfiltVm(metadata.dec.jointZero_index), 'filled', 'MarkerEdgeColor', 'k')
title(sprintf('Vm\n%s\n%s', basename, metadata.key.ID), 'Interpreter', 'none')
colormap(bluewhitered(256)), cb = colorbar;
cb.TickLabels = num2str((cb.Ticks+meanfiltVm(metadata.dec.jointZero_index))', '%2.1f');
cb.Label.String = '(mV)';
xlim([-15 15])
ylim([-15 15])
set(gca, 'TickDir', 'out')
export_fig(sprintf('scatter_Vm_%s_%s_mac.pdf', basename, metadata.key.ID))


figure; hold on; axis image
xlabel(' Left antenna (ipsi) - angular displ')
ylabel(' Right antenna (contra) - angular displ')
scatter(aL, aR, 200, meanFR-meanFR(metadata.dec.jointZero_index), 'filled', 'MarkerEdgeColor', 'k')
title(sprintf('firing rate\n%s\n%s', basename, metadata.key.ID), 'Interpreter', 'none')
colormap(bluewhitered(256)), cb = colorbar;
cb.TickLabels = num2str((cb.Ticks+meanFR(metadata.dec.jointZero_index))', '%2.1f');
cb.Label.String = '(Hz)';
xlim([-15 15])
ylim([-15 15])
set(gca, 'TickDir', 'out')
export_fig(sprintf('scatter_FR_%s_%s_mac.pdf', basename, metadata.key.ID))

%% to do: 4 plot explicit single antenna contributions

%% to do: 5 find all good trials and all good cells and repeat!!!!!!



