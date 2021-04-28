% load data for a givn experiment
clear
flyNum     = 223;
cellNum    = 1;
Nfullblocks2keep = [];
NfirstBlocksToExclude = [];
% basename = 'fly233_cell03'; % temp


% [metadataname, cellfolder] = uigetfile(sprintf('/Volumes/CurBioP&W18/project2_2piezo_ephys/p2/fly%3d_PP/fly%3d_cell%02d/metadata*', flyNum, flyNum, cellNum));
[metadataname, cellfolder] = uigetfile(sprintf('/Users/galileo/Dropbox (HMS)/p2/fly%3d_PP/fly%3d_cell%02d/metadata*', flyNum, flyNum, cellNum));

metadata = load(fullfile(cellfolder, metadataname));
flyfolder = metadata.experimentHandle.flyfolder;
basename = metadata.experimentHandle.basename;
flyNum = metadata.experimentHandle.flyNum;      % in case you picked a different one
cellNum = metadata.experimentHandle.cellNum;    % in case you picked a different one

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

clearvars -except voltage sensorLfull sensorRfull flyfolder timestamps cellNum flyNum basename metadata metadataname cellfolder Nfullblocks2keep NfirstBlocksToExclude




onsetsUsed = metadata.joint.allOnsetsPositions ( metadata.joint.allOnsetsPositions < length(timestamps) );
trialIndixes = metadata.joint.trialIndices(1:length(onsetsUsed));
repsPerBlock = length(metadata.dec.trialTypes) * metadata.userinput.N_withinPseudoBlock;
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


onsetsUsedFullBlock = onsetsUsed(indexFirstOnset: Nfullblocks*repsPerBlock);
trialIndicesFullBlock = metadata.joint.trialIndices(indexFirstOnset:length(onsetsUsedFullBlock));

    
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


% figure;
% ts_plot = timestamps>=0 & timestamps<100;
% plot(timestamps(ts_plot), sensorRfull(ts_plot))


% from these onsets I still want to discrd 150ms of ramping

%% PATCH FOR DETECTED ERROR IN STIM GENERATION - THERE MIGHT BE MORE
trialIndicesFullBlock(trialIndicesFullBlock==2) = 1;
trialIndicesFullBlock(trialIndicesFullBlock==3) = 1;
trialIndicesFullBlock(trialIndicesFullBlock==4) = 1;

trialIndicesFullBlock(trialIndicesFullBlock>4) = trialIndicesFullBlock(trialIndicesFullBlock>4)-3;

assert(length(unique(trialIndicesFullBlock)) == 84);


%% no matrix no more : color-coded scatter!
% 0. did I not filter spikes out?
    %         trialN = alltrials(i);
    %         onset = onsetsUsedFullBlock(trialN);
    %         firstI = onset - metadata.userinput.DAQ_fs*secs_baseline;
    %         lastI = onset + metadata.userinput.DAQ_fs*( secs_postStimOnset);
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
        
        
% 1. sort all the stimuli with unique trialType number.
% 2. extract trace from onset +0.150s to offset-0.150s
% 3. filter
clear perist
VOLTAGE_SR = metadata.userinput.DAQ_fs;
SPIKEDET_STDMULTIPL = 2;
delayBetweenSpikePeaksMs = 8; %5; %this really also depends on the threshold...


delayBetweenSpikePeaksPOINTS = delayBetweenSpikePeaksMs * VOLTAGE_SR/1e3;
% nTrialTypes = length(unique(metadata.dec.trialTypes));
nTrialTypes = length(unique(trialIndicesFullBlock));  %patch

durHalfRamp_sec = 0.150;

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


% check threshold now


%%
volts_thresh = baselineStart.delta_Vm;
volts_thresh(volts_thresh < SPIKE_THRESHOLD) = 0;
[~, baselineStart.locs] = findpeaks(volts_thresh, 'MinPeakDistance',delayBetweenSpikePeaksPOINTS); 
baselineStart.spikes = zeros(1, length(volts_thresh));
baselineStart.spikes(baselineStart.locs) = 1;
%OK!
              
        
for stim = 1:nTrialTypes
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

angL_dec = atand(dispL_dec/metadata.dec.armL);
angR_dec = atand(dispR_dec/metadata.dec.armR);

% displacementsL = armL * tand(angL);
% displacementsR = armR * tand(angR);
% piezo90_ratio_umperVolt = 9;

aL = metadata.dec.angL(4:end-3);
aR = metadata.dec.angR(4:end-3);

figure; hold on
plot(angL_dec); 
plot(aL)
legend('decoded', 'command')
title('L angular displacements')

figure; hold on
plot(angR_dec); 
plot(aR)
legend('decoded', 'command')
title('R angular displacements')


%% regress resting position response levels to be subtracted                        -- CHECK THIS!!!! seems way too noisy - do I actually use it????
trialListI = 1:length(trialIndicesFullBlock);
resting.PosI = find(trialIndicesFullBlock==metadata.dec.jointZero_index);
blockSize = length(metadata.dec.trialTypes);
% figure; hold on
% % plot(resting.PosI, zeros(size(resting.PosI)), '+r')
% xticks(0.5:blockSize:length(trialIndicesFullBlock)+0.5)
% ax = gca;
% ax.XGrid = 'on';
% ax.YGrid = 'off';
% ax.XTickLabel = [];
% % ax.YAxis.Visible = 'off';

% % LATER ADD PRE AND POST BASELINE AVERAGE POINTS TOO
% 
% 
% stims120.stimI = find(metadata.dec.intensities == 120);
% for i = 1:length(stims120.stimI)
%     stims120.medFiltTrace_avg(:,i) = mean(perist(stims120.stimI(i)).medFiltTrace, 2);
%     stims120.PosI(:,i) = find(trialIndicesFullBlock == stims120.stimI(i));
%     plot(stims120.PosI(:,i), stims120.medFiltTrace_avg(:,i), 'x-')
%     pause(0.25)
% end
% 
% resting.medFiltTrace_avg = mean(perist(metadata.dec.jointZero_index).medFiltTrace, 2);
% plot(resting.PosI, resting.medFiltTrace_avg, '+-r')



%% SAMPLE STIMULI
%patching:
% basically the scatter plot heere is correct, because angular
% displacements have been shifted and indeed they match the decoded ones
% all the other labels have not been shifted though - so they need extra
% care. e.g.:
% % aL = metadata.dec.angL(4:end-3);
% % aR = metadata.dec.angR(4:end-3);
directions = metadata.dec.directions(4:end-3);
intensities = metadata.dec.intensities(4:end-3);

max_leftIpsi_n70deg_120amp = find(directions(1:71) == -70 & intensities(1:71) == 120);
max_rightContra_70deg_120amp = find(directions(1:71) == 70 & intensities(1:71) == 120);
max_Head_0deg_120amp = find(directions(1:71) == 0 & intensities(1:71) == 120);
max_Tail_360deg_120amp = 84;

STIMS_USE = [ metadata.dec.jointZero_index, max_Tail_360deg_120amp, max_rightContra_70deg_120amp, max_Head_0deg_120amp, max_leftIpsi_n70deg_120amp];
STIM_LABS = {'none', 'max tail wind', 'max contra wind', 'max head wind' , 'max ipsi wind'};
% (ok)


%% plot all example traces (early, middle, late) per sample stimulus
Ntr2plot = Nfullblocks;

figure('WindowStyle', 'normal', 'Position', [1921 41 1024 1163])

[ha, pos] = tight_subplot(length(STIMS_USE), Ntr2plot, [0.025, 0.05]);
YLIM = [-70, -20];  %[-50, 10];

ha_idx = 1;
for stim = STIMS_USE
    trNumbers = round(linspace(1, sum(trialIndicesFullBlock==stim), Ntr2plot));
    blNumbers = round(linspace(1, Nfullblocks, Ntr2plot));
    for i_tr = 1 : Ntr2plot
        trN = trNumbers(i_tr);
        blN = blNumbers(i_tr);
        axes(ha(ha_idx));
        plot(T_ON, perist(stim).voltage(trN,:))
        hold on
        plot(T_ON(perist(stim).locs{trN}), -28*ones(size(perist(stim).locs{trN})), 'xg')
        
        if stim == STIMS_USE(1)
            if i_tr == 1
                title(sprintf('%s %s\nblock %d', metadata.key.ID(1:8), metadata.key.ID(10:end), blN))
            else
                title(sprintf('%d/%d',blN, Nfullblocks))
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

savefig(sprintf('sampleTrials_voltageTrace_%s_all_mac', metadata.key.ID))
export_fig(sprintf('sampleTrials_voltageTrace_%s_all_mac.pdf', metadata.key.ID))

%% % plot all traces with two preceding and two following stimuli, arranged vertically - skip - it takes too long now. Select example data at the end
figFolder = fullfile(flyfolder, sprintf('singleStimsTraces_%s_%s',basename, metadata.key.ID));
mkdir(figFolder);

for stim = 1:length(perist)
    
    figure('WindowStyle', 'normal', 'Position', [916     1   523   804])
    Ntr2plot = Nfullblocks;
    [ha, pos] = tight_subplot( Ntr2plot, 1, [0.025, 0.25]);
    YLIM = [-65, -15];
    ha_idx = 1;

    trNumbers = round(linspace(1, sum(trialIndicesFullBlock==stim), Ntr2plot));
    blNumbers = round(linspace(1, Nfullblocks, Ntr2plot));
    for i_tr = 1 : Ntr2plot
        trN = trNumbers(i_tr);
        blN = blNumbers(i_tr);

        axes(ha(ha_idx));
        plot(T_prepost2, perist(stim).prepost2voltage(trN,:))
        hold on
        
        axis tight
        ha(ha_idx).Box = 'Off';
        ylim(YLIM)
        ha(ha_idx).XTick = T_prepost2(prepost2_stimDelimitations);
        if ha_idx==1
            ha(ha_idx).YAxis.Visible = 'on';
            stringName = sprintf('%s_%s_stim%02d', metadata.key.ID(1:8), metadata.key.ID(10:end), stim);
            title(stringName, 'Interpreter', 'none')
        else
            ha(ha_idx).YAxis.Visible = 'off';
        end
        ha(ha_idx).XAxis.Visible = 'off';
        ha(ha_idx).XGrid = 'on';
        ha(ha_idx).YGrid = 'on';
        
        ha_idx = ha_idx+1;
    end
    
    export_fig(fullfile(figFolder, sprintf('%s.eps',stringName)))
    close

end


%%
Ntr2plot = 2;

figure('WindowStyle', 'normal', 'Position', [1921 41 1024 1163])

[ha, pos] = tight_subplot(length(STIMS_USE), Ntr2plot, [0.025, 0.05]);
YLIM = [-55, -5];

ha_idx = 1;
for stim = STIMS_USE
    trNumbers = round(linspace(1, sum(trialIndicesFullBlock==stim), Ntr2plot));
    blNumbers = round(linspace(1, Nfullblocks, Ntr2plot));
    for i_tr = 1 : Ntr2plot
        trN = trNumbers(i_tr);
        blN = blNumbers(i_tr);
        axes(ha(ha_idx));
        plot(T_ON, perist(stim).voltage(trN,:))
        hold on
        plot(T_ON(perist(stim).locs{trN}), -28*ones(size(perist(stim).locs{trN})), 'xg')
        
        if stim == STIMS_USE(1)
            if i_tr == 1
                title(sprintf('%s %s\nblock %d', metadata.key.ID(1:8), metadata.key.ID(10:end), blN))
            else
                title(sprintf('%d/%d',blN, Nfullblocks))
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

savefig(sprintf('sampleTrials_voltageTrace_%s_firstLast_mac', metadata.key.ID))
export_fig(sprintf('sampleTrials_voltageTrace_%s_firstLast_mac.pdf', metadata.key.ID))

%% direction tuning curves
for stim = 1:nTrialTypes
    filtVm{stim} = mean(perist(stim).medFiltTrace, 2);
    FR{stim} = sum(perist(stim).spikes,2);
    
    meanfiltVm(stim) =  mean(mean(perist(stim).medFiltTrace));
    meanFR(stim) = mean(sum(perist(stim).spikes,2));
end
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
title({basename; metadata.key.ID}, 'Interpreter', 'none')
f.WindowStyle = 'normal';
f.Position = [1922 516 1019 364];
savefig(sprintf('dirTuning_filteredVoltageRestSubtr_%s', metadata.key.ID))
export_fig(sprintf('dirTuning_filteredVoltageRestSubtr_%s.pdf', metadata.key.ID))


f = figure; hold on
plot(directions(I30), meanFR(I30)-meanFR(1)); 
plot(directions(I60), meanFR(I60)-meanFR(1));
plot(directions(I120), meanFR(I120)-meanFR(1));
legend('30', '60', '120')
legend boxoff
legend('Location', 'NorthWest')
xlabel('estimated wind direction')
ylabel('mean firing rate change with respect to rest')
title({basename; metadata.key.ID}, 'Interpreter', 'none')
f.WindowStyle = 'normal';
f.Position = [1922 516 1019 364];

savefig(sprintf('dirTuning_firingRateRestSubtr_%s_mac', metadata.key.ID))
export_fig(sprintf('dirTuning_firingRateRestSubtr_%s_mac.pdf', metadata.key.ID))

save(fullfile(flyfolder, sprintf('data_%s_%s.mat',metadata.experimentHandle.basename, metadata.key.ID)), ...
    'meanfiltVm', 'meanFR', 'filtVm', 'FR', 'metadata', 'perist', 'directions', 'intensities', 'aL', 'aR', 'trialIndicesFullBlock', 'prepost2_stimDelimitations', 'T_prepost2' )

cab(1)
%% new - check labels

figure; hold on; axis image
xlabel(' Left antenna (ipsi) - angular displ')
ylabel(' Right antenna (contra) - angular displ')
scatter(aL, aR, 120, meanfiltVm-meanfiltVm(1), 'filled', 'MarkerEdgeColor', 'k')
title(sprintf('Vm\n%s\n%s', metadata.experimentHandle.basename, metadata.key.ID), 'Interpreter', 'none')
colormap(bluewhitered(256)), cb = colorbar;
cb.TickLabels = num2str((cb.Ticks+meanfiltVm(1))', '%2.1f');
cb.Label.String = '(mV)';
export_fig(sprintf('scatter_Vm_%s_mac.pdf', metadata.key.ID))


figure; hold on; axis image
xlabel(' Left antenna (ipsi) - angular displ')
ylabel(' Right antenna (contra) - angular displ')
scatter(aL, aR, 120, meanFR-meanFR(1), 'filled', 'MarkerEdgeColor', 'k')
title(sprintf('firing rate\n%s\n%s', metadata.experimentHandle.basename, metadata.key.ID), 'Interpreter', 'none')
colormap(bluewhitered(256)), cb = colorbar;
cb.TickLabels = num2str((cb.Ticks+meanFR(1))', '%2.1f');
cb.Label.String = '(Hz)';
export_fig(sprintf('scatter_FR_%s_mac.pdf', metadata.key.ID))








    