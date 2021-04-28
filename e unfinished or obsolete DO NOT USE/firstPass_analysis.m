% load data for a givn experiment
flyNum     = 263;
cellNum    = 1;
% basename = 'fly233_cell03'; % temp


[metadataname, cellfolder] = uigetfile(sprintf('D:/Dropbox (HMS)/p2/fly%3d_PP/fly%3d_cell%02d/metadata*', flyNum, flyNum, cellNum));
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

cd(flyfolder)

%% convert signals and plot voltage

[ current , voltage, scaled ] = get_scaled_voltage_and_current_PP( data' );
% figure;
% plot(timestamps, scaled)


%%
XLIM = [100, 150];
indexes_x = timestamps >= XLIM(1) & timestamps <= XLIM(2);
figure;
plot(timestamps(indexes_x), voltage(indexes_x))
title(sprintf('%s - %s', basename, metadata.key.ID), 'Interpreter', 'none')
ylabel('Voltage (mV)')
xlabel('time (s)')
ylim([-60, 10])
% xlim(XLIM)
box off
set(gca, 'TickDir', 'out')
axV = gca;
savefig(sprintf('excerpt_voltageTrace_%s', metadata.key.ID))
% export_fig(sprintf('%s_%s_voltage.pdf', basename, metadata.key.ID))

% figure;
% plot(timestamps, current)
% axC = gca;
% ylabel('Current (pA)')
% xlabel('time (s)')
% ylim([-100, 100])
% xlim(XLIM)
% box off
% set(gca, 'TickDir', 'out')

% linkaxes([axV, axC], 'x')

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




figure; hold on
plot(timestamps, data(HWsettings.DATAch.PIEZO_SENSOR_RIGHT,:))
plot(timestamps(onsetsUsed), 5*ones(size(onsetsUsed)), 'xr')

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

f = figure; 
f.WindowStyle = 'normal';
f.Position = [-1117          42         186         898];
f.Color = [1,1,1];
clear periV
for i_st = 1:length(STIMS_USE)
    subplot(length(STIMS_USE),1,i_st); hold on
    stim = STIMS_USE(i_st);
    alltrials = find(trialIndicesFullBlock==stim);
    for i = 1:sum(trialIndicesFullBlock==stim)
        trialN = alltrials(i);
        onset = onsetsUsedFullBlock(trialN);
        firstI = onset - metadata.userinput.DAQ_fs*secs_baseline;
        lastI = onset + metadata.userinput.DAQ_fs*( secs_postStimOnset);
        periV(i,:) = voltage(firstI:lastI);
        if i == 1 && i_st == 1
            T = (timestamps(1:length(periV)) - secs_baseline);
        end
        plot(T, periV(i,:), '-b')
    end
    mPeriV = mean(periV(1:sum(trialIndicesFullBlock==stim),:)); %check it works - new
    plot( T, mPeriV, '-r', 'LineWidth', 2)
    ax = gca;
    ax.FontSize = 10;
    ax.FontName = 'Arial';
    ylim([-65, -30])
%     ylim([-40, 10])
    xlim([-secs_baseline, secs_postStimOnset])
    ylabel('mV')
    title(STIM_LABS{i_st})
end
xlabel('time (s)')
export_fig(sprintf('%s_%s_sampleTraces2.pdf', basename, metadata.key.ID ))
saveas(gcf,sprintf('%s_%s_sampleTraces2.eps', basename, metadata.key.ID ))


%% psth
[ cur_psth ] = calculate_psth_TO(voltage, VOLTAGE_SR, SPIKE_THRESHOLD, BIN_SIZE);


%% matrix
NCONDITS = length(metadata.dec.i2displ_decoder);
matrix = nan(NCONDITS);

% pre-build - stim[onset,offset] window
for stim = 1:NCONDITS^2
    alltrials = find(trialIndicesFullBlock==stim);
    clear periV
    for i = 1:sum(trialIndicesFullBlock==stim)
        trialN = alltrials(i);
        onset = onsetsUsedFullBlock(trialN);
        firstI = onset;
        lastI = onset + metadata.userinput.DAQ_fs*(metadata.userinput.singleStimDur);
        periV(i,:) = voltage(firstI:lastI);
    end
    mPeriV_joint(stim,:) = mean(periV);
end

% build matrix
for R = 1:NCONDITS % by rows
    for L = 1:NCONDITS % by cols
        stimN = find(metadata.dec.joint2L_decoder_i==L & metadata.dec.joint2R_decoder_i==R);
        matrix(R,L) = mean(mPeriV_joint(stimN,:));
    end
end
matrixBS = matrix - matrix(ceil(NCONDITS/2), ceil(NCONDITS/2));


figure;
imagesc(matrixBS)
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
export_fig(sprintf('%s_%s_matrix.pdf', basename, metadata.key.ID ))
%% show antennae
figure;
fr = 1;
subplot(1,3,1)
imshowpair(allFr(2).frames(:,:,fr), allFr(2).frames(:,:,2))
subplot(1,3,2)
imshowpair(allFr(1).frames(:,:,fr), allFr(1).frames(:,:,2))
subplot(1,3,3)
imshowpair(allFr(3).frames(:,:,fr), allFr(3).frames(:,:,2))

%% read in command
fidC = fopen(fullfile(cellfolder, metadata.key.logname), 'r');
command = fread(fidC, 'single');
fclose(fidC);
command = reshape(command, 2, []);


%% left piezo sensor and command
figure; hold on
DC_left = mean(data(HWsettings.DATAch.PIEZO_SENSOR_LEFT,1:20e4));
plot(timestamps, data(HWsettings.DATAch.PIEZO_SENSOR_LEFT,:))
plot(timestamps, DC_left + command(1,1:length(timestamps)))
xlim([20 60])
title(sprintf('left piezo sensor and command - %s', basename), 'Interpreter', 'none')

%% right piezo sensor and command
figure; hold on
DC_right = mean(data(HWsettings.DATAch.PIEZO_SENSOR_RIGHT,1:20e4));
plot(timestamps, data(HWsettings.DATAch.PIEZO_SENSOR_RIGHT,:))
plot(timestamps, DC_right + command(2,1:length(timestamps)))
xlim([20 60])
title(sprintf('right piezo sensor and command - %s', basename), 'Interpreter', 'none')

%%
 tic

    fprintf('\nRECORDING: elapsed time (s): %8.2f\t',toc)
    % start the loop
    i = 0;
    while i < 100       % Check if the loop has to be stopped
        pause(0.05)
        i = i+1;
        fprintf('%c',repmat(8,9,1)) ;   % clear up previous time
        fprintf('%8.2f\t',toc) ;        % display elapsed time
    end
    