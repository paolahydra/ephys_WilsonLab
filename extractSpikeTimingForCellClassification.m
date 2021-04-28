clear

flyCode = [ 1     5     7     3     4     6     8     9    10     2]; %for blinded sorting, if needed
saveDir = '/Users/galileo/Dropbox (HMS)/p2/isi_data';
isiGenData = matfile(fullfile(saveDir, 'isiGenData.mat'), 'Writable', true); %single z planes to their own stack, general infos
Nfly = length(isiGenData.FlyNum) + 1;
isiGenData.flyCode(1, Nfly) = flyCode(Nfly);
% % initialize it once:
% isiGenData.FlyNum(1,1) = 0;
% Nfly = 1;


% SET PARAMETERS
flyNum     = 279;
isiGenData.FlyNum(1, Nfly) = 279; %or use 2331 2332

cellNum    = 1;

isiGenData.fill(1, Nfly) = 1; %1 poster, 2 anterior, 0 no fill, 1.2 likely posterior

durHalfRamp_sec = 0.15;
Nfullblocks2keep = [];

SPIKEDET_STDMULTIPL = 6;
delayBetweenSpikePeaksMs = 8; %5; %this really also depends on the threshold...


% NfirstBlocksToExclude = [];   %I will not use it here, because that can only
% % be due to poor stimulus control. But we do not care about that for our
% % purpose here


[metadataname, cellfolder] = uigetfile(sprintf('/Users/galileo/Dropbox (HMS)/p2/fly%3d_PP/fly%3d_cell%02d/metadata*', flyNum, flyNum, cellNum));
metadata = load(fullfile(cellfolder, metadataname));
flyfolder = metadata.experimentHandle.flyfolder;
flyNum = metadata.experimentHandle.flyNum;      % in case you picked a different one
cellNum = metadata.experimentHandle.cellNum;    % in case you picked a different one
ephysdataname = dir(fullfile(cellfolder, sprintf('ephysData_%s_*.bin', metadata.key.ID)));
ephysdataname = ephysdataname(1).name;



HWsettings = sensor_settings_PP_olderFlies;
fid = fopen(fullfile(cellfolder, ephysdataname), 'r');
data = fread(fid, 'double' );
fclose(fid);
data = reshape(data, 10, []);       % 9 was used before 246   % use 10 after fly 246 (june 1 2018)
timestamps = (data(1,:))';
data = data(2:end,:);


% convert signals and plot voltage
[ ~ , voltage, ~ ] = get_scaled_voltage_and_current_PP( data' );

% % also include spikes shortly after break-in? NOT WORTH IT
% [paramdataname, ~] = uigetfile(sprintf('/Users/galileo/Dropbox (HMS)/p2/fly%3d_PP/fly%3d_cell%02d/parameters_%s*.bin', flyNum, flyNum, cellNum, metadata.key.ID));
% 
% fid = fopen(fullfile(cellfolder, paramdataname), 'r');
% data = fread(fid, 'double' );
% fclose(fid);
% data = reshape(data, 7, []);       % 9 was used before 246   % use 10 after fly 246 (june 1 2018)
% timestamps = (data(1,:))';
% data = data(2:end,:);
% [ ~ , voltage, ~ ] = get_scaled_voltage_and_current_PP( data' );


clearvars -except Nfly isiGenData saveDir SPIKEDET_STDMULTIPL delayBetweenSpikePeaksMs durHalfRamp_sec voltage flyfolder timestamps cellNum flyNum basename metadata metadataname cellfolder Nfullblocks2keep NfirstBlocksToExclude

%% find out last timepoint included in analysis. Anything after that may be not good signal (e.g. cell dying)
onsetsUsed = metadata.joint.allOnsetsPositions ( metadata.joint.allOnsetsPositions < length(timestamps) );
trialIndixes = metadata.joint.trialIndices(1:length(onsetsUsed));
% old flies:
% repsPerBlock = (length(metadata.dec.joint2L_decoder_i) + metadata.userinput.factor_jointZero - 1) * metadata.userinput.N_withinPseudoBlock;
% new flies:
repsPerBlock = length(metadata.dec.trialTypes) * metadata.userinput.N_withinPseudoBlock;


Nfullblocks = floor(length(onsetsUsed) / repsPerBlock)


if ~isempty(Nfullblocks2keep)
    if Nfullblocks >= Nfullblocks2keep
        Nfullblocks =  Nfullblocks2keep;
    end
end
onsetsUsedFullBlock = onsetsUsed(1: Nfullblocks*repsPerBlock);

%% crop voltage trace to the last onset, or to 15 minutes, whichever happens first
VOLTAGE_SR = metadata.userinput.DAQ_fs;
cropT_I = min(onsetsUsedFullBlock(end), 15*60*VOLTAGE_SR);
voltage = voltage(1:cropT_I);
timestamps = timestamps(1:cropT_I);

%caclulate threshold on baseline only, for consistency with previous spike
%detection
lastBaseline_I = onsetsUsedFullBlock(1) - durHalfRamp_sec*VOLTAGE_SR;
baselineStart.voltage = voltage(1:lastBaseline_I);
baselineStart.medFiltTrace = medfilt1(baselineStart.voltage, 0.08*VOLTAGE_SR, 'truncate');
baselineStart.delta_Vm = baselineStart.voltage - baselineStart.medFiltTrace; 
baselineStart.timestamps = timestamps(1:lastBaseline_I);

timestamps(end)/60

%% spike extraction
delayBetweenSpikePeaksPOINTS = delayBetweenSpikePeaksMs * VOLTAGE_SR/1e3;
SPIKE_THRESHOLD = std(baselineStart.delta_Vm)*SPIKEDET_STDMULTIPL;

figure; hold on; plot(baselineStart.timestamps, baselineStart.delta_Vm)
plot([baselineStart.timestamps(1), baselineStart.timestamps(end)], [SPIKE_THRESHOLD, SPIKE_THRESHOLD], '--r')
xlabel('time (sec)')
ylabel('mV')
title('high-pass filtered and centered voltage: starting baseline')


medFiltTrace = medfilt1(voltage, 0.08*VOLTAGE_SR, 'truncate');
delta_Vm = voltage - medFiltTrace; 

volts_thresh = delta_Vm;
volts_thresh(volts_thresh < SPIKE_THRESHOLD) = 0;
[~, locs] = findpeaks(volts_thresh, 'MinPeakDistance',delayBetweenSpikePeaksPOINTS); 
ts_spikes = timestamps(locs);

x = diff(ts_spikes);
[~,edges] = histcounts(log10(x));
figure; 
histogram(x,10.^edges)
set(gca, 'xscale','log')

%% save figure (for me) and data (for Tatsuo):
export_fig(fullfile(saveDir, sprintf('isi_hist_fly%d.pdf',isiGenData.FlyNum(1, Nfly))))
isiGenData.spike_ts(1, Nfly) = {ts_spikes};

%% combine two runs if too short?
ts_spikes1 = isiGenData.spike_ts(1, Nfly);
ts_spikes1 = ts_spikes1{1};
ts_spikes_zeroed = ts_spikes -ts_spikes(1); 
ts_spikes2 = cat(1, ts_spikes1, ts_spikes1(end) + mode(ts_spikes1) + ts_spikes_zeroed);

isiGenData.spike_ts(1, Nfly) = {ts_spikes2};


x = diff(ts_spikes2);
[~,edges] = histcounts(log10(x));
figure; 
histogram(x,10.^edges)
set(gca, 'xscale','log')


