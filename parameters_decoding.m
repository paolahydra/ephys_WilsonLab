% estimate parameters from parameter file 
% load data for a givn experiment
flyNum     = 223;
cellNum    = 1;

% [metadataname, cellfolder] = uigetfile(sprintf('D:/Dropbox (HMS)/p2/fly%3d_PP/fly%3d_cell%02d/metadata*', flyNum, flyNum, cellNum));
[metadataname, cellfolder] = uigetfile(sprintf('/Users/galileo/Dropbox (HMS)/p2/fly%3d_PP/fly%3d_cell%02d/metadata*', flyNum, flyNum, cellNum));

metadata = load(fullfile(cellfolder, metadataname));
flyfolder = metadata.experimentHandle.flyfolder;

% basename = metadata.experimentHandle.basename;
% flyNum = metadata.experimentHandle.flyNum;      % in case you picked a different one
% cellNum = metadata.experimentHandle.cellNum;    % in case you picked a different one -- IN PRINCIPLE
basename = sprintf('fly%3d_cell%02d',flyNum, cellNum);  % -- I recycled metadata here

parametersname = dir(fullfile(cellfolder, sprintf('parameters_%s_*.bin', metadata.key.ID)));
parametersname = parametersname(end).name;


HWsettings = sensor_settings_PP;

fid = fopen(fullfile(cellfolder, parametersname), 'r');
data = fread(fid, 'double' );
fclose(fid);

data = reshape(data, 7, []);        % 8 channels + 1 timestamps
timestamps = (data(1,:))';
data = data(2:end,:);
data = data';  %for consistency with TO's code

cd(cellfolder)

%% plot all voltage and current curves temporarily
[ current, voltage, ~ ] = get_scaled_voltage_and_current_PP( data );

figure; plot(timestamps, voltage); ylabel('mV'); axv = gca;
% figure; plot(timestamps, current); ylabel('pA'); axc = gca;

%% save snippets of voltage after break in
restVoltagePlotsFolder = '/Users/galileo/Dropbox (HMS)/p2/restingVoltagePlots';

XLIM = round(get(gca, 'xlim')); % choose 2 seconds by the end of param recording, after break-in.
YLIM = [-60, 10];
set(gca, 'xlim', XLIM, 'ylim', YLIM)

% export_fig(fullfile(restVoltagePlotsFolder, sprintf('%s_restingVoltage.eps', basename)) )
export_fig(fullfile(restVoltagePlotsFolder, sprintf('%s_restingVoltage_immedAfterBreakin.eps', basename)) )

%% break it up based on mode
mode_trace = data(:,HWsettings.DATAch.MODE_A_DAQ_AI);
figure; plot(timestamps, mode_trace);
I0_epoques = a; %use the first one second to determine Vm 
%shortly before I can take Raccess
cClamp_epoques = mode_trace > 0 & mode_trace < 3.5;
VClamp_epoques = mode_trace > 3.5 & mode_trace < 7;




%%
figure; plot(timestamps(1:4001), data( 1:4001, HWsettings.DATAch.VOLTAGE_NON_SCALED_A_DAQ_AI ))
hold on

[ current, voltage, scaled ] = get_scaled_voltage_and_current_PP( data );
figure; hold on
plot(timestamps(1:4001), scaled(1:4001)); ylabel('scaled, [pA]')
plot(timestamps(1:4001), current(1:4001));